%% ============================================================
%  3D Synthetic Sensor Data for KF Testing (NED frame)
%% ============================================================
close all;
clear; clc;

% -------------------------------
% Time
% --------------------------------
model_update_rate = 500; % Hz
GPS_update_rate   = 10;  % Hz
Baro_update_rate  = 10; % Hz
Accel_update_rate = 10; % Hz

dt = 1 / model_update_rate;        % 1/f
GPS_dt = 1 / GPS_update_rate;   % new GPS info/ x runs of model update
Baro_dt = 1 / Baro_update_rate;   % new GPS info/ x runs of model update
Accel_dt = 1 / Accel_update_rate;   % new GPS info/ x runs of model update

T  = 60;          % seconds
t  = (0:dt:T)';
N  = length(t);

% -------------------------------
% Reference LLA
% --------------------------------
lat0 = 51.507351;     % deg
lon0 = -0.127758;    % deg
h0   = 25;            % m

% -------------------------------
% TRUE ACCELERATION (NED)
% --------------------------------
a_true = zeros(N,3);

% North accel pulse
a_true(t > 5 & t < 10, 1) = 0.3;

% East accel pulse
a_true(t > 20 & t < 25, 2) = -0.2;

% Down accel pulse
a_true(t > 35 & t < 40, 3) = 0.5;

% -------------------------------
% TRUE VELOCITY (NED)
% --------------------------------
v_true = cumtrapz(t, a_true);

% Add constant cruise velocity
v_true(:,1) = v_true(:,1);    % 5 m/s North
v_true(:,2) = v_true(:,2);    % 2 m/s East

% -------------------------------
% TRUE POSITION (NED)
% --------------------------------
p_true = cumtrapz(t, v_true);   % [x_N, y_E, z_D]

% ============================================================
% SENSOR MODELS
% ============================================================

% Indices for each sensor rate (assumes exact integer ratio)
k_gps   = round(GPS_dt   / dt);     % 0.1 / 0.002 = 50
k_baro  = round(Baro_dt  / dt);     % 0.004 / 0.002 = 2
k_accel = round(Accel_dt / dt);     % 0.004 / 0.002 = 2

idx_gps   = 1:k_gps:N;
idx_baro  = 1:k_baro:N;
idx_accel = 1:k_accel:N;

t_gps   = t(idx_gps);
t_baro  = t(idx_baro);
t_accel = t(idx_accel);

assert(abs(k_gps*dt - GPS_dt) < 1e-12,  'GPS rate not integer multiple of base dt');
assert(abs(k_baro*dt - Baro_dt) < 1e-12,'Baro rate not integer multiple of base dt');
assert(abs(k_accel*dt - Accel_dt) < 1e-12,'Accel rate not integer multiple of base dt');

% -------------------------------
% Accelerometer (Low-G)
% --------------------------------
sigma_accel = 0.05; % m/s^2

a_true_accel   = a_true(idx_accel, :);                          % Mx3
accel_meas     = a_true_accel + sigma_accel*randn(numel(idx_accel),3);

accel_ts = timeseries(accel_meas, t_accel);
accel_ts.Name = 'accel_NED';

% -------------------------------
% Barometer (Down position only)
% --------------------------------
sigma_baro = 1.0; % m

pz_true_baro = p_true(idx_baro, 3);                             % Mx1 (Down)
baro_meas    = pz_true_baro + sigma_baro*randn(numel(idx_baro),1);

baro_ts = timeseries(baro_meas, t_baro);
baro_ts.Name = 'baro_pz';

% -------------------------------
% GPS (LLA + NED velocity)
% --------------------------------
sigma_gps_pos = 1.5; % m
sigma_gps_vel = 0.1; % m/s

% Noisy NED position and velocity sampled at GPS instants
p_gps = p_true(idx_gps,:) + sigma_gps_pos*randn(numel(idx_gps),3);
v_gps = v_true(idx_gps,:) + sigma_gps_vel*randn(numel(idx_gps),3);

% Convert to LLA at GPS instants
lat_gps = lat0 + p_gps(:,1) / 111320;
lon_gps = lon0 + p_gps(:,2) ./ (111320 * cosd(lat0));
h_gps   = h0  - p_gps(:,3);

gps_meas = [lat_gps, lon_gps, h_gps, v_gps(:,1), v_gps(:,2), v_gps(:,3)]; % Mx6

gps_ts = timeseries(gps_meas, t_gps);
gps_ts.Name = 'gps_LLA_NEDvel';

fprintf("Base samples:  %d at %.1f Hz\n", N, model_update_rate);
fprintf("Accel samples: %d at %.1f Hz\n", accel_ts.Length, Accel_update_rate);
fprintf("Baro samples:  %d at %.1f Hz\n", baro_ts.Length,  Baro_update_rate);
fprintf("GPS samples:   %d at %.1f Hz\n", gps_ts.Length,   GPS_update_rate);

% ============================================================
%  RUN SIMULINK + PLOT RESULTS
% ============================================================

% --- Make sure model is on path / in current folder
modelName = "Kalman_tester";

% Optional: load sensor data into base workspace (From Workspace reads here)
assignin('base', 'accel_ts', accel_ts);
assignin('base', 'baro_ts',  baro_ts);
assignin('base', 'gps_ts',   gps_ts);
assignin('base', 'dt',       dt);   % if your model uses dt block from workspace

% --- Load model
load_system(modelName);

% --- Configure stop time to match our generated data
set_param(modelName, 'StopTime', num2str(T));

% (Recommended) set solver to discrete if you want strict dt stepping
% Uncomment if needed:
% set_param(modelName, 'SolverType', 'Fixed-step');
% set_param(modelName, 'Solver', 'discrete');
% set_param(modelName, 'FixedStep', num2str(dt));

% --- Run simulation
simOut = sim(modelName);

%%
%% -------------------------------
% Extract measured/estimated state
% -------------------------------

ts = simOut.measured_x;          % timeseries

t_est = ts.Time;                 % 6001x1
x_est = squeeze(ts.Data).';      % (9x6001) -> (6001x9)

% Now:
% x_est(k,:) = [pN pE pD vN vE vD aN aE aD] at time t_est(k)


%% -------------------------------
% Interpolate true signals to estimator timebase (for clean overlay)
% -------------------------------
p_true_i = interp1(t, p_true, t_est, 'linear', 'extrap');
v_true_i = interp1(t, v_true, t_est, 'linear', 'extrap');
a_true_i = interp1(t, a_true, t_est, 'linear', 'extrap');

% Assume state order: [pN pE pD vN vE vD aN aE aD]
p_est = x_est(:,1:3);
v_est = x_est(:,4:6);
a_est = x_est(:,7:9);

%% -------------------------------
% Plot: Position
% -------------------------------
figure(1);
plot(t_est, p_true_i(:,1), t_est, p_est(:,1)); grid on;
xlabel('Time (s)'); ylabel('North position (m)');
legend('True','Estimated');

figure(2);
plot(t_est, p_true_i(:,2), t_est, p_est(:,2)); grid on;
xlabel('Time (s)'); ylabel('East position (m)');
legend('True','Estimated');

figure(3);
plot(t_est, p_true_i(:,3), t_est, p_est(:,3)); grid on;
xlabel('Time (s)'); ylabel('Down position (m)');
legend('True','Estimated');

%% -------------------------------
% Plot: Velocity
% -------------------------------
figure(4);
plot(t_est, v_true_i(:,1), t_est, v_est(:,1)); grid on;
xlabel('Time (s)'); ylabel('v_N (m/s)');
legend('True','Estimated');

figure(5);
plot(t_est, v_true_i(:,2), t_est, v_est(:,2)); grid on;
xlabel('Time (s)'); ylabel('v_E (m/s)');
legend('True','Estimated');

figure(6);
plot(t_est, v_true_i(:,3), t_est, v_est(:,3)); grid on;
xlabel('Time (s)'); ylabel('v_D (m/s)');
legend('True','Estimated');

%% -------------------------------
% Plot: Acceleration
% -------------------------------
figure(7);
plot(t_est, a_true_i(:,1), t_est, a_est(:,1)); grid on;
xlabel('Time (s)'); ylabel('a_N (m/s^2)');
legend('True','Estimated');

figure(8);
plot(t_est, a_true_i(:,2), t_est, a_est(:,2)); grid on;
xlabel('Time (s)'); ylabel('a_E (m/s^2)');
legend('True','Estimated');

figure(9);
plot(t_est, a_true_i(:,3), t_est, a_est(:,3)); grid on;
xlabel('Time (s)'); ylabel('a_D (m/s^2)');
legend('True','Estimated');

fprintf("Simulation complete. Plots generated.\n");