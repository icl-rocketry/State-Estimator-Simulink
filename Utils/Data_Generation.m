%% ============================================================
%  3D Synthetic Sensor Data for KF Testing (NED frame)
%% ============================================================
close all;
clear; clc;

% -------------------------------
% Time
% --------------------------------
model_update_rate = 500; % Hz
GPS_update_rate   = 500;  % Hz
Baro_update_rate  = 500; % Hz
Accel_update_rate = 500; % Hz
Gyro_update_rate = 500; % Hz

dt = 1 / model_update_rate;        % 1/f
GPS_dt = 1 / GPS_update_rate;
Baro_dt = 1 / Baro_update_rate;
Accel_dt = 1 / Accel_update_rate;
Gyro_dt = 1 / Gyro_update_rate;

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
% TRUE ATTITUDE + TRUE GYRO
% Robust against no-motion / empty-find / NaN propagation
% ============================================================

eps_a = 1e-4;
tau_att = 0.25;
alpha = dt / (tau_att + dt);

% Ensure column time vector and consistent length
t = t(:);
N = length(t);

% Magnitude of acceleration
a_mag = vecnorm(a_true, 2, 2);

% Unit acceleration direction
u = zeros(N,3);
idx_move = isfinite(a_mag) & (a_mag > eps_a);

u(idx_move,:) = a_true(idx_move,:) ./ a_mag(idx_move);

% Desired yaw/pitch
yaw_des   = NaN(N,1);
pitch_des = NaN(N,1);

yaw_des(idx_move)   = atan2(u(idx_move,2), u(idx_move,1));

u3 = u(:,3);
u3 = max(-1, min(1, u3));    % clamp
pitch_des(idx_move) = -asin(u3(idx_move));

% If no motion ever occurs, default to zero attitude
if any(idx_move)
    first_idx = find(idx_move,1,'first');

    yaw_des   = fillmissing(yaw_des,   'previous');
    pitch_des = fillmissing(pitch_des, 'previous');

    yaw_des(1:first_idx-1)   = yaw_des(first_idx);
    pitch_des(1:first_idx-1) = pitch_des(first_idx);
else
    yaw_des(:) = 0;
    pitch_des(:) = 0;
end

% Unwrap yaw
yaw_des = unwrap(yaw_des);

% Roll excitation
roll_amp_deg = 20;
roll_freq_hz = 0.15;
roll_des = deg2rad(roll_amp_deg) * sin(2*pi*roll_freq_hz*t);

% Smooth
roll  = zeros(N,1);
pitch = zeros(N,1);
yaw   = zeros(N,1);

roll(1)  = roll_des(1);
pitch(1) = pitch_des(1);
yaw(1)   = yaw_des(1);

for k = 2:N
    roll(k)  = roll(k-1)  + alpha*(roll_des(k)  - roll(k-1));
    pitch(k) = pitch(k-1) + alpha*(pitch_des(k) - pitch(k-1));
    yaw(k)   = yaw(k-1)   + alpha*(yaw_des(k)   - yaw(k-1));
end

% Derivatives
rolldot  = [diff(roll)/dt; 0];
pitchdot = [diff(pitch)/dt; 0];
yawdot   = [diff(yaw)/dt; 0];

% Full Euler-rate to body-rate mapping
wx_true = rolldot - sin(pitch).*yawdot;
wy_true = cos(roll).*pitchdot + sin(roll).*cos(pitch).*yawdot;
wz_true = -sin(roll).*pitchdot + cos(roll).*cos(pitch).*yawdot;

omega_true = [wx_true, wy_true, wz_true];

% Final cleanup guard
omega_true(~isfinite(omega_true)) = 0;
%% -------------------------------
% Convert True Euler Angles -> Quaternion
% scalar-first convention: [q0 q1 q2 q3]
% -------------------------------
q_true = zeros(N,4);

for k = 1:N

    cr = cos(roll(k)/2);
    sr = sin(roll(k)/2);

    cp = cos(pitch(k)/2);
    sp = sin(pitch(k)/2);

    cy = cos(yaw(k)/2);
    sy = sin(yaw(k)/2);

    q0 = cr*cp*cy + sr*sp*sy;
    q1 = sr*cp*cy - cr*sp*sy;
    q2 = cr*sp*cy + sr*cp*sy;
    q3 = cr*cp*sy - sr*sp*cy;

    q_true(k,:) = [q0 q1 q2 q3];
end

% Optional but good practice: normalize
q_true = q_true ./ vecnorm(q_true,2,2);
% ============================================================
% SENSOR MODELS
% ============================================================

% Indices for each sensor rate (assumes exact integer ratio)
k_gps   = round(GPS_dt   / dt);
k_baro  = round(Baro_dt  / dt);
k_accel = round(Accel_dt / dt);
k_gyro  = round(Gyro_dt  / dt);

idx_gps   = 1:k_gps:N;
idx_baro  = 1:k_baro:N;
idx_accel = 1:k_accel:N;
idx_gyro  = 1:k_gyro:N;

t_gps   = t(idx_gps);
t_baro  = t(idx_baro);
t_accel = t(idx_accel);
t_gyro  = t(idx_gyro);

assert(abs(k_gps*dt - GPS_dt) < 1e-12,  'GPS rate not integer multiple of base dt');
assert(abs(k_baro*dt - Baro_dt) < 1e-12,'Baro rate not integer multiple of base dt');
assert(abs(k_accel*dt - Accel_dt) < 1e-12,'Accel rate not integer multiple of base dt');
assert(abs(k_gyro*dt - Gyro_dt) < 1e-12,'Gyro rate not integer multiple of base dt');

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

% -------------------------------
% Gyroscope (sample true body rates + add noise)
% --------------------------------
sigma_gyro = 0.005; % rad/s (tune)

omega_true_gyro = omega_true(idx_gyro,:);   % sample at gyro instants (here = accel instants)
gyro_meas = omega_true_gyro + sigma_gyro*randn(numel(idx_gyro),3);

gyro_ts = timeseries(gyro_meas, t_gyro);
gyro_ts.Name = 'gyro_pqr_radps';

% ============================================================
%  DISPLAY SENSOR INFO
% ============================================================
fprintf("Base samples:  %d at %.1f Hz\n", N, model_update_rate);
fprintf("Accel samples: %d at %.1f Hz\n", accel_ts.Length, Accel_update_rate);
fprintf("Baro samples:  %d at %.1f Hz\n", baro_ts.Length,  Baro_update_rate);
fprintf("GPS samples:   %d at %.1f Hz\n", gps_ts.Length,   GPS_update_rate);
fprintf("Gyro samples:  %d at %.1f Hz\n", gyro_ts.Length,  Gyro_update_rate);

% ============================================================
%  RUN SIMULINK + PLOT RESULTS
% ============================================================

% --- Make sure model is on path / in current folder
modelName = "Kalman_tester";

% Optional: load sensor data into base workspace (From Workspace reads here)
assignin('base', 'gyro_ts', gyro_ts);
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

%% -------------------------------
% Extract measured/estimated state
% -------------------------------

ts = simOut.measured_x;          % timeseries

t_est = ts.Time;
x_est = squeeze(ts.Data).';      

% 12-state order:
% x_est(k,:) = [pN pE pD vN vE vD aN aE aD roll pitch yaw]

%% -------------------------------
% Interpolate true signals to estimator timebase (for clean overlay)
% -------------------------------
p_true_i = interp1(t, p_true, t_est, 'linear', 'extrap');
v_true_i = interp1(t, v_true, t_est, 'linear', 'extrap');
a_true_i = interp1(t, a_true, t_est, 'linear', 'extrap');

p_est = x_est(:,1:3);
v_est = x_est(:,4:6);
a_est = x_est(:,7:9);
quart_est = x_est(:,10:13);

q0_est  = quart_est(:,1);
q1_est = quart_est(:,2);
q2_est   = quart_est(:,3);
q3_est   = quart_est(:,4);
%% -------------------------------
% Plot formatting defaults
% -------------------------------
lw = 1.8;
fs = 12;

%% -------------------------------
% Plot: Position
% -------------------------------
figure(1); clf;
plot(t_est, p_true_i(:,1), 'b', 'LineWidth', lw); hold on;
plot(t_est, p_est(:,1), 'r--', 'LineWidth', lw);
grid on; box on;
xlabel('Time (s)', 'FontSize', fs);
ylabel('North Position (m)', 'FontSize', fs);
title('North Position: True vs Estimated', 'FontSize', fs+1);
legend('True','Estimated','Location','best');
set(gca,'FontSize',fs);

figure(2); clf;
plot(t_est, p_true_i(:,2), 'b', 'LineWidth', lw); hold on;
plot(t_est, p_est(:,2), 'r--', 'LineWidth', lw);
grid on; box on;
xlabel('Time (s)', 'FontSize', fs);
ylabel('East Position (m)', 'FontSize', fs);
title('East Position: True vs Estimated', 'FontSize', fs+1);
legend('True','Estimated','Location','best');
set(gca,'FontSize',fs);

figure(3); clf;
plot(t_est, p_true_i(:,3), 'b', 'LineWidth', lw); hold on;
plot(t_est, p_est(:,3), 'r--', 'LineWidth', lw);
grid on; box on;
xlabel('Time (s)', 'FontSize', fs);
ylabel('Down Position (m)', 'FontSize', fs);
title('Down Position: True vs Estimated', 'FontSize', fs+1);
legend('True','Estimated','Location','best');
set(gca,'FontSize',fs);

%% -------------------------------
% Plot: Velocity
% -------------------------------
figure(4); clf;
plot(t_est, v_true_i(:,1), 'b', 'LineWidth', lw); hold on;
plot(t_est, v_est(:,1), 'r--', 'LineWidth', lw);
grid on; box on;
xlabel('Time (s)', 'FontSize', fs);
ylabel('v_N (m/s)', 'FontSize', fs);
title('North Velocity: True vs Estimated', 'FontSize', fs+1);
legend('True','Estimated','Location','best');
set(gca,'FontSize',fs);

figure(5); clf;
plot(t_est, v_true_i(:,2), 'b', 'LineWidth', lw); hold on;
plot(t_est, v_est(:,2), 'r--', 'LineWidth', lw);
grid on; box on;
xlabel('Time (s)', 'FontSize', fs);
ylabel('v_E (m/s)', 'FontSize', fs);
title('East Velocity: True vs Estimated', 'FontSize', fs+1);
legend('True','Estimated','Location','best');
set(gca,'FontSize',fs);

figure(6); clf;
plot(t_est, v_true_i(:,3), 'b', 'LineWidth', lw); hold on;
plot(t_est, v_est(:,3), 'r--', 'LineWidth', lw);
grid on; box on;
xlabel('Time (s)', 'FontSize', fs);
ylabel('v_D (m/s)', 'FontSize', fs);
title('Down Velocity: True vs Estimated', 'FontSize', fs+1);
legend('True','Estimated','Location','best');
set(gca,'FontSize',fs);

%% -------------------------------
% Plot: Acceleration
% -------------------------------
figure(7); clf;
plot(t_est, a_true_i(:,1), 'b', 'LineWidth', lw); hold on;
plot(t_est, a_est(:,1), 'r--', 'LineWidth', lw);
grid on; box on;
xlabel('Time (s)', 'FontSize', fs);
ylabel('a_N (m/s^2)', 'FontSize', fs);
title('North Acceleration: True vs Estimated', 'FontSize', fs+1);
legend('True','Estimated','Location','best');
set(gca,'FontSize',fs);

figure(8); clf;
plot(t_est, a_true_i(:,2), 'b', 'LineWidth', lw); hold on;
plot(t_est, a_est(:,2), 'r--', 'LineWidth', lw);
grid on; box on;
xlabel('Time (s)', 'FontSize', fs);
ylabel('a_E (m/s^2)', 'FontSize', fs);
title('East Acceleration: True vs Estimated', 'FontSize', fs+1);
legend('True','Estimated','Location','best');
set(gca,'FontSize',fs);

figure(9); clf;
plot(t_est, a_true_i(:,3), 'b', 'LineWidth', lw); hold on;
plot(t_est, a_est(:,3), 'r--', 'LineWidth', lw);
grid on; box on;
xlabel('Time (s)', 'FontSize', fs);
ylabel('a_D (m/s^2)', 'FontSize', fs);
title('Down Acceleration: True vs Estimated', 'FontSize', fs+1);
legend('True','Estimated','Location','best');
set(gca,'FontSize',fs);

%% -------------------------------
% Plot: Quaternions
% -------------------------------

% True quaternion components from your model
% Adjust these names to match your workspace variables
q0_true = q_true(:,1);
q1_true = q_true(:,2);
q2_true = q_true(:,3);
q3_true = q_true(:,4);

% Interpolate true quaternion onto estimator time base
q0_true_i = interp1(t, q0_true, t_est, 'linear', 'extrap');
q1_true_i = interp1(t, q1_true, t_est, 'linear', 'extrap');
q2_true_i = interp1(t, q2_true, t_est, 'linear', 'extrap');
q3_true_i = interp1(t, q3_true, t_est, 'linear', 'extrap');

figure(10); clf;
plot(t_est, q0_true_i, 'b', 'LineWidth', lw); hold on;
plot(t_est, q0_est, 'r--', 'LineWidth', lw);
grid on; box on;
xlabel('Time (s)', 'FontSize', fs);
ylabel('q_0', 'FontSize', fs);
title('Quaternion q_0: True vs Estimated', 'FontSize', fs+1);
legend('True','Estimated','Location','best');
set(gca,'FontSize',fs);

figure(11); clf;
plot(t_est, q1_true_i, 'b', 'LineWidth', lw); hold on;
plot(t_est, q1_est, 'r--', 'LineWidth', lw);
grid on; box on;
xlabel('Time (s)', 'FontSize', fs);
ylabel('q_1', 'FontSize', fs);
title('Quaternion q_1: True vs Estimated', 'FontSize', fs+1);
legend('True','Estimated','Location','best');
set(gca,'FontSize',fs);

figure(12); clf;
plot(t_est, q2_true_i, 'b', 'LineWidth', lw); hold on;
plot(t_est, q2_est, 'r--', 'LineWidth', lw);
grid on; box on;
xlabel('Time (s)', 'FontSize', fs);
ylabel('q_2', 'FontSize', fs);
title('Quaternion q_2: True vs Estimated', 'FontSize', fs+1);
legend('True','Estimated','Location','best');
set(gca,'FontSize',fs);

figure(13); clf;
plot(t_est, q3_true_i, 'b', 'LineWidth', lw); hold on;
plot(t_est, q3_est, 'r--', 'LineWidth', lw);
grid on; box on;
xlabel('Time (s)', 'FontSize', fs);
ylabel('q_3', 'FontSize', fs);
title('Quaternion q_3: True vs Estimated', 'FontSize', fs+1);
legend('True','Estimated','Location','best');
set(gca,'FontSize',fs);

fprintf("Simulation complete. Plots generated.\n");