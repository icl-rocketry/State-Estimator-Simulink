function [s_r, sdot_r, sddot_r] = trapezoidalRefBlockXYZ(t, T)
%#codegen
% 3-D trapezoidal reference with smooth polynomial blending between
% acceleration, cruise and deceleration phases – and between consecutive
% motion stages (0→via and via→final).
%

%% 1. Shorthand aliases ----------------------------------------------------
x0=T.x0; xi=T.xi; xf=T.xf;
y0=T.y0; yi=T.yi; yf=T.yf;
z0=T.z0; zi=T.zi; zf=T.zf;
%% ---- 2. Pre-allocate ----------------------------------------------------
s_r     = zeros(3,1);
sdot_r  = zeros(3,1);
sddot_r = zeros(3,1);
% -------------------------------------------------------------------------
% Configure the key-points of the motion here:
% -------------------------------------------------------------------------
%            [x;   y;   z]
 p0 = [x0;    y0;   z0];   % start position
 pi = [xi;    yi;   zi];   % via-point (intermediate target)
 pf = [xf;    yf;   zf];   % final position

 t1=T.t1;  t2=T.t2;
 dt_blend = T.dt_blend;          % polynomial blend window at each transition [s]
 tblend_frac = T.frac_tb;
% -------------------------------------------------------------------------

% Pre-allocate outputs (x,y,z)
 s_r     = zeros(3,1);
 sdot_r  = zeros(3,1);
 sddot_r = zeros(3,1);
% 0b) Step-mode override (any axis) -------------------------------
if T.useStep
    % use initial start pos from your traj struct:
    p0 = [T.x0; T.y0; T.z0];
    % if time < stepTime, stay at p0, else jump to stepTarget
    if t < T.stepTime
        s_r = p0;
    else
        s_r = T.stepTarget;
    end
    sdot_r  = zeros(3,1);
    sddot_r = zeros(3,1);
    return;
end
%=================== 1) Final stop blend (after second stage) =============
 if t >= (t1 + t2 - dt_blend/2)
     t_end = t1 + t2;
     alpha = min((t - (t_end - dt_blend/2)) / dt_blend, 1);
     sigma = 6*alpha^5 - 15*alpha^4 + 10*alpha^3;   % smooth-step 5-th order

     tau2 = t - t1;                                         % local time in stage-2
     [s_pre, sdot_pre, sddot_pre] = localTrapezoidalStage(tau2, pi, pf, t2, dt_blend);

     s_post     = pf;                       % stationary target
     sdot_post  = zeros(3,1);
     sddot_post = zeros(3,1);

     s_r     = (1 - sigma)*s_pre     + sigma*s_post;
     sdot_r  = (1 - sigma)*sdot_pre  + sigma*sdot_post;
     sddot_r = (1 - sigma)*sddot_pre + sigma*sddot_post;
     return;
 end

%=================== 2) Blend between stage-1 and stage-2 =================
 if t > t1 - dt_blend/2 && t < t1 + dt_blend/2
     alpha = (t - (t1 - dt_blend/2)) / dt_blend;
     sigma = 6*alpha^5 - 15*alpha^4 + 10*alpha^3;

     tau1 = t;                 % local time in stage-1
     tau2 = t - t1;            % local time in stage-2

     [s1, sdot1, sddot1] = localTrapezoidalStage(tau1, p0, pi, t1, dt_blend);
     [s2, sdot2, sddot2] = localTrapezoidalStage(tau2, pi, pf, t2, dt_blend);

     s_r     = (1 - sigma)*s1     + sigma*s2;
     sdot_r  = (1 - sigma)*sdot1  + sigma*sdot2;
     sddot_r = (1 - sigma)*sddot1 + sigma*sddot2;
     return;
 end

%=================== 3) Regular evaluation inside one stage ===============
 if t <= t1
     tau = t;   pA = p0;  pB = pi;  tf = t1;
 else
     tau = t - t1;  pA = pi;  pB = pf;  tf = t2;
 end
 [s_r, sdot_r, sddot_r] = localTrapezoidalStage(tau, pA, pB, tf, dt_blend);
end

%=========================================================================
function [s_r, sdot_r, sddot_r] = localTrapezoidalStage(tau, pA, pB, tf, dt_blend) %#codegen
% Vector-form trapezoidal motion with quintic-polynomial blends at the
% cruise boundaries.  Valid for column vectors pA, pB (size = 3×1).

 tb =0.2 * tf;                     % acceleration / deceleration time
 if tb > tf/2;  tb = tf/2;  end
 tc = tf - 2*tb;                    % constant-velocity slice duration

 D   = pB - pA;                     % total displacement (vector)
 den = tb^2 + tb*tc;                % shared denominator
 a   = D / den;                     % constant acceleration (vector)
 V   = a * tb;                      % cruise velocity (vector)

 s_r     = zeros(3,1);
 sdot_r  = zeros(3,1);
 sddot_r = zeros(3,1);

 % ---- 3.1  Acceleration phase ------------------------------------------
 if tau <= tb - dt_blend/2
     sddot_r =  a;
     sdot_r  =  a * tau;
     s_r     =  pA + 0.5*a * tau^2;

 % ---- 3.2  Blend  (acc → cruise) ---------------------------------------
 elseif tau <= tb + dt_blend/2
     tau_b  = tau - (tb - dt_blend/2);
     alpha  = tau_b / dt_blend;
     sigma  = 6*alpha^5 - 15*alpha^4 + 10*alpha^3;

     % edge-A (still accelerating)
     sddot_a = a;
     sdot_a  = a * tau;
     s_a     = pA + 0.5*a * tau^2;

     % edge-C (already at cruise)
     sddot_c = zeros(3,1);
     sdot_c  = V;
     s_c     = pA + 0.5*a * tb^2 + V * (tau - tb);

     s_r     = (1 - sigma).*s_a     + sigma.*s_c;
     sdot_r  = (1 - sigma).*sdot_a  + sigma.*sdot_c;
     sddot_r = (1 - sigma).*sddot_a + sigma.*sddot_c;

 % ---- 3.3  Constant-velocity phase -------------------------------------
 elseif tau <= tb + tc - dt_blend/2
     sddot_r = zeros(3,1);
     sdot_r  = V;
     s_r     = pA + 0.5*a * tb^2 + V * (tau - tb);

 % ---- 3.4  Blend  (cruise → dec) ---------------------------------------
 elseif tau <= tb + tc + dt_blend/2
     tau_b  = tau - (tb + tc - dt_blend/2);
     alpha  = tau_b / dt_blend;
     sigma  = 6*alpha^5 - 15*alpha^4 + 10*alpha^3;

     % edge-C (still cruising)
     sddot_c = zeros(3,1);
     sdot_c  = V;
     s_c     = pA + 0.5*a * tb^2 + V * (tau - tb);

     % edge-D (already decelerating)
     tau2 = tau - (tb + tc);
     p3   = pA + 0.5*a * tb^2 + V * tc;
     sddot_d = -a;
     sdot_d  = V - a * tau2;
     s_d     = p3 + V * tau2 - 0.5*a * tau2.^2;

     s_r     = (1 - sigma).*s_c     + sigma.*s_d;
     sdot_r  = (1 - sigma).*sdot_c  + sigma.*sdot_d;
     sddot_r = (1 - sigma).*sddot_c + sigma.*sddot_d;

 % ---- 3.5  Deceleration phase ------------------------------------------
 else
     tau2 = tau - (tb + tc);
     p3   = pA + 0.5*a * tb^2 + V * tc;
     sddot_r = -a;
     sdot_r  = V - a * tau2;
     s_r     = p3 + V * tau2 - 0.5*a * tau2.^2;
 end
end