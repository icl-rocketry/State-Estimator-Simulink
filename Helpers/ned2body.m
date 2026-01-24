function v_b = ned2body(v_n, quat)
%NED2BODY Rotate a vector from NED frame to body frame using quaternion
%
% Inputs:
%   v_n   : [3x1] vector in NED or world frame
%   quat  : [1x4] quaternion [w x y z] (body wrt NED)
% Output:
%   v_b   : [3x1] vector in body frame

arguments
    v_n (3,1) double
    quat (1,4) double
end

R_bn = quat2rotm(quat);    % Aerospace Toolbox
v_b  = R_bn * v_n;
end
