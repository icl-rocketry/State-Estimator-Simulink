function v_b = ned2body(v_n, quat)
%NED2BODY Rotate a vector from NED frame to body frame using quaternion
%
% Inputs:
%   v_n   : [3x1] vector in NED or world frame
%   quat  : [4x1] quaternion vector [w x y z] (NED to body)
% Output:
%   v_b   : [3x1] vector in body frame

arguments
    v_n (3,1) double
    quat (4,1) double
end

R_nb = quat2rotm(quat');    % Aerospace Toolbox
v_b  = R_nb * v_n;
end
