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

qw = quat(1); qx = quat(2); qy = quat(3); qz = quat(4);

R_bn = [1 - 2*(qy^2 + qz^2),     2*(qx*qy - qz*qw),     2*(qx*qz + qy*qw);
        2*(qx*qy + qz*qw), 1 - 2*(qx^2 + qz^2),     2*(qy*qz - qx*qw);
        2*(qx*qz - qy*qw),     2*(qy*qz + qx*qw), 1 - 2*(qx^2 + qy^2)];

v_b  = R_bn' * v_n;
end
