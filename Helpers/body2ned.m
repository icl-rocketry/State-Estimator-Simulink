function v_n = body2ned(v_b, quat)
%BODY2NED Rotate vector from body frame to NED frame
%
% quat: [w x y z], representing NED → body rotation

arguments
    v_b (3,1) double
    quat (4,1) double
end

R_bn = quat2rotm(quat');   % body → NED

v_n = R_bn * v_b;

end