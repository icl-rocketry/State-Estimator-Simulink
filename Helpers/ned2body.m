function v_b = ned2body(v_n, quat)
%NED2BODY Rotate vector from NED into body frame.
%
% quat: [w;x;y;z], representing BODY -> NED rotation.
% quat must be normalised.

arguments
    v_n (3,1) double
    quat (4,1) double
end

R_body_to_ned = quat2rotm(quat.');

v_b = R_body_to_ned.' * v_n;

end