function v_n = body2ned(v_b, quat)
%BODY2NED Rotate vector from body frame into NED.
%
% quat: [w;x;y;z], representing BODY -> NED rotation.
% quat must be normalised.

arguments
    v_b (3,1) double
    quat (4,1) double
end

R_body_to_ned = quat2rotm(quat.');

v_n = R_body_to_ned * v_b;

end