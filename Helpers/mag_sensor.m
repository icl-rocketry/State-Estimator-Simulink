function mag_output = mag_sensor(orient_pos)
% orient_pos = [quaternion(4); position_ned(3)]
% mag_output = [field_body_mG; field_ned_mG]

quaternion = orient_pos(1:4);
position_down = orient_pos(7);

% Location used by the existing model
latitude_deg = 51.49897;
longitude_deg = -0.17430;

% Assumes the NED origin has zero ellipsoidal height
home_ellipsoid_height = 0;
ellipsoid_height = home_ellipsoid_height - position_down;

decimal_year = decyear(2025,2,24);

% WMM output is [North;East;Down] in nanotesla
[field_ned,~,~,~,~] = wrldmagm( ...
    ellipsoid_height, ...
    latitude_deg, ...
    longitude_deg, ...
    decimal_year, ...
    '2025');

field_ned = field_ned(:);

% Convert nT to mG
field_ned_mG = field_ned ./ 100;

% Normalize quaternion
quaternion = quaternion(:);
q_norm = norm(quaternion);

if q_norm < 1e-12
    quaternion = [1;0;0;0];
else
    quaternion = quaternion ./ q_norm;
end

% Rotate NED magnetic field into body axes
field_body_mG = ned2body(field_ned_mG, quaternion);

mag_output = [field_body_mG;
              field_ned_mG];

end