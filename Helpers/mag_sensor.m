function B_b = mag_sensor(orient_pos)

% Magnetic field setup (runs before simulation)

height = orient_pos(7);  % or from initial condition
lat = 51.49897;
lon = -0.17430;

time = decyear(2025,2,24);   %

[B_w, ~, ~, ~, ~] = wrldmagm(height, lat, lon, time, '2025');

% Ensure column vectors
Orientation = orient_pos(1:4);
Orientation = Orientation(:);
B_w = B_w(:);

% Convert from nT to mG
B_w = B_w ./ 100;

% Normalize quaternion
q = Orientation / norm(Orientation);

% Quaternion inverse (for unit quaternion)
q_inv = [q(1); -q(2:4)];

% Convert vector to pure quaternion
B_quat = [0; B_w];

% Perform rotation: B_b = q^{-1} * B_w * q, converting Earth magnetic field in NED frame to body frame
B_rot = quatmultiply(quatmultiply(q_inv.', B_quat.'), q.'); 

% Extract vector part
B_b = [B_rot(2:4).'; B_w];

end