function true_lla = gps_converter(position)
% position = [N; E; D] in meters (NED frame, D positive down)
% true_lla = [lat; lon; h]

north = position(1);
east  = position(2);
down  = position(3);

lat0 = 51.507351;
lon0 = -0.127758;
h0   = 25;

spheroid = wgs84Ellipsoid("meter");
[lat, lon, h] = ned2geodetic(north, east, down, lat0, lon0, h0, spheroid);

true_lla = [lat; lon; h];
end