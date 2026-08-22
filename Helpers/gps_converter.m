function true_lla = gps_converter(position)
% position = [N; E; D] in meters (NED frame, D positive down)
% true_lla = [lat; lon; h]

north = position(1);
east  = position(2);
down  = position(3);

lat0 = 51.49897;
lon0 = -0.17430;
h0   = 0;

spheroid = wgs84Ellipsoid("meter");
[lat, lon, h] = ned2geodetic(north, east, down, lat0, lon0, h0, spheroid);

true_lla = [lat; lon; h];
end