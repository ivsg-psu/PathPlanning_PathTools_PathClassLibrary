clear, clc
% Plot out the routes and determine the starting and ending index
load("IntersectionGPSData.mat");
start_idx = 1;
break_idx = 100000;
latitude = GPS.Latitude(start_idx:break_idx);
longitude = GPS.Longitude(start_idx:break_idx);

plot(longitude, latitude)
