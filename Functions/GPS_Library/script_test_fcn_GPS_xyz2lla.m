% script_test_fcn_GPS_xyz2lla.m
% Tests fcn_GPS_xyz2lla

% Revision history:
%   2021_01_14:
%       -- wrote the code

close all
clc
clear all %#ok<CLALL> % Clear any old variables

%% Test case 1: basic call
point_XYZ = [1016990.6723, -4727731.8284, 4145310.3562];
point_LLA = fcn_GPS_xyz2lla(point_XYZ);
fprintf(1, 'INPUT:  Point in ECEF coordinate system is ( %.4f, %.4f, %.4f ) \n', point_XYZ)
fprintf(1, 'OUTPUT: Point in Geodetic coordinate system is ( %.4f, %.4f, %.4f ) \n', point_LLA)

%% Test case 2: NaN are not valid input [ERROR]
point_XYZ = [1016990.6723, -4727731.8284, NaN];
point_LLA = fcn_GPS_xyz2lla(point_XYZ); %#ok<NASGU>

%% Test case 3: String/Char are not valid input [ERROR]
point_XYZ = [1016990.6723, -4727731.8284, 'a'];
point_LLA = fcn_GPS_xyz2lla(point_XYZ); %#ok<NASGU>

%% Test case 4: incorrect size of input in first dimension [ERROR]
point_XYZ = [1016990.6723, -4727731.8284, 4145310.3562;...
             977196.6307, -4756316.1339, 4122279.0630];
point_LLA = fcn_GPS_xyz2lla(point_XYZ); %#ok<NASGU>

%% Test case 5: incorrect size of input in second dimension [ERROR]
point_XYZ = [1016990.6723, -4727731.8284];
point_LLA = fcn_GPS_xyz2lla(point_XYZ);
