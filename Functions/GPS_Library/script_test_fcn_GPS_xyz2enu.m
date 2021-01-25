% script_test_fcn_GPS_xyz2enu.m
% Tests fcn_GPS_xyz2enu

% Revision history:
%   2021_01_14:
%       -- wrote the code

close all
clc
clear all %#ok<CLALL> % Clear any old variables

%% Test case 1: basic call
point_XYZ = [1016990.6723, -4727731.8284, 4145310.3562];
reference_LLA = [40.7934, -77.8600, 351.7392];
point_ENU = fcn_GPS_xyz2enu(point_XYZ, reference_LLA);
fprintf(1, 'INPUT:  Point in ECEF coordinate system is ( %.4f, %.4f, %.4f ) \n', point_XYZ)
fprintf(1, 'OUTPUT: Point in ENU coordinate system is ( %.4f, %.4f, %.4f ) \n', point_ENU)

%% Test case 2: NaN are not valid input [ERROR]
point_XYZ = [1016990.6723, -4727731.8284, NaN];
reference_LLA = [40.7934, -77.8600, 351.7392];
point_ENU = fcn_GPS_xyz2enu(point_XYZ, reference_LLA); %#ok<NASGU>

%% Test case 3: String/Char are not valid input [ERROR]
point_XYZ = [1016990.6723, -4727731.8284, 'a'];
reference_LLA = [40.7934, -77.8600, 351.7392];
point_ENU = fcn_GPS_xyz2enu(point_XYZ, reference_LLA); %#ok<NASGU>

%% Test case 4: incorrect size of input in first dimension [ERROR]
point_XYZ = [1016990.6723, -4727731.8284, 4145310.3562;...
             977196.6307, -4756316.1339, 4122279.0630];
reference_LLA = [40.7934, -77.8600, 351.7392];
point_ENU = fcn_GPS_xyz2enu(point_XYZ, reference_LLA); %#ok<NASGU>

%% Test case 5: incorrect size of input in second dimension [ERROR]
point_XYZ = [1016990.6723, -4727731.8284];
reference_LLA = [40.7934, -77.8600, 351.7392];
point_ENU = fcn_GPS_xyz2enu(point_XYZ, reference_LLA);
