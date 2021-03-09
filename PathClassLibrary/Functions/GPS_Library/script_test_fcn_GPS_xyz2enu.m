% script_test_fcn_GPS_xyz2enu.m
% Tests fcn_GPS_xyz2enu

% Revision history:
%   2021_01_14:
%       -- wrote the code

close all
clc
clear all %#ok<CLALL> % Clear any old variables

%% Test case 1: basic call
path_XYZ = [1016990.6723, -4727731.8284, 4145310.3562];
reference_LLA = [40.7934, -77.8600, 351.7392];
path_ENU = fcn_GPS_xyz2enu(path_XYZ, reference_LLA);
fprintf(1, 'INPUT:  Point in ECEF coordinate system is ( %.4f, %.4f, %.4f ) \n', path_XYZ)
fprintf(1, 'OUTPUT: Point in ENU coordinate system is ( %.4f, %.4f, %.4f ) \n', path_ENU)

%% Test case 2: NaN are not valid input [ERROR]
path_XYZ = [1016990.6723, -4727731.8284, NaN];
reference_LLA = [40.7934, -77.8600, 351.7392];
path_ENU = fcn_GPS_xyz2enu(path_XYZ, reference_LLA); %#ok<NASGU>

%% Test case 3: String/Char are not valid input [ERROR]
path_XYZ = [1016990.6723, -4727731.8284, 'a'];
reference_LLA = [40.7934, -77.8600, 351.7392];
path_ENU = fcn_GPS_xyz2enu(path_XYZ, reference_LLA); %#ok<NASGU>

%% Test case 4: Incorrect size of input in second dimension [ERROR]
path_XYZ = [1016990.6723, -4727731.8284];
reference_LLA = [40.7934, -77.8600, 351.7392];
path_ENU = fcn_GPS_xyz2enu(path_XYZ, reference_LLA); %#ok<NASGU>

%% Test case 5: basic call with just two points
path_XYZ = [1016990.6723, -4727731.8284, 4145310.3562;...
             977196.6307, -4756316.1339, 4122279.0630];
reference_LLA = [40.7934, -77.8600, 351.7392];
fig_num = 12345;
path_ENU = fcn_GPS_xyz2enu(path_XYZ, reference_LLA, fig_num); %#ok<NASGU>

%% Test case 6: basic call with multiple points
path_XYZ = readmatrix('sample_path_XYZ_data.csv');
reference_LLA = [40.7934, -77.8600, 351.7392];
fig_num = 12346;
path_ENU = fcn_GPS_xyz2enu(path_XYZ, reference_LLA, fig_num);
