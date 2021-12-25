% script_test_fcn_GPS_lla2xyz.m
% Tests fcn_GPS_lla2xyz

% Revision history:
%   2021_01_14:
%       -- wrote the code

close all
clc
clear all %#ok<CLALL> % Clear any old variables

%% Test case 1: basic call
path_LLA = [40.7934, -77.8600, 351.7392];
path_XYZ = fcn_GPS_lla2xyz(path_LLA);
fprintf(1, 'INPUT:  Point in Geodetic coordinate system is ( %.4f, %.4f, %.4f ) \n', path_LLA)
fprintf(1, 'OUTPUT: Point in ECEF coordinate system is ( %.4f, %.4f, %.4f ) \n', path_XYZ)

%% Test case 2: NaN are not valid input [ERROR]
path_LLA = [40.7934, -77.8600, NaN];
path_XYZ = fcn_GPS_lla2xyz(path_LLA); %#ok<NASGU>

%% Test case 3: String/Char are not valid input [ERROR]
path_LLA = [40.7934, -77.8600, 'a'];
path_XYZ = fcn_GPS_lla2xyz(path_LLA); %#ok<NASGU>

%% Test case 4: Incorrect size of input in second dimension [ERROR]
path_LLA = [40.7934, -77.8600];
path_XYZ = fcn_GPS_lla2xyz(path_LLA); %#ok<NASGU>

%% Test case 5: latitude is not in the range [ERROR]
path_LLA = [100.7934, -77.8600, 351.7392];
path_XYZ = fcn_GPS_lla2xyz(path_LLA); %#ok<NASGU>

%% Test case 6: longitude is not in the range [ERROR]
path_LLA = [40.7934, -187.8600, 351.7392];
path_XYZ = fcn_GPS_lla2xyz(path_LLA); %#ok<NASGU>

%% Test case 7: basic call with just two points
path_LLA = [40.7934, -77.8600, 351.7392;...
            40.52, -78.39, 355];
fig_num = 12345;
path_XYZ = fcn_GPS_lla2xyz(path_LLA, fig_num); %#ok<NASGU>

%% Test case 8: basic call with multiple points
path_LLA = readmatrix('sample_path_LLA_data.csv');
fig_num = 12346;
path_XYZ = fcn_GPS_lla2xyz(path_LLA, fig_num);
