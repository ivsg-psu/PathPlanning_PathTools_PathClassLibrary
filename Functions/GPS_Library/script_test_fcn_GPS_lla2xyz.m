% script_test_fcn_GPS_lla2xyz.m
% Tests fcn_GPS_lla2xyz

% Revision history:
%   2021_01_14:
%       -- wrote the code

close all
clc
clear all %#ok<CLALL> % Clear any old variables

%% Test case 1: basic call
point_LLA = [40.7934, -77.8600, 351.7392];
point_XYZ = fcn_GPS_lla2xyz(point_LLA);
fprintf(1, 'INPUT:  Point in Geodetic coordinate system is ( %.4f, %.4f, %.4f ) \n', point_LLA)
fprintf(1, 'OUTPUT: Point in ECEF coordinate system is ( %.4f, %.4f, %.4f ) \n', point_XYZ)

%% Test case 2: NaN are not valid input [ERROR]
point_LLA = [40.7934, -77.8600, NaN];
point_XYZ = fcn_GPS_lla2xyz(point_LLA); %#ok<NASGU>

%% Test case 3: String/Char are not valid input [ERROR]
point_LLA = [40.7934, -77.8600, 'a'];
point_XYZ = fcn_GPS_lla2xyz(point_LLA); %#ok<NASGU>

%% Test case 4: incorrect size of input in first dimension [ERROR]
point_LLA = [40.7934, -77.8600, 351.7392;...
             40.52, -78.39, 355];
point_XYZ = fcn_GPS_lla2xyz(point_LLA); %#ok<NASGU>

%% Test case 5: incorrect size of input in second dimension [ERROR]
point_LLA = [40.7934, -77.8600];
point_XYZ = fcn_GPS_lla2xyz(point_LLA); %#ok<NASGU>

%% Test case 6: latitude is not in the range [ERROR]
point_LLA = [100.7934, -77.8600, 351.7392];
point_XYZ = fcn_GPS_lla2xyz(point_LLA); %#ok<NASGU>

%% Test case 7: longitude is not in the range [ERROR]
point_LLA = [40.7934, -187.8600, 351.7392];
point_XYZ = fcn_GPS_lla2xyz(point_LLA);
