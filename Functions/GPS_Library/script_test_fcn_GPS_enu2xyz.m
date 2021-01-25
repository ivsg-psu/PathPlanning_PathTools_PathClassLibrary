% script_test_fcn_GPS_enu2xyz.m
% Tests fcn_GPS_enu2xyz

% Revision history:
%   2021_01_14:
%       -- wrote the code

close all
clc
clear all %#ok<CLALL> % Clear any old variables

%% Test case 1: basic call
point_ENU = [-44915.4256, -30226.1806, -226.4526];
reference_LLA = [40.7934, -77.8600, 351.7392];
point_XYZ = fcn_GPS_enu2xyz(point_ENU, reference_LLA);
fprintf(1, 'INPUT:  Point in ENU coordinate system is ( %.4f, %.4f, %.4f ) \n', point_ENU)
fprintf(1, 'OUTPUT: Point in ECEF coordinate system is ( %.4f, %.4f, %.4f ) \n', point_XYZ)

%% Test case 2: NaN are not valid input [ERROR]
point_ENU = [-44915.4256, -30226.1806, NaN];
reference_LLA = [40.7934, -77.8600, 351.7392];
point_XYZ = fcn_GPS_enu2xyz(point_ENU, reference_LLA); %#ok<NASGU>

%% Test case 3: String/Char are not valid input [ERROR]
point_ENU = [-44915.4256, -30226.1806, 'a'];
reference_LLA = [40.7934, -77.8600, 351.7392];
point_XYZ = fcn_GPS_enu2xyz(point_ENU, reference_LLA); %#ok<NASGU>

%% Test case 4: incorrect size of input in first dimension [ERROR]
point_ENU = [-44915.4256, -30226.1806, -226.4526;...
             0, 0, 0];
reference_LLA = [40.7934, -77.8600, 351.7392];
point_XYZ = fcn_GPS_enu2xyz(point_ENU, reference_LLA); %#ok<NASGU>

%% Test case 5: incorrect size of input in second dimension [ERROR]
point_ENU = [-44915.4256, -30226.1806];
reference_LLA = [40.7934, -77.8600, 351.7392];
point_XYZ = fcn_GPS_enu2xyz(point_ENU, reference_LLA);
