% script_test_fcn_GPS_enu2llaPath.m
% Tests fcn_GPS_enu2llaPath

% Revision history:
%   2021_01_14:
%       -- wrote the code

close all
clc
clear all %#ok<CLALL> % Clear any old variables

%% Test case 1: basic call with just two points
path_ENU = [-44915.4256, -30226.1806, -226.4526;...
             0, 0, 0];
reference_LLA = [140.7934, -77.8600, 351.7392];
fig_num = 12345;
path_LLA = fcn_GPS_enu2llaPath(path_ENU, reference_LLA, fig_num); %#ok<NASGU>

%% Test case 2: basic call with just multiple points
path_ENU = readmatrix('sample_path_ENU_data.csv');
reference_LLA = [40.7934, -77.8600, 351.7392];
fig_num = 12346;
path_LLA = fcn_GPS_enu2llaPath(path_ENU, reference_LLA, fig_num);