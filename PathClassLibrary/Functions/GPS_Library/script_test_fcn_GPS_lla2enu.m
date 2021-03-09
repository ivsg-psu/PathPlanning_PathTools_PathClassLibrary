% script_test_fcn_GPS_lla2enu.m
% Tests fcn_GPS_lla2enu

% Revision history:
%   2021_01_14:
%       -- wrote the code

close all
clc
clear all %#ok<CLALL> % Clear any old variables

%% Test case 1: basic call with just two points
path_LLA = [40.7934, -77.8600, 351.7392;...
            40.52, -78.39, 355];
reference_LLA = [40.7934, -77.8600, 351.7392];
fig_num = 12345;
path_ENU = fcn_GPS_lla2enu(path_LLA, reference_LLA, fig_num); %#ok<NASGU>

%% Test case 2: basic call with just multiple points
path_LLA = readmatrix('sample_path_LLA_data.csv');
reference_LLA = [40.7934, -77.8600, 351.7392];
fig_num = 12346;
path_ENU = fcn_GPS_lla2enu(path_LLA, reference_LLA, fig_num);