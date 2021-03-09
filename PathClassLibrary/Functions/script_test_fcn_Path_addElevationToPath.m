% script_test_fcn_Path_addElevationToPath.m
% This is a script to exercise the function: fcn_Path_addElevationToPath.m
% This function was written on 2021_03_06 by Satya Prasad, szm888@psu.edu


close all;

%% BASIC example 1
point = [0.5 0.2; 1.4 1.3];
reference_elevated_path = [0 0 0.1; 1 0 0.2; 2 0 0.3; 2 1 0.4];

fignum = 111;
elevated_path = fcn_Path_addElevationToPath(point, reference_elevated_path, fignum); %#ok<NASGU>

%% BASIC example 1.2 - works
point = [0.5 0.2; 1.4 1.3]; % Define the query point
reference_elevated_path = [0 0 0.1; 0.5 0.2 0.2; 0.9 0.9 0.3; 1.5 0.6 0.4; 3 0 0.5]; % Define an XY path
fignum = 112; % Define the figure number

% Snap the point onto the path
elevated_path = fcn_Path_addElevationToPath(point, reference_elevated_path, fignum);
