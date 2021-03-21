% script_test_fcn_Path_addElevationToPath.m
% This is a script to exercise the function: fcn_Path_addElevationToPath.m
% This function was written on 2021_03_06 by Satya Prasad, szm888@psu.edu

% Revision history:
% 2021_03_20 - by S. Brennan
% - Added another example that is a bit more clearly "above" the query!
 
close all;

%% BASIC example 1
point = [0.5 0.2; 1.4 1.3];
reference_elevated_path = [0 0 0.1; 1 0 0.2; 2 0 0.3; 2 1 0.4];

fignum = 111;
elevated_path = fcn_Path_addElevationToPath(point, reference_elevated_path, fignum); %#ok<NASGU>

%% BASIC example 1.2 - works
point = [0.5 0.2; 1.4 1.3]; % Define the query point as an XY
reference_elevated_path = [0 0 0.1; 0.5 0.2 0.2; 0.9 0.9 0.3; 1.5 0.6 0.4; 3 0 0.5]; % Define an XYZ path
fignum = 112; % Define the figure number

% Snap the point onto the path
elevated_path = fcn_Path_addElevationToPath(point, reference_elevated_path, fignum);

%% BASIC example 1.3 - works
point = [0.5 0.2; 1.4 1.3]; % Define the query point as an XY
reference_elevated_path = [0 0 0.1; 0.25 0.2 0.2; 0.9 0.9 0.3; 1.1 1.1 0.4; 2.3 2.7 0.5]; % Define an XYZ path
fignum = 113; % Define the figure number

% Snap the point onto the path
elevated_path = fcn_Path_addElevationToPath(point, reference_elevated_path, fignum);
