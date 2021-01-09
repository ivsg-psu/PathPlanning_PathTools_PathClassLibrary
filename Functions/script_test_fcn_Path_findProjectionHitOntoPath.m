% script_test_fcn_Path_findProjectionHitOntoPath
% This is a script to exercise the function: fcn_Path_findProjectionHitOntoPath.m
% This function was written on 2020_11_10 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Modification history:
%      2020_12_31
%      -- updated for new argument list
%      2021_01_08
%      -- updated comments
%     2020_01_09
%     -- added more comments during clean-up

clc
close all

%% Simple test 1 - a simple intersection
fprintf(1,'Simple intersection result: \n');
path = [0 10; 10 10];
sensor_vector_start = [2 1]; 
sensor_vector_end   = [5 15];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 2 - no intersections
fprintf(1,'No intersection result: \n');
path = [-4 10; 2 10];
sensor_vector_start = [0 0]; 
sensor_vector_end   = [5 12];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 3 - multiple intersections
fprintf(1,'Multiple intersections result: \n');
path = [0 10; 10 10; 0 6; 10 6; 0 2];
sensor_vector_start = [0 0]; 
sensor_vector_end   = [5 12];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 4 - intersection through a vertex
fprintf(1,'Intersection through a vertex result: \n');
path = [0 5; 4 5; 8 2];
sensor_vector_start = [4 0]; 
sensor_vector_end   = [4 8];
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 5 - intersection at start of sensor
fprintf(1,'Intersection at start of sensor result: \n');
path = [0 5; 4 5; 8 2];
sensor_vector_start = [4 5]; 
sensor_vector_end   = [4 8];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 6 - intersection at end of sensor
fprintf(1,'Intersection at end of sensor result: \n');
path = [0 5; 4 5; 8 2];
sensor_vector_start = [4 0]; 
sensor_vector_end   = [4 5];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 7 - intersection beyond a sensor's range with flag
fprintf(1,'Intersection beyond sensor range result: \n');
path = [0 5; 4 5; 8 2];
sensor_vector_start = [4 0]; 
sensor_vector_end   = [4 2];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);


fig_debugging = 2344;
flag_search_type = 1;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

% Test the negative condition
sensor_vector_start = [4 6]; 
sensor_vector_end   = [4 8];
fig_debugging = 2345;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);


fig_debugging = 2346;
flag_search_type = 1;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%%
function print_results(distance,location)
fprintf(1,'Distance \t Location X \t Location Y \n');
fprintf(1,'%.3f \t\t %.3f \t\t\t %.3f\n',distance,location(:,1),location(:,2));
end
