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
%      2021_01_23
%      -- Added flag_search_type = 2 option, to allow multiple cross points
%      to be returned
%      -- Fixed bug with partially overlapping vectors not returning a
%      result
%      -- Added path segment output so that we can ID which segment was hit
%      2021_01_24
%      -- Fixed bug with overlapping colinear where two path segments
%      identified when there is only one

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

%% Simple test 2 - no intersections, returns NaN
% No intersection result: 
% Distance 	 Location X 	 Location Y 
% NaN 		 NaN 			 NaN

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

%% Simple test 3 - multiple intersections, returns only the first one
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


%% Simple test 7 - identically overlapping colinear 
fprintf(1,'identically overlapping colinear result: \n');
path = [0 10; 10 10];
sensor_vector_start = [0 10]; 
sensor_vector_end   = [10 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);
print_more_results(distance,location,path_segments);

%% Simple test 8 - partially overlapping colinear 1
fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10];
sensor_vector_start = [-2 10]; 
sensor_vector_end   = [10 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 9 - partially overlapping colinear 1
fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10];
sensor_vector_start = [-2 10]; 
sensor_vector_end   = [5 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 10 - partially overlapping colinear 1
fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10];
sensor_vector_start = [3 10]; 
sensor_vector_end   = [5 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 11 - partially overlapping colinear 1
fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10];
sensor_vector_start = [3 10]; 
sensor_vector_end   = [15 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 12 - super overlapping colinear 1
fprintf(1,'Super overlapping colinear result: \n');
path = [0 10; 10 10];
sensor_vector_start = [-3 10]; 
sensor_vector_end   = [15 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 13 - end overlapping colinear 1
fprintf(1,'End overlapping colinear result: \n');
path = [0 10; 10 10];
sensor_vector_start = [-3 10]; 
sensor_vector_end   = [0 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 14 - end overlapping colinear 2
fprintf(1,'End overlapping colinear result: \n');
path = [0 10; 10 10];
sensor_vector_start = [10 10]; 
sensor_vector_end   = [15 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);






%% Simple test 27 - identically overlapping colinear BACKWARDS
fprintf(1,'identically overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];
sensor_vector_end = [0 10]; 
sensor_vector_start   = [10 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 28 - partially overlapping colinear 1 BACKWARDS
fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];
sensor_vector_end = [-2 10]; 
sensor_vector_start   = [10 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 29 - partially overlapping colinear 1 BACKWARDS
fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];
sensor_vector_end = [-2 10]; 
sensor_vector_start   = [5 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 30 - partially overlapping colinear 1 BACKWARDS
fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];
sensor_vector_end = [3 10]; 
sensor_vector_start   = [5 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 31 - partially overlapping colinear 1 BACKWARDS
fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];
sensor_vector_end = [3 10]; 
sensor_vector_start   = [15 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 32 - super overlapping colinear 1 BACKWARDS
fprintf(1,'Super overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];
sensor_vector_end = [-3 10]; 
sensor_vector_start   = [15 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);


%% Simple test 33 - end overlapping colinear 1 BACKWARDS
fprintf(1,'End overlapping colinear result: \n');
path = [0 10; 10 10];
sensor_vector_end = [-3 10]; 
sensor_vector_start   = [0 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 34 - end overlapping colinear 2 BACKWARDS
fprintf(1,'End overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];
sensor_vector_end = [10 10]; 
sensor_vector_start   = [15 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);



%% Simple test 15 - non overlapping colinear 1
fprintf(1,'Non overlapping colinear result: \n');
path = [0 10; 10 10];
sensor_vector_start = [-3 10]; 
sensor_vector_end   = [-1 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Simple test 15 - non overlapping colinear 2
fprintf(1,'Non overlapping colinear result: \n');
path = [0 10; 10 10];
sensor_vector_start = [13 10]; 
sensor_vector_end   = [15 10];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);


%% Advanced test 1 - intersection beyond a sensor's range with flag
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

% Test showing that a sensor pointing away from a path "hits" the path with
% a negative distance
fig_debugging = 2346;
flag_search_type = 1;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%   __  __       _ _   _ _    _ _ _   
%  |  \/  |     | | | (_) |  | (_) |  
%  | \  / |_   _| | |_ _| |__| |_| |_ 
%  | |\/| | | | | | __| |  __  | | __|
%  | |  | | |_| | | |_| | |  | | | |_ 
%  |_|  |_|\__,_|_|\__|_|_|  |_|_|\__|
%                                     
%                                     


%% Advanced test 2 - multiple intersections
fprintf(1,'Single intersections reporting only first result: \n');
path = [0 10; 10 10; 0 6; 10 6; 0 2];
sensor_vector_start = [0 0]; 
sensor_vector_end   = [5 12];
fig_debugging = 23487;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

fprintf(1,'Multiple intersections reporting all results: \n');
path = [0 10; 10 10; 0 6; 10 6; 0 2];
sensor_vector_start = [0 0]; 
sensor_vector_end   = [5 12];
fig_debugging = 23488;
flag_search_type = 2;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced test 3 - multiple intersections possible, but no hits
fprintf(1,'Multiple intersections possible but no hits, reporting all results: \n');
path = [0 10; 10 10; 0 6; 10 6; 0 2];
sensor_vector_start = [0 0]; 
sensor_vector_end   = [0.5 1.2];
fig_debugging = 23499;
flag_search_type = 2;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced test 4 - multiple intersections possible, but few hits
fprintf(1,'Multiple intersections possible but few hits, reporting all results: \n');
path = [0 10; 10 10; 0 6; 10 6; 0 2];
sensor_vector_start = [0 0]; 
sensor_vector_end   = [2.5 6];
fig_debugging = 1010;
flag_search_type = 2;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);



%   __  __       _ _   _ _    _ _ _    ____                 _                   _             
%  |  \/  |     | | | (_) |  | (_) |  / __ \               | |                 (_)            
%  | \  / |_   _| | |_ _| |__| |_| |_| |  | |_   _____ _ __| | __ _ _ __  _ __  _ _ __   __ _ 
%  | |\/| | | | | | __| |  __  | | __| |  | \ \ / / _ \ '__| |/ _` | '_ \| '_ \| | '_ \ / _` |
%  | |  | | |_| | | |_| | |  | | | |_| |__| |\ V /  __/ |  | | (_| | |_) | |_) | | | | | (_| |
%  |_|  |_|\__,_|_|\__|_|_|  |_|_|\__|\____/  \_/ \___|_|  |_|\__,_| .__/| .__/|_|_| |_|\__, |
%                                                                  | |   | |             __/ |
%                                                                  |_|   |_|            |___/ 


%% Advanced Multihit Overlapping test - identically overlapping colinear 
fprintf(1,'identically overlapping colinear result: \n');
path = [0 10; 10 10];
sensor_vector_start = [0 10]; 
sensor_vector_end   = [10 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);
print_more_results(distance,location,path_segments);


%% Advanced Multihit Overlapping  test 8 - partially overlapping colinear 1
fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10];
sensor_vector_start = [-2 10]; 
sensor_vector_end   = [10 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 9 - partially overlapping colinear 1
fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10];
sensor_vector_start = [-2 10]; 
sensor_vector_end   = [5 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 10 - partially overlapping colinear 1
fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10];
sensor_vector_start = [3 10]; 
sensor_vector_end   = [5 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 11 - partially overlapping colinear 1
fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10];
sensor_vector_start = [3 10]; 
sensor_vector_end   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 12 - super overlapping colinear 1
fprintf(1,'Super overlapping colinear result: \n');
path = [0 10; 10 10];
sensor_vector_start = [-3 10]; 
sensor_vector_end   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 13 - end overlapping colinear 1
fprintf(1,'End overlapping colinear result: \n');
path = [0 10; 10 10];
sensor_vector_start = [-3 10]; 
sensor_vector_end   = [0 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 14 - end overlapping colinear 2
fprintf(1,'End overlapping colinear result: \n');
path = [0 10; 10 10];
sensor_vector_start = [10 10]; 
sensor_vector_end   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);






%% Advanced Multihit Overlapping  test 27 - identically overlapping colinear BACKWARDS
fprintf(1,'identically overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];
sensor_vector_end = [0 10]; 
sensor_vector_start   = [10 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 28 - partially overlapping colinear 1 BACKWARDS
fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];
sensor_vector_end = [-2 10]; 
sensor_vector_start   = [10 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 29 - partially overlapping colinear 1 BACKWARDS
fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];
sensor_vector_end = [-2 10]; 
sensor_vector_start   = [5 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 30 - partially overlapping colinear 1 BACKWARDS
fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];
sensor_vector_end = [3 10]; 
sensor_vector_start   = [5 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 31 - partially overlapping colinear 1 BACKWARDS
fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];
sensor_vector_end = [3 10]; 
sensor_vector_start   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 32 - super overlapping colinear 1 BACKWARDS
fprintf(1,'Super overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];
sensor_vector_end = [-3 10]; 
sensor_vector_start   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);


%% Advanced Multihit Overlapping  test 33 - end overlapping colinear 1 BACKWARDS
fprintf(1,'End overlapping colinear result: \n');
path = [0 10; 10 10];
sensor_vector_end = [-3 10]; 
sensor_vector_start   = [0 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 34 - end overlapping colinear 2 BACKWARDS
fprintf(1,'End overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10];
sensor_vector_end = [10 10]; 
sensor_vector_start   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);



%% Advanced Multihit Overlapping  test 15 - non overlapping colinear 1
fprintf(1,'Non overlapping colinear result: \n');
path = [0 10; 10 10];
sensor_vector_start = [-3 10]; 
sensor_vector_end   = [-1 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

%% Advanced Multihit Overlapping  test 15 - non overlapping colinear 2
fprintf(1,'Non overlapping colinear result: \n');
path = [0 10; 10 10];
sensor_vector_start = [13 10]; 
sensor_vector_end   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);



%   __  __       _ _   _ _    _ _ _   __  __       _ _ _   _____      _   _     
%  |  \/  |     | | | (_) |  | (_) | |  \/  |     | (_) | |  __ \    | | | |    
%  | \  / |_   _| | |_ _| |__| |_| |_| \  / |_   _| |_| |_| |__) |_ _| |_| |__  
%  | |\/| | | | | | __| |  __  | | __| |\/| | | | | | | __|  ___/ _` | __| '_ \ 
%  | |  | | |_| | | |_| | |  | | | |_| |  | | |_| | | | |_| |  | (_| | |_| | | |
%  |_|  |_|\__,_|_|\__|_|_|  |_|_|\__|_|  |_|\__,_|_|_|\__|_|   \__,_|\__|_| |_|
%                                                                               
%                                                                               


%% Advanced Multihit Overlapping test - identically overlapping colinear 
fprintf(1,'identically overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];
sensor_vector_start = [0 10]; 
sensor_vector_end   = [10 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 8 - partially overlapping colinear 1
fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10]; 
sensor_vector_start = [-2 10]; 
sensor_vector_end   = [10 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 9 - partially overlapping colinear 1
fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10]; 
sensor_vector_start = [-2 10]; 
sensor_vector_end   = [5 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 10 - partially overlapping colinear 1
fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10]; 
sensor_vector_start = [3 10]; 
sensor_vector_end   = [5 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 11 - partially overlapping colinear 1
fprintf(1,'Partially overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10]; 
sensor_vector_start = [3 10]; 
sensor_vector_end   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 12 - super overlapping colinear 1
fprintf(1,'Super overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];  
sensor_vector_start = [-3 10]; 
sensor_vector_end   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 13 - end overlapping colinear 1
fprintf(1,'End overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];  
sensor_vector_start = [-3 10]; 
sensor_vector_end   = [0 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 14 - end overlapping colinear 2
fprintf(1,'End overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];  
sensor_vector_start = [10 10]; 
sensor_vector_end   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);






%% Advanced Multihit Overlapping  test 27 - identically overlapping colinear BACKWARDS
fprintf(1,'identically overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];  
sensor_vector_end = [0 10]; 
sensor_vector_start   = [10 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 28 - partially overlapping colinear 1 BACKWARDS
fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];  
sensor_vector_end = [-2 10]; 
sensor_vector_start   = [10 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 29 - partially overlapping colinear 1 BACKWARDS
fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];  
sensor_vector_end = [-2 10]; 
sensor_vector_start   = [5 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 30 - partially overlapping colinear 1 BACKWARDS
fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];  
sensor_vector_end = [3 10]; 
sensor_vector_start   = [5 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 31 - partially overlapping colinear 1 BACKWARDS
fprintf(1,'Partially overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];  
sensor_vector_end = [3 10]; 
sensor_vector_start   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 32 - super overlapping colinear 1 BACKWARDS
fprintf(1,'Super overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];  
sensor_vector_end = [-3 10]; 
sensor_vector_start   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);


%% Advanced Multihit Overlapping  test 33 - end overlapping colinear 1 BACKWARDS
fprintf(1,'End overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];  
sensor_vector_end = [-3 10]; 
sensor_vector_start   = [0 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 34 - end overlapping colinear 2 BACKWARDS
fprintf(1,'End overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];  
sensor_vector_end = [10 10]; 
sensor_vector_start   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);



%% Advanced Multihit Overlapping  test 15 - non overlapping colinear 1
fprintf(1,'Non overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];  
sensor_vector_start = [-3 10]; 
sensor_vector_end   = [-1 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);

%% Advanced Multihit Overlapping  test 15 - non overlapping colinear 2
fprintf(1,'Non overlapping colinear result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];  
sensor_vector_start = [13 10]; 
sensor_vector_end   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);


%%
function print_results(distance,location)
fprintf(1,'Distance \t Location X \t Location Y \n');
if ~isempty(distance)
    for i_result = 1:length(distance(:,1))
        fprintf(1,'%.3f \t\t %.3f \t\t\t %.3f\n',distance(i_result),location(i_result,1),location(i_result,2));
    end
end
end

%%
function print_more_results(distance,location,path_segments)
fprintf(1,'Distance \t Location X \t Location Y \t PathSegment \n');
if ~isempty(distance)
    for i_result = 1:length(distance(:,1))
        fprintf(1,'%.3f \t\t %.3f \t\t\t %.3f \t\t %.0d\n',distance(i_result),location(i_result,1),location(i_result,2),path_segments(i_result));
    end
end
end
