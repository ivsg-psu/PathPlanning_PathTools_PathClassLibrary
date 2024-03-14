% script_test_fcn_Path_findProjectionHitOntoPath
% This is a script to exercise the function: fcn_Path_findProjectionHitOntoPath.m
% This function was written on 2020_11_10 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Modification history:
% 2020_12_31 - S. Brennan 
% -- updated for new argument list 
% 
% 2021_01_08 - S. Brennan
% -- updated comments
%
% 2020_01_09 - S. Brennan
% -- added more comments during clean-up
%
% 2021_01_23 - S. Brennan
% -- Added flag_search_type = 2 option, to allow multiple cross points to
% be returned
% -- Fixed bug with partially overlapping vectors not returning a result
% -- Added path segment output so that we can ID which segment was hit
%
% 2021_01_24 - S. Brennan
% -- Fixed bug with overlapping colinear where two path segments identified
% when there is only one
%
% 2021_01_24 - S. Brennan
% -- Added assertions to force checking
% -- Added test cases for flag=3 flag=4 options


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

assert(isequal(round(distance,4),9.2043));
assert(isequal(round(location,4),[3.9286,10.0000]));

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

assert(isnan(distance));
assert(all(isnan(location)));

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

assert(isequal(round(distance,4),2.6));
assert(isequal(round(location,4),[1.0000    2.4000]));


%% Simple test 4 - intersection through a vertex
fig_debugging = 345454;

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

assert(isequal(round(distance,4),5));
assert(isequal(round(location,4),[4 5]));


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

assert(isequal(round(distance,4),0));
assert(isequal(round(location,4),[4 5]));

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

assert(isequal(round(distance,4),5));
assert(isequal(round(location,4),[4 5]));


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

assert(isequal(round(distance,4),0));
assert(isequal(round(location,4),[0 10]));


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


assert(isequal(round(distance,4),2));
assert(isequal(round(location,4),[0 10]));

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

assert(isequal(round(distance,4),2));
assert(isequal(round(location,4),[0 10]));

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


assert(isequal(round(distance,4),0));
assert(isequal(round(location,4),[3 10]));

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

assert(isequal(round(distance,4),0));
assert(isequal(round(location,4),[3 10]));

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


assert(isequal(round(distance,4),3));
assert(isequal(round(location,4),[0 10]));

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

assert(isequal(round(distance,4),3));
assert(isequal(round(location,4),[0 10]));

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

assert(isequal(round(distance,4),0));
assert(isequal(round(location,4),[10 10]));


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

assert(isequal(round(distance,4),0));
assert(isequal(round(location,4),[10 10]));

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

assert(isequal(round(distance,4),0));
assert(isequal(round(location,4),[10 10]));

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

assert(isequal(round(distance,4),0));
assert(isequal(round(location,4),[5 10]));

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

assert(isequal(round(distance,4),0));
assert(isequal(round(location,4),[5 10]));

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

assert(isequal(round(distance,4),5));
assert(isequal(round(location,4),[10 10]));

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

assert(isequal(round(distance,4),5));
assert(isequal(round(location,4),[10 10]));


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

assert(isequal(round(distance,4),0));
assert(isequal(round(location,4),[0 10]));

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


assert(isequal(round(distance,4),5));
assert(isequal(round(location,4),[10 10]));


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

assert(isnan(distance));
assert(all(isnan(location)));

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

assert(isnan(distance));
assert(all(isnan(location)));

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

assert(isnan(distance));
assert(all(isnan(location)));


fig_debugging = 2344;
flag_search_type = 1;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);


assert(isequal(round(distance,4),5));
assert(isequal(round(location,4),[4 5]));

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


assert(isnan(distance));
assert(all(isnan(location)));

% Test showing that a sensor pointing away from a path "hits" the path with
% a negative distance
fig_debugging = 2346;
flag_search_type = 1;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

assert(isequal(round(distance,4),-1));
assert(isequal(round(location,4),[4 5]));

%% Multi-hit tests
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

assert(isequal(round(distance,4),2.6));
assert(isequal(round(location,4),[1 2.4]));

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

assert(isequal(round(distance',4),[10.8333    7.8000    6.5000    2.6000]));
assert(isequal(round(location,4),[    4.1667   10.0000
    3.0000    7.2000
    2.5000    6.0000
    1.0000    2.4000]));

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

assert(isempty(distance));
assert(isempty(location));

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

assert(isequal(round(distance,4),[6.5; 2.6]));
assert(isequal(round(location,4),[2.5000    6.0000
    1.0000    2.4000]));


%% Multi-hit overlapping tests
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

assert(isequal(round(distance,4),[0; 10]));
assert(isequal(round(location,4),[     0    10
    10    10]));


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

assert(isequal(round(distance,4),[2; 12]));
assert(isequal(round(location,4),[     0    10
    10    10]));

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

assert(isequal(round(distance,4),[2; 7]));
assert(isequal(round(location,4),[0   10.0000
    5.0000   10.0000]));

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

assert(isequal(round(distance,4),[0; 2]));
assert(isequal(round(location,4),[     3    10
     5    10]));

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

assert(isequal(round(distance,4),[0; 7]));
assert(isequal(round(location,4),[     3    10
    10    10]));

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

assert(isequal(round(distance,4),[3; 13]));
assert(isequal(round(location,4),[0    10
    10    10]));

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

assert(isequal(round(distance,4),3));
assert(isequal(round(location,4),[ 0    10]));

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

assert(isequal(round(distance,4),0));
assert(isequal(round(location,4),[10 10]));


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

assert(isequal(round(distance,4),[0; 10]));
assert(isequal(round(location,4),[    10    10
     0    10]));

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

assert(isequal(round(distance,4),[0; 10]));
assert(isequal(round(location,4),[    10    10
     0    10]));

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

assert(isequal(round(distance,4),[0; 5]));
assert(isequal(round(location,4),[    5    10
     0    10]));

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

assert(isequal(round(distance,4),[0; 2]));
assert(isequal(round(location,4),[    5    10
     3    10]));

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

assert(isequal(round(distance,4),[5; 12]));
assert(isequal(round(location,4),[    10    10
     3    10]));

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

assert(isequal(round(distance,4),[5; 15]));
assert(isequal(round(location,4),[    10    10
     0    10]));

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

assert(isequal(round(distance,4),0));
assert(isequal(round(location,4),[    0    10]));

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

assert(isequal(round(distance,4),5));
assert(isequal(round(location,4),[10    10]));


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

assert(isempty(distance));
assert(isempty(location));

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

assert(isempty(distance));
assert(isempty(location));


%% multiple hits with multiple paths
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

assert(isequal(round(distance,4),[     0
    10
    10]));
assert(isequal(round(location,4),[     0    10
    10    10
    10    10]));

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

assert(isequal(round(distance,4),[     2
    12
    12]));
assert(isequal(round(location,4),[     0    10
    10    10
    10    10]));

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

assert(isequal(round(distance,4),[    2.0000
    7.0000]));
assert(isequal(round(location,4),[         0   10.0000
    5.0000   10.0000]));

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

assert(isequal(round(distance,4),[0; 2]));
assert(isequal(round(location,4),[     3    10
     5    10]));

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

assert(isequal(round(distance,4),[     0
     7
    11
    11
     7
    12]));
assert(isequal(round(location,4),[     3    10
    10    10
    14    10
    14    10
    10    10
    15    10]));

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

assert(isequal(round(distance,4),[     3
    13
    17
    17
    13
    18]));
assert(isequal(round(location,4),[     0    10
    10    10
    14    10
    14    10
    10    10
    15    10]));

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

assert(isequal(round(distance,4),3));
assert(isequal(round(location,4),[0 10]));

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

assert(isequal(round(distance,4),[     0
     0
     4
     4
     5]));
assert(isequal(round(location,4),[    10    10
    10    10
    14    10
    14    10
    15    10]));

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

assert(isequal(round(distance,4),[          0
     0
    10]));
assert(isequal(round(location,4),[    10    10
    10    10
     0    10]));

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

assert(isequal(round(distance,4),[     0
     0
    10]));
assert(isequal(round(location,4),[    10    10
    10    10
     0    10]));

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

assert(isequal(round(distance,4),[     0
     5]));
assert(isequal(round(location,4),[
     5    10
     0    10]));

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

assert(isequal(round(distance,4),[0; 2]));
assert(isequal(round(location,4),[     5    10
     3    10]));

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

assert(isequal(round(distance,4),[     5
     5
     1
     0
    12
     1]));
assert(isequal(round(location,4),[   10.0000   10.0000
   10.0000   10.0000
   14.0000   10.0000
   15.0000   10.0000
    3.0000   10.0000
   14.0000   10.0000]));

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

assert(isequal(round(distance,4),[     5
     5
     1
     0
    15
     1]));
assert(isequal(round(location,4),[    10    10
    10    10
    14    10
    15    10
     0    10
    14    10]));


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

assert(isequal(round(distance,4),0));
assert(isequal(round(location,4),[0 10]));

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

assert(isequal(round(distance,4),[     5
     5
     1
     0
     1]));
assert(isequal(round(location,4),[    10    10
    10    10
    14    10
    15    10
    14    10]));

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

assert(isempty(distance));
assert(isempty(location));

%% Advanced Multihit Overlapping  test 15 - partially non-overlapping colinear 2
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

assert(isequal(round(distance,4),[     1
     1
     2]));
assert(isequal(round(location,4),[    14    10
    14    10
    15    10]));

%% flag equal to 3 and 4 tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ______ _                          ____
% |  ____| |               ______   |___ \
% | |__  | | __ _  __ _   |______|    __) |
% |  __| | |/ _` |/ _` |   ______    |__ <
% | |    | | (_| | (_| |  |______|   ___) |
% |_|    |_|\__,_|\__, |            |____/
%                  __/ |
%                 |___/
%                  _
%                 | |
%   __ _ _ __   __| |
%  / _` | '_ \ / _` |
% | (_| | | | | (_| |
%  \__,_|_| |_|\__,_|
%
%
%  ______ _                         _  _
% |  ____| |               ______  | || |
% | |__  | | __ _  __ _   |______| | || |_
% |  __| | |/ _` |/ _` |   ______  |__   _|
% | |    | | (_| | (_| |  |______|    | |
% |_|    |_|\__,_|\__, |              |_|
%                  __/ |
%                 |___/
%
%
%   ___ __ _ ___  ___  ___
%  / __/ _` / __|/ _ \/ __|
% | (_| (_| \__ \  __/\__ \
%  \___\__,_|___/\___||___/
%
% http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Flag%20%20%3D%20%203%0Aand%0AFlag%20%20%3D%204%0Acases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simple test - Extend wall vector (for flag_search_type = 3)
fig_debugging = 34341;
fprintf(1,'Simple test result, flag = 3, long sensor: \n');
path = [-4 10; 2 10];
sensor_vector_start = [0 0];
sensor_vector_end   = 2*[4 6];

flag_search_type =3;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

assert(isequal(round(distance,4),12.0185));
assert(isequal(round(location,4),[6.6667   10.0000]));
assert(isequal(path_segments,1));

%% Simple test - No intersection (for flag_search_type = 3)
fig_debugging = 34342;
fprintf(1,'Simple test result, flag = 3, short sensor: \n');
path = [-4 10; 2 10];
sensor_vector_start = [0 0];
sensor_vector_end   = [4 6];

flag_search_type =3;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

assert(isnan(distance));
assert(all(isnan(location)));
assert(isnan(path_segments));

%% Simple test - Extend both the vectors (for flag_search_type = 4)
fig_debugging = 34343;
fprintf(1,'Simple test result, flag = 4, long sensor: \n');
path = [-4 10; 2 10];
sensor_vector_start = [0 0];
sensor_vector_end   = 2*[4 6];

flag_search_type =4;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

assert(isequal(round(distance,4),12.0185));
assert(isequal(round(location,4),[6.6667   10.0000]));
assert(isequal(path_segments,1));

%% Simple test - Extend both the vectors (for flag_search_type = 4)
fig_debugging = 34344;
fprintf(1,'Simple test result, flag = 4, short sensor: \n');
path = [-4 10; 2 10];
sensor_vector_start = [0 0];
sensor_vector_end   = [4 6];

flag_search_type =4;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

assert(isequal(round(distance,4),12.0185));
assert(isequal(round(location,4),[6.6667   10.0000]));
assert(isequal(path_segments,1));

%% Parallel test for flag_search_type = 3
fig_debugging = 3434666;
fprintf(1,'Parallel test for flag=3: \n');
path = [0 0; 10 0];
sensor_vector_start = [0 1];
sensor_vector_end   = [10 1];

flag_search_type =3;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

assert(isnan(distance));
assert(all(isnan(location)));
assert(isnan(path_segments));

%% Parallel test for flag_search_type = 4
fig_debugging = 3434777;
fprintf(1,'Parallel test for flag=4: \n');
path = [0 0; 10 0];
sensor_vector_start = [0 1];
sensor_vector_end   = [10 1];

flag_search_type =4;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

assert(isnan(distance));
assert(all(isnan(location)));
assert(isnan(path_segments));


%% Multi-hit test - flag_search_type = 3
fig_debugging = 343433;
fprintf(1,'No intersection result: \n');
path = [-4 0; -2 0; 0 -2; 2 0; 4 0];
sensor_vector_start = [0 0];
sensor_vector_end   = [0 4];

flag_search_type = 3;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

assert(isequal(round(distance,4),0));
assert(isequal(round(location,4),[0   0]));
assert(isequal(path_segments,1));

%% Multi-hit test - flag_search_type = 4
fig_debugging = 343434;
fprintf(1,'No intersection result: \n');
path = [-4 0; -2 0; 0 -2; 2 0; 4 0];
sensor_vector_start = [0 0];
sensor_vector_end   = [0 4];

flag_search_type = 4;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

assert(isequal(round(distance,4),0));
assert(isequal(round(location,4),[0   0]));
assert(isequal(path_segments,1));

%% Multi-hit test - flag_search_type = 3
fig_debugging = 343437;
fprintf(1,'No intersection result: \n');
path = [-4 0; -2 0; 0 -2; 2 0; 4 0];
sensor_vector_start = [0 1];
sensor_vector_end   = [0 4];

flag_search_type = 3;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

assert(isnan(distance));
assert(all(isnan(location)));
assert(isnan(path_segments));


%% Multi-hit test - flag_search_type = 4
fig_debugging = 343438;
fprintf(1,'No intersection result: \n');
path = [-4 0; -2 0; 0 -2; 2 0; 4 0];
sensor_vector_start = [0 1];
sensor_vector_end   = [0 4];

flag_search_type = 4;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

assert(isequal(round(distance,4),-1));
assert(isequal(round(location,4),[0   0]));
assert(isequal(path_segments,1));

%% Advanced test 1 - intersection beyond a sensor's range with flag (for flag_search_type = 4)
fprintf(1,'Intersection beyond sensor range result: \n');
path = [0 5; 4 5; 8 2];

sensor_vector_start = [4 0];
sensor_vector_end   = [4 2];
fig_debugging = 34345;

flag_search_type =0;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
title('Intersection beyond sensor range result, Flag = 0')
print_results(distance,location);

assert(isnan(distance));
assert(all(isnan(location)));
assert(isnan(path_segments));

fig_debugging = 34346;

flag_search_type = 4;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
title('Intersection beyond sensor range result, Flag = 4')
print_results(distance,location);

assert(isequal(round(distance,4),5));
assert(isequal(round(location,4),[4 5]));
assert(isequal(path_segments,1));

% Test the negative condition
sensor_vector_start = [4 6];
sensor_vector_end   = [4 8];
fig_debugging = 34347;

flag_search_type = 0;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
title('Intersection beyond sensor range result, Flag = 0')
print_results(distance,location);

assert(isnan(distance));
assert(all(isnan(location)));
assert(isnan(path_segments));

fig_debugging = 34348;
flag_search_type = 4;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
title('Intersection beyond sensor range result, Flag = 4')
print_results(distance,location);

assert(isequal(round(distance,4),-1));
assert(isequal(round(location,4),[4 5]));
assert(isequal(path_segments,1));


%% Test of fast implementation mode 
% NOT IMPLEMENTED YET
% 
% % Perform the calculation in slow mode
% fig_num = [];
% REPS = 100; minTimeSlow = Inf; 
% tic;
% for i=1:REPS
%     tstart = tic;
%     [distance,location] = ...
%     fcn_geometry_findIntersectionOfSegments(...
%     wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
%     flag_search_type, (fig_num));
%     telapsed = toc(tstart);
%     minTimeSlow = min(telapsed,minTimeSlow);
% end
% averageTimeSlow = toc/REPS;
% 
% % Perform the operation in fast mode
% fig_num = -1;
% minTimeFast = Inf; nsum = 10;
% tic;
% for i=1:REPS
%     tstart = tic;
%     [distance,location] = ...
%     fcn_geometry_findIntersectionOfSegments(...
%     wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
%     flag_search_type, (fig_num));
%     telapsed = toc(tstart);
%     minTimeFast = min(telapsed,minTimeFast);
% end
% averageTimeFast = toc/REPS;
% 
% fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_fitVectorToNPoints:\n');
% fprintf(1,'N repetitions: %.0d\n',REPS);
% fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
% fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
% fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
% fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
% fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
% fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);


%%
function print_results(distance,location)
fprintf(1,'\nDistance \t Location X \t Location Y \n');
if ~isempty(distance)
    for i_result = 1:length(distance(:,1))
        fprintf(1,'%.3f \t\t %.3f \t\t\t %.3f\n',distance(i_result),location(i_result,1),location(i_result,2));
    end
end
end

%%
function print_more_results(distance,location,path_segments)
fprintf(1,'\nDistance \t Location X \t Location Y \t PathSegment \n');
if ~isempty(distance)
    for i_result = 1:length(distance(:,1))
        fprintf(1,'%.3f \t\t %.3f \t\t\t %.3f \t\t %.0d\n',distance(i_result),location(i_result,1),location(i_result,2),path_segments(i_result));
    end
end
end
