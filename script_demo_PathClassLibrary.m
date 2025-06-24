% This is a demonstration script to show the primary functionality of the
% path class library.

% Revisions:
% 2022_08_19 S. Brennan, sbrennan@psu.edu
% -- added clearing variable data on an example, as it gave errors
% -- fixed errors with the station-averaging example
% 2023_04_22 S. Brennan, sbrennan@psu.edu
% -- updated loading conditions
% -- improved comments and README.MD
% 2023_06_05 S. Brennan, sbrennan@psu.edu""
% -- cleaned up the workspace codes to use functions, work on MacOS
% 2023_08_25 to 2-23_09_06 by S. Brennan
% -- added examples of XY to St conversions, and vice versa
% -- added tools to calculate the centerline of paths
% 2024_03_14 by S. Brennan
% -- added the fix to the intersection calculation code to allow
% intersection options where the path segments are extended, not just the
% sensor segment
% -- deprecated fcn_Path_snapPointOntoNearestPath, replaced it with fcn_Path_snapPointToPathViaVectors
% -- fixed 2D / 3D bugs in addElevation - this is still not working
% correctly, but at least not throwing errors (it is close to correct, but
% need assertion and careful testing)
% -- fixed many of the scripts that were failing during the automated
% script testing.
% 2024_09_26 - Sean Brennan
% -- Updated function fcn_INTERNAL_clearUtilitiesFromPathAndFolders
% 2025_06_14 by S. Brennan
% -- fixed small bug in fcn_Path_findProjectionHitOntoPath
% -- updated DebugTools_v2024_12_18
% -- updated script_test_all_functions to latest version
% -- rewrote fcn_Path_findProjectionHitOntoPath to be fcn_Path_findSensorHitOnWall

% TO-DO:
% 2024_03_14 - S. Brennan 
% - need to add environmental variable methods for debugging.
% See for example the Geometry class. This will require that EVERY function
% will be fixed.
% - Need to use the Debug libraray, not Path library, to check inputs.
% - fix addElevation
% 2024_05_15 - Aneesh Batchu
% -- Found a bug in "fcn_Path_findProjectionHitOntoPath". A test case to
% demonstrate the BUG was written in "script_test_fcn_Path_findProjectionHitOntoPath"

%% Prep the workspace
close all
clc

%% Dependencies and Setup of the Code
% The code requires several other libraries to work, namely the following
% 
% * DebugTools - the repo can be found at: https://github.com/ivsg-psu/Errata_Tutorials_DebugTools
% * PathClassLibrary - the repo can be found at: https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary
% * GPS - this is the library that converts from ENU to/from LLA
% * GetUserInputs - this 
% List what libraries we need, and where to find the codes for each
clear library_name library_folders library_url

ith_library = 1;
library_name{ith_library}    = 'DebugTools_v2024_12_18';
library_folders{ith_library} = {'Functions','Data'};
library_url{ith_library}     = 'https://github.com/ivsg-psu/Errata_Tutorials_DebugTools/archive/refs/tags/DebugTools_v2024_12_18.zip';

ith_library = ith_library+1;
library_name{ith_library}    = 'GPSClass_v2023_06_29';
library_folders{ith_library} = {'Functions'};
library_url{ith_library}     = 'https://github.com/ivsg-psu/FieldDataCollection_GPSRelatedCodes_GPSClass/archive/refs/tags/GPSClass_v2023_06_29.zip';


% ith_library = ith_library+1;
% library_name{ith_library}    = 'GetUserInputPath_v2023_02_01';
% library_folders{ith_library} = {''};
% library_url{ith_library}     = 'https://github.com/ivsg-psu/PathPlanning_PathTools_GetUserInputPath/blob/main/Releases/GetUserInputPath_v2023_02_01.zip?raw=true';
% 
% ith_library = ith_library+1;
% library_name{ith_library}    = 'AlignCoordinates_2023_03_29';
% library_folders{ith_library} = {'Functions'};
% library_url{ith_library}     = 'https://github.com/ivsg-psu/PathPlanning_GeomTools_AlignCoordinates/blob/main/Releases/AlignCoordinates_2023_03_29.zip?raw=true';


%% Clear paths and folders, if needed
if 1==0
    clear flag_PathClass_Folders_Initialized;
    fcn_INTERNAL_clearUtilitiesFromPathAndFolders;
end

%% Set environment flags for input checking
% These are values to set if we want to check inputs or do debugging
% setenv('MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS','1');
% setenv('MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG','1');
setenv('MATLABFLAG_PATHCLASS_FLAG_CHECK_INPUTS','1');
setenv('MATLABFLAG_PATHCLASS_FLAG_DO_DEBUG','0');


%% Do we need to set up the work space?
if ~exist('flag_PathClass_Folders_Initialized','var')
    this_project_folders = {'Functions','Data'};
    fcn_INTERNAL_initializeUtilities(library_name,library_folders,library_url,this_project_folders);  
    flag_PathClass_Folders_Initialized = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See http://patorjk.com/software/taag/#p=display&f=Big&t=Getting%20Started
% 
%    _____      _   _   _                _____ _             _           _ 
%   / ____|    | | | | (_)              / ____| |           | |         | |
%  | |  __  ___| |_| |_ _ _ __   __ _  | (___ | |_ __ _ _ __| |_ ___  __| |
%  | | |_ |/ _ \ __| __| | '_ \ / _` |  \___ \| __/ _` | '__| __/ _ \/ _` |
%  | |__| |  __/ |_| |_| | | | | (_| |  ____) | || (_| | |  | ||  __/ (_| |
%   \_____|\___|\__|\__|_|_| |_|\__, | |_____/ \__\__,_|_|   \__\___|\__,_|
%                                __/ |                                     
%                               |___/                                      
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Welcome to the PathClass library!')

%% Show how to load sample paths, fcn_Path_fillSamplePaths  
% Call the function to fill in an array of "path" type
paths_array = fcn_Path_fillSamplePaths;

% We can even save one of these as a single "path"
single_path = paths_array{1};

%% Allow the user to self-select a path, fcn_Path_fillPathViaUserInputs
% fcn_Path_fillPathViaUserInputs
% A function for the user to click on the figure to generate XY path.
% Points are collected and plotted until the user double clicks. If the
% user right-clicks anywhere in the plot, the last point is deleted. Once
% the user double-clicks, the results are output from the function.
%
% FORMAT: 
%
%      pathXY = fcn_Path_fillPathViaUserInputs(fig_num)

%% Show how to calculate the relative angle changes between path segments, fcn_Path_calcDiffAnglesBetweenPathSegments
% Avoids the use of atan function since this breaks near -180 degree point

% Pick first path as reference_traversal structure
paths_array = fcn_Path_fillSamplePaths;
paths_to_check = paths_array{1};
fig_num = 11111; 
diff_angles = fcn_Path_calcDiffAnglesBetweenPathSegments(paths_to_check,fig_num);


%% Show how to calculate the yaw angles along a path, fcn_Path_calcYawFromPathSegments

% Basic call with one path
fig_num = 22222;
yaw_angles = fcn_Path_calcYawFromPathSegments(paths_to_check,fig_num); %#ok<*NASGU>


% Multiple paths
for i_path = 1:length(paths_array)
    % Pick first path as reference_traversal structure
    path_to_check = paths_array{i_path};
    
    % Pick first path as reference_traversal structure
    traversal_to_check = fcn_Path_convertPathToTraversalStructure(paths_array{i_path});
    all_traversals.traversal{i_path} = traversal_to_check;
    
    fig_num = 33333;
    yaw_angles = fcn_Path_calcYawFromPathSegments(path_to_check,fig_num);
end

% Plot the results? (Note: they are plotted below as well)
if 1==1
    fig_num = 333331;
    fcn_Path_plotTraversalsYaw(all_traversals,fig_num);
    fig_num = 333332;
    fcn_Path_plotTraversalsXY(all_traversals,fig_num);
end

%% Show how to convert paths to traversals, fcn_Path_convertPathToTraversalStructure
fig_num = 44444;
traversal = fcn_Path_convertPathToTraversalStructure(single_path,fig_num);
simple_example.traversal{1} = traversal;

fcn_Path_plotTraversalsXY(simple_example,fig_num);
xlabel('X [m]');
ylabel('Y [m]');

%% Show how to plot a traversal, fcn_Path_convertPathToTraversalStructure
fig_num = 55555;

fcn_Path_plotTraversalsXY(simple_example,fig_num);
xlabel('X [m]');
ylabel('Y [m]');

%% Show how to plot the yaw angle of traversal, fcn_Path_plotTraversalsYaw

% Fill in some dummy data
paths_array = fcn_Path_fillSamplePaths;
 
clear data

% Preallocate
data.traversal{length(paths_array)} = [];

% Convert paths into traversals
for i_traveral = 1:length(paths_array)
    traversal = fcn_Path_convertPathToTraversalStructure(paths_array{i_traveral});
    data.traversal{i_traveral} = traversal;
end

fig_num = 666661;
fcn_Path_plotTraversalsXY(data,fig_num);


% Next, specify the figure number to show that it will NOT auto-label the
% axes if figure is already given and it puts the plots into this figure.
fig_num = 666662;
fcn_Path_plotTraversalsYaw(data,fig_num);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See http://patorjk.com/software/taag/#p=display&f=Big&t=Basic%20Path%20Operations
%   ____            _        _____      _   _        ____                       _   _                 
%  |  _ \          (_)      |  __ \    | | | |      / __ \                     | | (_)                
%  | |_) | __ _ ___ _  ___  | |__) |_ _| |_| |__   | |  | |_ __   ___ _ __ __ _| |_ _  ___  _ __  ___ 
%  |  _ < / _` / __| |/ __| |  ___/ _` | __| '_ \  | |  | | '_ \ / _ \ '__/ _` | __| |/ _ \| '_ \/ __|
%  | |_) | (_| \__ \ | (__  | |  | (_| | |_| | | | | |__| | |_) |  __/ | | (_| | |_| | (_) | | | \__ \
%  |____/ \__,_|___/_|\___| |_|   \__,_|\__|_| |_|  \____/| .__/ \___|_|  \__,_|\__|_|\___/|_| |_|___/
%                                                         | |                                         
%                                                         |_|                                         
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Show how to find intersection of a path onto another path, fcn_Path_findProjectionHitOntoPath
% This is the core function of the entire library
% There are several flag options:
%
%            0: return distance and location of first intersection only if
%            the given sensor_vector overlaps the path (this is the
%            default). Returns NaNs if no intersections are found.
%
%            1: return distane and location of first intersection if any
%            projection of the sensor vector, in any direction, hits the
%            path (in other words, if there is any intersection). Note that
%            distance returned will be negative if the nearest intersection
%            is in the opposite direction of the given sensor vector.
%
%            2: returns distances and locations as M x 1 and M x 2 vectors
%            respectively, where the M rows represent all the detected
%            intersections.
%


% Simple test 3 - multiple intersections
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
% print_results(distance,location);

% Now we set the flag to 1, which searches in any direction
% Test showing that a sensor pointing away from a path "hits" the path with
% a negative distance
fig_debugging = 2346;
flag_search_type = 1;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
% print_results(distance,location);

% Now we set the flag to 2, which will give ALL the intersections
% GIVES:
% Multiple intersections reporting all results: 
% Distance 	 Location X 	 Location Y 
% 10.833 		 4.167 			 10.000
% 7.800 		 3.000 			 7.200
% 6.500 		 2.500 			 6.000
% 2.600 		 1.000 			 2.400

% Advanced test 2 - multiple intersections
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
% print_results(distance,location);

% If the sensor and path overlap, thus producing infininite intersections,
% only the first and last overlap points are given
% Advanced Multihit Overlapping  test 10 - partially overlapping colinear 1
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
% print_results(distance,location);
% % Gives:
% Partially overlapping colinear result: 
% Distance 	 Location X 	 Location Y 
% 0.000 		 3.000 			 10.000
% 2.000 		 5.000 			 10.000

% The algorithm is quite robust to overlapping conditions
% Advanced Multihit Overlapping  test 32 - super overlapping colinear 1 BACKWARDS
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
% print_more_results(distance,location,path_segments);
% GIVES:
% Super overlapping colinear BACKWARDS result: 
% Distance 	 Location X 	 Location Y 	 PathSegment 
% 5.000 		 10.000 			 10.000 		 1
% 5.000 		 10.000 			 10.000 		 2
% 1.000 		 14.000 			 10.000 		 3
% 0.000 		 15.000 			 10.000 		 4
% 15.000 		 0.000 		     	 10.000 		 1
% 1.000 		 14.000 			 10.000 		 4

%% Finding the intersection of traversals, fcn_Path_findIntersectionsBetweenTraversals
% fcn_Path_findIntersectionsBetweenTraversals
% Given two traversals, finds the intersection points between the traverals
% and returns the results as points, station coordinates in 1, and station
% coordinates in 2

% BASIC example 5 - two crossings of 1 onto 2
fig_num = 5;

% Create a dummy path and convert it to traversal_1
traversal_1_path = [0 -1; 2 2;  5 -3];  
traversal_1 = fcn_Path_convertPathToTraversalStructure(traversal_1_path);

% Create a dummy path and convert it to traversal_2
traversal_2_path = [0 0; 7 2];  
traversal_2 = fcn_Path_convertPathToTraversalStructure(traversal_2_path);


% Calculate the function output
[intersection_points,...
    s_coordinates_in_traversal_1,...
    s_coordinates_in_traversal_2] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    traversal_1,...
    traversal_2, ...
    fig_num);


% BASIC example 14 - Traversal 2 is on top of traversal 1 for two areas
fig_num = 14;

% Create a dummy path and convert it to traversal_1
traversal_1_path = [0 0; 8 2];  
traversal_1 = fcn_Path_convertPathToTraversalStructure(traversal_1_path);

% Create a dummy path and convert it to traversal_2
traversal_2_path = [1 1; 2 0.5; 4 1; 5 3; 6 1.5; 7 1.75; 9 1];  
traversal_2 = fcn_Path_convertPathToTraversalStructure(traversal_2_path);


% Calculate the function output
[intersection_points,...
    s_coordinates_in_traversal_1,...
    s_coordinates_in_traversal_2] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    traversal_1,...
    traversal_2, ...
    fig_num);


%% Show how to snap a point onto a path, fcn_Path_snapPointToPathViaVectors
% This can be tricky when the point is not on a perpendicular projection of
% any path segment. Examples are given in the documentation PPT.

% BASIC example 1.2 - works
point = [1.4 1.3]; % Define the query point
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 1.5 0.6; 3 0]; % Define an XY path
fignum = 112; % Define the figure number
 
% Snap the point onto the path
[closest_path_point,s_coordinate,...
    first_path_point_index,second_path_point_index,...
    percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(point, pathXY,fignum);
 
% Print results to the workspace
fprintf(1,'Figure: %d\n',fignum);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);


%% Show how to snap a point onto a traversal, fcn_Path_snapPointOntoNearestTraversal
point = [0.5 0.2];
pathXY = [0 0; 1 0; 2 0; 2 1];
traversal = fcn_Path_convertPathToTraversalStructure(pathXY);

fignum = 111;
[closest_path_point,s_coordinate,path_point_yaw,first_path_point_index,...
    second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestTraversal(point, traversal,fignum);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fignum, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See http://patorjk.com/software/taag/#p=display&f=Big&t=Basic%20Path%20Operations
% 
%   _______   _                     _                ____                       _   _                 
%  |__   __| (_)                   (_)              / __ \                     | | (_)                
%     | |_ __ _ _ __ ___  _ __ ___  _ _ __   __ _  | |  | |_ __   ___ _ __ __ _| |_ _  ___  _ __  ___ 
%     | | '__| | '_ ` _ \| '_ ` _ \| | '_ \ / _` | | |  | | '_ \ / _ \ '__/ _` | __| |/ _ \| '_ \/ __|
%     | | |  | | | | | | | | | | | | | | | | (_| | | |__| | |_) |  __/ | | (_| | |_| | (_) | | | \__ \
%     |_|_|  |_|_| |_| |_|_| |_| |_|_|_| |_|\__, |  \____/| .__/ \___|_|  \__,_|\__|_|\___/|_| |_|___/
%                                            __/ |        | |                                         
%                                           |___/         |_|                                         
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Show how to keep only part of a traversal, fcn_Path_findTraversalStationSegment
% This can be tricky when the station coordinates are outside of those
% within the traveral. See the test function for examples.

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
for i_Path = 1:1  % length(paths)
    traversal = fcn_Path_convertPathToTraversalStructure(...
        paths_array{i_Path});
    data.traversal{i_Path} = traversal;
end

% Plot the results?
if 1==1
    fig_num = 12;
    fcn_Path_plotTraversalsYaw(data,fig_num);

    fig_num = 13;
    fcn_Path_plotTraversalsXY(data,fig_num);
end

    
% BASIC example 1
s_coord_start = 10;
s_coord_end   = 100;
fignum = 111;
[traversal_segment1,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, fignum);
fprintf(1,...
    'Figure: %d, flag_outside_start is: %d, flag_outside_end is: %d \n',...
    fignum, flag_outside_start,flag_outside_end);
title('Normal query - Test case #1');


%% The following function can remove "pinch points", fcn_Path_removePinchPointInTraversal
% Given a traversal with a pinch point - an area where the traversal
% suddenly bends back on itself before continuing - this function removes
% the pinch point

% Pinch points are where a path crosses back onto itself, creating a loop.
% The minimum traversal along the path to go from start to end is to NOT go
% through the loop, but to jump at the "pinch point", thereby avoiding the
% loop or self-crossing area. The resulting path is one without "pinch
% points". This function is very useful to clean up weird polytopes that
% are self-crossing when traversing around an edge, or to remove weird jogs
% in data such as when projecting orthogonally.

% Advanced test case 3: show a pinch point in practice
fig_num = 1111;
figure(fig_num);
clf;
axis equal

% Clear any old variables
clear all_traversals

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});
all_traversals.traversal{1} = reference_traversal;


% Plot the results? (Note: they are plotted below as well)
if 1==0
    fig_num = 12;
    fcn_Path_plotTraversalsYaw(all_traversals,fig_num);
    fig_num = 13;
    fcn_Path_plotTraversalsXY(all_traversals,fig_num);
end

% Grab the "curve" of the path
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1}(13:20,:));
offsets = (0:1:10)'; 
flag_rounding_type = 1;
offset_traversals = fcn_Path_fillOffsetTraversalsAboutTraversal(reference_traversal, offsets,flag_rounding_type, fig_num);

% Fill in an array of "fixed" traversals
clear fixed_traversals
for ith_traversal = 1:length(offset_traversals.traversal)
    traversal_with_pinch_point = offset_traversals.traversal{ith_traversal};
    [traversal_no_pinch_point] = ...
        fcn_Path_removePinchPointInTraversal(...
        traversal_with_pinch_point);
    fixed_traversals.traversal{ith_traversal} = traversal_no_pinch_point; 
end

% Plot the results
fixed_fig_num = 2222;
figure(fixed_fig_num);
clf;
axis equal
hold on;
plot(reference_traversal.X,reference_traversal.Y,'b','Linewidth',3);
fcn_Path_plotTraversalsXY(fixed_traversals,fixed_fig_num)



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See http://patorjk.com/software/taag/#p=display&f=Big&t=Projection%20of%20Paths
% 
%   _____           _           _   _                      __   _____      _   _         
%  |  __ \         (_)         | | (_)                    / _| |  __ \    | | | |        
%  | |__) | __ ___  _  ___  ___| |_ _  ___  _ __     ___ | |_  | |__) |_ _| |_| |__  ___ 
%  |  ___/ '__/ _ \| |/ _ \/ __| __| |/ _ \| '_ \   / _ \|  _| |  ___/ _` | __| '_ \/ __|
%  | |   | | | (_) | |  __/ (__| |_| | (_) | | | | | (_) | |   | |  | (_| | |_| | | \__ \
%  |_|   |_|  \___/| |\___|\___|\__|_|\___/|_| |_|  \___/|_|   |_|   \__,_|\__|_| |_|___/
%                 _/ |                                                                   
%                |__/                                                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Project orthogonal vectors at user-defined stations, fcn_Path_findOrthogonalTraversalVectorsAtStations
% Given a central traversal and a set of stations along that traversal,
% finds the unit normal vector on the central traveral at each station
% point.

% Set up data
close all
central_path = [0 0; 1 1; 2 0];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
stations = [0; 1; 2^0.5-0.1; 2^0.5; 2^0.5+.1; 2; central_traversal.Station(end)];

% AVERAGING example 1 - default setting
flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

fig_num = 11;  % Define the figure

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalTraversalVectorsAtStations(...
    stations,central_traversal,flag_rounding_type,fig_num);
title('Vertex projection via prior segment (default, flag=1)');

% AVERAGING example 2 - use following segment
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use averagae projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

fig_num = 12;  % Define the figure

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalTraversalVectorsAtStations(...
    stations,central_traversal,flag_rounding_type,fig_num);
title('Vertex projection via following segment (flag=2)');

% AVERAGING example 3 - use average of both segments
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

fig_num = 13;  % Define the figure

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalTraversalVectorsAtStations(...
    stations,central_traversal,flag_rounding_type,fig_num);title('Vertex projection via averaging prior and following segment at vertex (flag=3)');

% AVERAGING example 4 - use average always
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

fig_num = 14;  % Define the figure

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalTraversalVectorsAtStations(...
    stations,central_traversal,flag_rounding_type,fig_num);title('Vertex projection via averaging everywhere (flag=4)');


%% Find how close one traversal is to another at specific station points, fcn_Path_findOrthogonalHitFromTraversalToTraversal
% fcn_Path_findOrthogonalHitFromTraversalToTraversals
% Given a central traversal and a set of stations along that traversal,
% finds the location on nearby traversals that are closest to the central
% traveral at each station point. Closest is defined via an orthogonal
% projection (or modifications of orthogonal projections) from the central
% traversal outward toward nearby traversals. Both positive and negative
% projections are included. Positive projections are those that, in the
% cross-product between the station direction and sensor
% projection, have a positive result. If a distance is in the positive
% direction, it is reported as positive. In the negative direction, it is
% reported as negative.

%% BASIC example 1 - parallel lines, query is in middle area
stations = 1; % Define the station

% Create a dummy central path and convert it to a traversal
central_path = [0 0; 2 0];  
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);

% Define a "nearby" path and convert it to a traversal
nearby_path = [0 4; 2 4];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);

% Set default values
flag_rounding_type = 3;
search_radius = 5;
fig_num = 1;

% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,...
    central_traversal,nearby_traversal,...
    flag_rounding_type,search_radius,fig_num);

% Make sure function worked
assert(isequal(round(closest_path_point,4),[1     4]));
assert(isequal(round(distances,4),[4])); %#ok<*NBRAK>


% BASIC example 1.5 - parallel lines, negative
stations = 1; % Define the station

% Create a dummy central path and convert it to a traversal
central_path = [0 0; 2 0];  
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);

% Define a "nearby" path and convert it to a traversal
nearby_path = [0 -4; 2 -4];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);

% Set default values
flag_rounding_type = 3;
search_radius = 5;
fig_num = 1;

% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,...
    central_traversal,nearby_traversal,...
    flag_rounding_type,search_radius,fig_num);

% Make sure function worked
assert(isequal(round(closest_path_point,4),[1     -4]));
assert(isequal(round(distances,4),[-4]));




% Example showing effect of flags
% Set up data
central_path = [0 0; 1 1; 2 0];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [-1 0.5; 0.5 2; 1.5 2; 3 0.5];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);
stations = [0; 1; 2^0.5-0.1; 2^0.5; 2^0.5+.1; 2; central_traversal.Station(end)];
search_radius = 1.5; % Distance to search for nearby segments

% AVERAGING example 1 - default setting
fig_num = 11;
flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num); %#ok<*ASGLU>
title('Vertex projection via prior segment (default, flag=1)');

% AVERAGING example 2 - use following segment
fig_num = 12;
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
title('Vertex projection via following segment (flag=2)');

% AVERAGING example 3 - use average of both segments
fig_num = 13;
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
title('Vertex projection via averaging prior and following segment at vertex (flag=3)');

% AVERAGING example 4 - use average always
fig_num = 14;
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
title('Vertex projection via averaging everywhere (flag=4)');

%% Find closest distance of many traversals to one central traversal, fcn_Path_findOrthoScatterFromTraversalToTraversals
% Given a central traversal and a set of stations along that traversal,
% finds the locations on the nearby traversals that are closest to the central
% traveral at each station point. Closest is defined via an orthogonal
% projection (or modifications of orthogonal projections) from the central
% traversal outward toward nearby traversals. Note that this function is
% essentially the multi-traversal version of:
% fcn_Path_findOrthogonalHitFromTraversalToTraversal

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
for i_Path = 1:length(paths_array)
    traversal = ...
        fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    all_traversals.traversal{i_Path} = traversal;
end

% Plot the results? (Note: they are plotted below as well)
if 1==0
    fig_num = 12;
    fcn_Path_plotTraversalsYaw(all_traversals,fig_num);
    fig_num = 13;
    fcn_Path_plotTraversalsXY(all_traversals,fig_num);
end


% Test case 1: basic call
reference_traversal = all_traversals.traversal{2};
reference_station_points = (0:10:reference_traversal.Station(end))';
flag_rounding_type = 3; % Use average of projections at end points
search_radius = 7;
fig_num = 1;

[closestXs, closestYs, closestDistances] = fcn_Path_findOrthoScatterFromTraversalToTraversals( reference_station_points, reference_traversal, all_traversals, flag_rounding_type,search_radius,fig_num); 
   
figure(11);
histogram([closestDistances(:,1);closestDistances(:,3)],30);
title('Histogram of all orthogonal distance projections');

%% Find offset traversals from a traversal, fcn_Path_fillOffsetTraversalsAboutTraversal
% fills in an array of traversals about a reference traversal at
% user-defined offset distances.
%
% Clear any old variables
clear all_traversals

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Test case 4: show how "pinching" can happen
fig_num = 4;
% Grab the "curve" of the path
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1}(13:20,:));
offsets = (-10:1:10)'; 
offset_traversals = fcn_Path_fillOffsetTraversalsAboutTraversal(reference_traversal, offsets,fig_num);
axis equal;

%% Convert one traversal's XY coordinates into SY coordinates using a reference traversal
% fcn_Path_convertTraversalXYtoSy
% Given a reference traversal and a set of stations along that traversal,
% finds the location on each nearby traversal that is closest to the central
% traveral at each station point. Closest is defined via an orthogonal
% projection (or modifications of orthogonal projections) from the central
% traversal outward toward nearby traversals. The results are then
% presented as an array of lateral offsets.
%
% NOTE: this is just a reformulation of the function:
% fcn_Path_findOrthogonalHitFromTraversalToTraversals
% The main difference is the final plotting result, but otherwise it uses
% the same functionality.

clear all_traversals
close all;

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
for i_Path = 1:length(paths_array)
    traversal = ...
        fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    all_traversals.traversal{i_Path} = traversal;
end

% Plot the results? (Note: they are plotted below as well)
if 1==0
    %     fig_num = 12;
    %     fcn_Path_plotTraversalsYaw(all_traversals,fig_num);
    fig_num = 13;
    fcn_Path_plotTraversalsXY(all_traversals,fig_num);
    title('Original paths in XY');
end

% Test case 1: basic call
reference_traversal = all_traversals.traversal{2};
reference_station_points = (0:20:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 40;
fig_num = 1;

[closestXs, closestYs, closestDistances] = ...
       fcn_Path_convertTraversalXYtoSy(...
       reference_station_points, reference_traversal, all_traversals,...
       flag_rounding_type,search_radius,fig_num); %#ok<*ASGLU>

subplot(1,2,2); axis normal;

%% Advanced testing example of fcn_Path_convertTraversalXYtoSy
% Set up data
% close all
reference_path = [0 0; 1 1; 2 0];
reference_traversal = fcn_Path_convertPathToTraversalStructure(reference_path);
% stations = [0; 1; 2^0.5-0.1; 2^0.5; 2^0.5+.1; 2; central_traversal.Station(end)];
reference_station_points_1st_half = linspace(0,sum(reference_path(2,:).^2,2).^0.5,11)'; 
reference_station_points_2nd_half = linspace(sum(reference_path(2,:).^2,2).^0.5,reference_traversal.Station(end),10)';
reference_station_points = [reference_station_points_1st_half(1:end-1,:);reference_station_points_2nd_half];

clear data
% Load a test path that is challenging for this reference path
test_path = fcn_Path_fillSamplePaths(4);
test_traversal.traversal{1} = fcn_Path_convertPathToTraversalStructure(test_path);

flag_rounding_type = 3;
search_radius = 40;
fig_num = 33;

[closestXs, closestYs, closestDistances] = ...
    fcn_Path_convertTraversalXYtoSy(...
    reference_station_points, reference_traversal, test_traversal,...
    flag_rounding_type,search_radius,fig_num);

% Repeat with different projection type
flag_rounding_type = 4;
search_radius = 40;
fig_num = 44;

[closestXs, closestYs, closestDistances] = ...
    fcn_Path_convertTraversalXYtoSy(...
    reference_station_points, reference_traversal, test_traversal,...
    flag_rounding_type,search_radius,fig_num);

   
%% Find random traversals (with varying smoothness) about traversal:
% fcn_Path_fillRandomTraversalsAboutTraversal
% fills in random traversals about a reference traversal. Points are
% generated via orthogonal projection using random normal distribution with
% either a default variance or optional user-defined variance. The station
% points can also be user-specified as randomly distributed uniformly, or
% default to the station points in the reference_traversal if no optional
% inputs are given. The first and last stations are forced to be the same
% as the reference_traversal to prevent the route from randomly becoming
% shorter with repeated calls to this function.

% Clear any old variables
clear all_traversals

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});
all_traversals.traversal{1} = reference_traversal;


% Plot the results? (Note: they are plotted below as well)
if 1==0
    fig_num = 12;
    fcn_Path_plotTraversalsYaw(all_traversals,fig_num);
    fig_num = 13;
    fcn_Path_plotTraversalsXY(all_traversals,fig_num);
end

% Test case 1: basic call for one trajectory
random_traversals = fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal);

figure(1212);
plot(random_traversals.traversal{1}.X,random_traversals.traversal{1}.Y,'r.-','Linewidth',3);


% Test case 4: show effects of spatial smoothness with many trajectories
%      random_traversals = ...
%      fcn_Path_fillRandomTraversalsAboutTraversal(...
%            reference_traversal,...
%            (std_deviation),...
%            (num_trajectories),...
%            (num_points),...
%            (flag_generate_random_stations),...
%            (spatial_smoothness),...
%            (fig_num));
emtpy_value = [];
flag_generate_random_stations = 0;
num_trajectories = 5;

fig_num = 41;
spatial_smoothness = 2;  % Units are meters
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title('Spatial smoothness: 2 meters (generates warning)');

fig_num = 42;
spatial_smoothness = 5;  % Units are meters
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('Spatial smoothness: %.0d meters',spatial_smoothness));


fig_num = 43;
spatial_smoothness = 10;  % Units are meters
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('Spatial smoothness: %.0d meters',spatial_smoothness));

fig_num = 44;
spatial_smoothness = 15;  % Units are meters
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('Spatial smoothness: %.0d meters',spatial_smoothness));

fig_num = 45;
spatial_smoothness = 25;  % Units are meters
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('Spatial smoothness: %.0d meters',spatial_smoothness));


% Test case 5: show effects of standard deviation
%      random_traversals = ...
%      fcn_Path_fillRandomTraversalsAboutTraversal(...
%            reference_traversal,...
%            (std_deviation),...
%            (num_trajectories),...
%            (num_points),...
%            (flag_generate_random_stations),...
%            (spatial_smoothness),...
%            (fig_num));
emtpy_value = [];
flag_generate_random_stations = 0;
num_trajectories = 5;
spatial_smoothness = 7;  % Units are meters

fig_num = 51;
std_deviation = 1;  % Units are meters
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    std_deviation,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('Standard deviation: %.0d meters',std_deviation));

fig_num = 52;
std_deviation = 2;  % Units are meters
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    std_deviation,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('Standard deviation: %.0d meters',std_deviation));

fig_num = 53;
std_deviation = 5;  % Units are meters
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    std_deviation,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('Standard deviation: %.0d meters',std_deviation));


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See http://patorjk.com/software/taag/#p=display&f=Big&t=Statistical%20Analysis%20of%20Paths
%    _____ _        _   _     _   _           _                        _           _              __   _____      _   _         
%   / ____| |      | | (_)   | | (_)         | |     /\               | |         (_)            / _| |  __ \    | | | |        
%  | (___ | |_ __ _| |_ _ ___| |_ _  ___ __ _| |    /  \   _ __   __ _| |_   _ ___ _ ___    ___ | |_  | |__) |_ _| |_| |__  ___ 
%   \___ \| __/ _` | __| / __| __| |/ __/ _` | |   / /\ \ | '_ \ / _` | | | | / __| / __|  / _ \|  _| |  ___/ _` | __| '_ \/ __|
%   ____) | || (_| | |_| \__ \ |_| | (_| (_| | |  / ____ \| | | | (_| | | |_| \__ \ \__ \ | (_) | |   | |  | (_| | |_| | | \__ \
%  |_____/ \__\__,_|\__|_|___/\__|_|\___\__,_|_| /_/    \_\_| |_|\__,_|_|\__, |___/_|___/  \___/|_|   |_|   \__,_|\__|_| |_|___/
%                                                                         __/ |                                                 
%                                                                        |___/                                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plotting a band around a traversal, fcn_Path_plotTraversalXYWithUpperLowerBands
% Plots a traversal with a band defined by an upper and lower traversal.
% All traversals must have the same data length. 

% Test case 1: basic call for one trajectory
fig_num = 23333;
middle_path = [0 0; 1 1; 4 0; 5 0.5];
upper_path = [0 2; 1 2; 3.9 3; 4.9 1.5];
lower_path = [0 -1; 0.9 0; 4.1 -3; 5.1 0];
middle_traversal = fcn_Path_convertPathToTraversalStructure(middle_path);
upper_traversal = fcn_Path_convertPathToTraversalStructure(upper_path);
lower_traversal = fcn_Path_convertPathToTraversalStructure(lower_path);

fcn_Path_plotTraversalXYWithUpperLowerBands( middle_traversal, upper_traversal, lower_traversal, fig_num);

%% Approximating the lateral variance of a traversal, fcn_Path_calcSingleTraversalStandardDeviation
% fcn_Path_calcSingleTraversalStandardDeviation
% calculates the standard deviation in the offsets of a single traversal by
% analyzing the variance in angles along a reference_traversal, then
% multiplying these by the average segment length in the reference
% traversal. The resulting standard deviation approximates the
% variance in lateral offset that occurs at the end of each segment, versus
% a line projected from the previous segment.
% Clear any old variables
clear all_traversals

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});


% Test case 1: basic call for one trajectory
fig_num = 252525;
std_deviation = fcn_Path_calcSingleTraversalStandardDeviation(reference_traversal,fig_num);

%% Plotting variance bands about a traversal, fcn_Path_plotTraversalXYWithVarianceBands
% fcn_Path_plotTraversalXYWithVarianceBands
% Plots a traversal with a variance band around the path

% Clear any old variables
clear all_traversals

% Fill in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;

% Test case 4: advanced call for multiple trajectories
fig_num = 444444;
std_deviation = 2;
for i_Path = 1:length(paths)
    reference_traversal = fcn_Path_convertPathToTraversalStructure(paths{i_Path});
    fcn_Path_plotTraversalXYWithVarianceBands(reference_traversal,...
        std_deviation,fig_num);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See http://patorjk.com/software/taag/#p=display&f=Big&t=Path%20Averaging%20Methods
%   _____      _   _                                        _               __  __      _   _               _     
%  |  __ \    | | | |         /\                           (_)             |  \/  |    | | | |             | |    
%  | |__) |_ _| |_| |__      /  \__   _____ _ __ __ _  __ _ _ _ __   __ _  | \  / | ___| |_| |__   ___   __| |___ 
%  |  ___/ _` | __| '_ \    / /\ \ \ / / _ \ '__/ _` |/ _` | | '_ \ / _` | | |\/| |/ _ \ __| '_ \ / _ \ / _` / __|
%  | |  | (_| | |_| | | |  / ____ \ V /  __/ | | (_| | (_| | | | | | (_| | | |  | |  __/ |_| | | | (_) | (_| \__ \
%  |_|   \__,_|\__|_| |_| /_/    \_\_/ \___|_|  \__,_|\__, |_|_| |_|\__, | |_|  |_|\___|\__|_| |_|\___/ \__,_|___/
%                                                      __/ |         __/ |                                        
%                                                     |___/         |___/                                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;

% Fill in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
for i_Path = 1:length(paths)
    traversal = fcn_Path_convertPathToTraversalStructure(paths{i_Path});
    data.traversal{i_Path} = traversal;
end

% Plot the results?
if 1==0
    fig_num = 12;
    fcn_Path_plotTraversalsYaw(data,fig_num);
    fig_num = 13;
    fcn_Path_plotTraversalsXY(data,fig_num);
end

% %% Averaging by station, fcn_Path_findAverageTraversalViaStationAlignment
% % fcn_Path_findAverageTraversalViaStationAlignment
% % finds the average traversal of several traversals by averaging the
% % lateral offsets at the same stations.
% 
% [aligned_Data_ByStation,mean_Data] = ...
%     fcn_Path_findAverageTraversalViaStationAlignment(data);
% 
% % Plot the final XY result of mean station
% path_points_fig = 11111;
% fcn_Path_plotTraversalsXY(data,path_points_fig);
% hold on
% plot(mean_Data.mean_xEast,mean_Data.mean_yNorth,'Linewidth',4);
% title('Original paths and final average path via station averaging')
% xlabel('X [m]')
% ylabel('Y [m]')


%% Averaging by closest point, fcn_Path_findAverageTraversalViaClosestPoint
% fcn_Path_findAverageTraversalViaClosestPoint
% finds the average of several traversals by taking a reference traversal
% (or, if one is not given, using the traversal with longest number of
% points) and for each point in the traversal finding the nearest point in
% other traversals. The nearest point is determined by the "snap" of the
% reference traversal vertices to the closest location of each path. Thus,
% the resulting projection is orthogonal to each individual nearby path
% (and not usually orthogonal to the reference path).

path_average_final2 = fcn_Path_findAverageTraversalViaClosestPoint(data);

% Plot the final XY result of closest point
path_points_fig = 22222;
fcn_Path_plotTraversalsXY(data,path_points_fig);
hold on
plot(path_average_final2.X,path_average_final2.Y,'Linewidth',4);
title('Original paths and final average path via closest point averaging')
xlabel('X [m]')
ylabel('Y [m]')

%% Averaging by orthogonal projection, fcn_Path_findAverageTraversalViaOrthoProjection
% fcn_Path_findAverageTraversalViaOrthoProjection
% finds the average traversal of several traversals by taking a reference
% traversal or, if of a referemce traversal is not given, it uses as a
% reference the traversal with longest number of points. 
%
% As additional outputs, for each point in the traversal, this function
% also finds the intersection point in other traversals via orthogonal
% projection, and saves these points into arrays to denote the closest X
% and Y coordinates, and the distances. These X, Y, and distance arrays
% have M columns, one for each traversal, and N rows, one for each station
% in the reference traversal.

path_average_final3 = fcn_Path_findAverageTraversalViaOrthoProjection(data);

% Plot the final XY result of orthogonal
path_points_fig = 33333;
fcn_Path_plotTraversalsXY(data,path_points_fig);
hold on
plot(path_average_final3.X,path_average_final3.Y,'Linewidth',4);
title('Original paths and final average path via orthogonal projections')
xlabel('X [m]')
ylabel('Y [m]')


%% Plot the final XY results of all three
path_points_fig = 123;
figure(path_points_fig);
clf;
hold on
% plot(mean_Data.mean_xEast,mean_Data.mean_yNorth,'Linewidth',4);
plot(path_average_final2.X,path_average_final2.Y,'Linewidth',4);
plot(path_average_final3.X,path_average_final3.Y,'Linewidth',4);
fcn_Path_plotTraversalsXY(data,path_points_fig);
title('Original paths and final average paths');
legend('Average Station','Closest point','Orthogonal projection','Paths')
xlabel('X [m]')
ylabel('Y [m]')

%% The following are sub-functions that are used in the averaging methods


%% Find closest points on other traversals to a traversal, fcn_Path_findClosestPointsToTraversal
% [closestXs,closestYs,closestZs,closestYaws] = fcn_Path_findClosestPointsToTraversal(...
%    reference_traversal,data,varargin)
% This function finds the projection point on each traversal from a
% reference path by snapping each vertex of the reference traversal onto
% the nearby traversals. Each "snap" calculates the nearest point within
% each other traversal by measuring the orthogonal distance from the nearby
% traversal, to the point on the reference_traversal.

% Setup
clear data
central_path = [-2 1; 1 4; 3 2; 5 2; 6 3; 7 2];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [-1 0.5; 0.5 2; 1.5 2; 3 0.5; 4 3; 5 4; 6 3; 7 1];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);
data.traversal{1} = nearby_traversal;

% MULTICROSS example 1 - default setting
fig_num = 13425;
flag_yaw = 1;
flag_3D = 0;
[closestXs,closestYs,closestZs] = ...
    fcn_Path_findClosestPointsToTraversal(central_traversal,data,flag_yaw,flag_3D,fig_num);


%% Finding the traversal with the most data, fcn_Path_findTraversalWithMostData
% finds the traversal index with the most amount of data (determined as the
% most elements in the X array)

% Fill in some dummy data
paths_array = fcn_Path_fillSamplePaths;

% Convert paths into traversals
for i_traveral = 1:length(paths_array)
    traversal = fcn_Path_convertPathToTraversalStructure(paths_array{i_traveral});
    data.traversal{i_traveral} = traversal;
end

index_of_longest = fcn_Path_findTraversalWithMostData(data);
fprintf(1,'The longest path of the %.0d paths was path %.0d with %.0d elements\n',...
    length(data.traversal),...
    index_of_longest,...
    length(data.traversal{index_of_longest}.X));

%% Generating a new traversal via station resampling, fcn_Path_newTraversalByStationResampling
% fcn_Path_newTraversalByStationResampling
% creates a new traversal by resampling a given traversal at given station
% points. 
%
% Note: if the stations are intended to align in space between the
% input_traversal and new_traversal traversals, then the first station
% point must be zero.
%
% If the stations are outside the station range of the input traversal,
% then extraploation is used to extend the input_traversal linearly
% outward. This can result in bad data if the path is not approximately
% linear at the endpoints.

% Basic example 1 - start at zero
% Fill in sample paths (as a starter)
basic_path = [0 0; 10 0; 20 0];
input_traversal = fcn_Path_convertPathToTraversalStructure(basic_path);

fig_num = 23444;

% Redecimate the traversal at 1-meter increments
interval = 1;
new_stations    = (0:interval:5)';
new_traversal = fcn_Path_newTraversalByStationResampling(input_traversal, new_stations, fig_num);

% ADVANCED EXAMPLE
% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;
input_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

fig_num = 23445;

% Redecimate the traversal at 1-meter increments
interval = 10;
new_stations    = (0:interval:input_traversal.Station(end))';
new_traversal = fcn_Path_newTraversalByStationResampling(input_traversal, new_stations, fig_num);


%% Remove forward/backward jogs from paths, fcn_Path_cleanPathFromForwardBackwardJogs
% fcn_Path_cleanPathFromForwardBackwardJogs
% Finds and removes situations where the path is jumping forward and
% backward. This is detected by finding situations where the angle between
% segments is more than a threshold (currently pi/4), and then taking these
% segments, and the one before and after, and removing them. It then
% re-scans the path (up to 3 times) to again check for these situations.

% Example 1: Basic call 
fig_num = 1;
path_with_jogs = [0 0; 1 1; 2 2.2; 3.3 3; 2.5 2.7; 3.5 3.6; 5 5];
clean_path = fcn_Path_cleanPathFromForwardBackwardJogs...
    (path_with_jogs,fig_num);

plot(path_with_jogs(:,1), path_with_jogs(:,2),'.-','Linewidth',2,'Markersize',25);
plot(clean_path(:,1), clean_path(:,2),'.-','Linewidth',2,'Markersize',25);
title('Original path with jogs and cleaned path')
xlabel('X [m]')
ylabel('Y [m]')

%% Demonstration of fcn_Path_findCenterlineVoteFromTraversalToTraversal
% This function finds the center projected from one traversal toward
% another
from_path = [0 0; 1 1; 2 1; 3 4];
to_path   = from_path + ones(length(from_path(:,1)),1)*[0 1];
from_traversal =  fcn_Path_convertPathToTraversalStructure(from_path);
to_traversal =  fcn_Path_convertPathToTraversalStructure(to_path);
flag_rounding_type = 1;
search_radius = 10;
flag_project_full_distance = 0; % Set to 1 to do a full distance projection, not halfway
fig_num = fig_num+1;

[centerline_points_projected,unit_vectors_orthogonal] = ...
    fcn_Path_findCenterlineVoteFromTraversalToTraversal(...
    from_traversal,to_traversal,(flag_rounding_type),(search_radius),(flag_project_full_distance), (fig_num));


%% Demonstration of fcn_Path_findCenterPathBetweenTwoPaths
% This function finds the center projected from one traversal toward
% another
first_path = [0 0; 1 1; 2 1; 3 2];
second_path   = [0.5 1.5; 1.5 2.1; 4 6];

flag_rounding_type = 1;
search_radius = 10;
fig_num = fig_num+1;

figure(fig_num);
clf;
hold on;
grid on;
axis equal;

center_path = ...
    fcn_Path_findCenterPathBetweenTwoPaths(...
    first_path,second_path,(flag_rounding_type),(search_radius),(fig_num));


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See http://patorjk.com/software/taag/#p=display&f=Big&t=Elevated%20Paths
%   ______ _                 _           _   _____      _   _         
%  |  ____| |               | |         | | |  __ \    | | | |        
%  | |__  | | _____   ____ _| |_ ___  __| | | |__) |_ _| |_| |__  ___ 
%  |  __| | |/ _ \ \ / / _` | __/ _ \/ _` | |  ___/ _` | __| '_ \/ __|
%  | |____| |  __/\ V / (_| | ||  __/ (_| | | |  | (_| | |_| | | \__ \
%  |______|_|\___| \_/ \__,_|\__\___|\__,_| |_|   \__,_|\__|_| |_|___/
%                                                                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add elevation to a path by finding 2 nearest neighbors, fcn_Path_addElevationToPath 
% fcn_Path_addElevationToPath
% Adds elevation to the 2-dimensional 'path' based on the nearest neighbors
% in a 3-dimesional 'reference_elevated_path'
% 
% FORMAT:
%
%      elevated_path = fcn_Path_addElevationToPath(path, ...
%                      reference_elevated_path, (fig_num))

% BASIC example 1.3
point = [0.5 0.2; 1.4 1.3]; % Define the query point as an XY
reference_elevated_path = [0 0 0.1; 0.25 0.2 0.2; 0.9 0.9 0.3; 1.1 1.1 0.4; 2.3 2.7 0.5]; % Define an XYZ path
fignum = 113; % Define the figure number

% Snap the point onto the path
% elevated_path = fcn_Path_addElevationToPath(point, reference_elevated_path, fignum);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=XYZ%20%20to%20STH%20%0Aand%20%0ASTH%20%20to%20XYZ%20%20%20%20%0AConversions%0A
%  __   ____     ________   _           _____ _______ _    _     
%  \ \ / /\ \   / /___  /  | |         / ____|__   __| |  | |    
%   \ V /  \ \_/ /   / /   | |_ ___   | (___    | |  | |__| |    
%    > <    \   /   / /    | __/ _ \   \___ \   | |  |  __  |    
%   / . \    | |   / /__   | || (_) |  ____) |  | |  | |  | |    
%  /_/ \_\   |_|  /_____|   \__\___/  |_____/   |_|  |_|  |_|    
%                  | |                                           
%    __ _ _ __   __| |                                           
%   / _` | '_ \ / _` |                                           
%  | (_| | | | | (_| |                                           
%   \__,_|_| |_|\__,_|                                           
%    _____ _______ _    _    _         __   ____     ________    
%   / ____|__   __| |  | |  | |        \ \ / /\ \   / /___  /    
%  | (___    | |  | |__| |  | |_ ___    \ V /  \ \_/ /   / /     
%   \___ \   | |  |  __  |  | __/ _ \    > <    \   /   / /      
%   ____) |  | |  | |  | |  | || (_) |  / . \    | |   / /__     
%  |_____/   |_|  |_|  |_|   \__\___/  /_/ \_\   |_|  /_____|    
%   / ____|                            (_)                       
%  | |     ___  _ ____   _____ _ __ ___ _  ___  _ __  ___        
%  | |    / _ \| '_ \ \ / / _ \ '__/ __| |/ _ \| '_ \/ __|       
%  | |___| (_) | | | \ V /  __/ |  \__ \ | (_) | | | \__ \       
%   \_____\___/|_| |_|\_/ \___|_|  |___/_|\___/|_| |_|___/       
%                                                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                               
%% fcn_Path_convertXY2St
XY_points = [-2 -1; -1 0; -0.5 0.4; 0 0; 0.5 -0.5; 1 -0.4];
referencePath = [-3 -3; -1 -0.5; 0.5 0; 3 3];
flag_snap_type = 3;

St_points_XY = fcn_Path_convertXY2St(referencePath,XY_points, flag_snap_type);
St_points_ref = fcn_Path_convertXY2St(referencePath,referencePath, flag_snap_type);


fig_num = 999;
figure(fig_num);
clf;

subplot(1,2,1);
hold on;
grid on;
axis equal;

plot(XY_points(:,1),XY_points(:,2),'b.-','LineWidth',3,'MarkerSize',20)
plot(referencePath(:,1),referencePath(:,2),'r.-','LineWidth',3,'MarkerSize',20)
title('XY coordinates');

subplot(1,2,2);
hold on;
grid on;
axis equal;
plot(St_points_XY(:,1),St_points_XY(:,2),'b.-','LineWidth',3,'MarkerSize',20)
plot(St_points_ref(:,1),St_points_ref(:,2),'r.-','LineWidth',3,'MarkerSize',20)
title('St coordinates');


%% fcn_Path_convertSt2XY
St_points = [2 -1; 3 0; 3.5 0.4; 4 0; 4.5 -0.5; 5 -0.4];
referencePath = [-3 -3; -1 -0.5; 0.5 0; 3 3];
flag_snap_type = 1;

St_points_ref   = fcn_Path_convertXY2St(referencePath,referencePath, flag_snap_type);
XY_points_from_St = fcn_Path_convertSt2XY(referencePath,St_points, flag_snap_type);


fig_num = 999;
figure(fig_num);
clf;

subplot(1,2,1);
hold on;
grid on;
axis equal;

plot(St_points(:,1),St_points(:,2),'b.-','LineWidth',3,'MarkerSize',20)
plot(St_points_ref(:,1),St_points_ref(:,2),'r.-','LineWidth',3,'MarkerSize',20)
title('St coordinates');

subplot(1,2,2);
hold on;
grid on;
axis equal;
plot(XY_points_from_St(:,1),XY_points_from_St(:,2),'b.-','LineWidth',3,'MarkerSize',20)
plot(referencePath(:,1),referencePath(:,2),'r.-','LineWidth',3,'MarkerSize',20)
title('XY coordinates');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See http://patorjk.com/software/taag/#p=display&f=Big&t=Miscellaneous
%   __  __ _              _ _                                  
%  |  \/  (_)            | | |                                 
%  | \  / |_ ___  ___ ___| | | __ _ _ __   ___  ___  _   _ ___ 
%  | |\/| | / __|/ __/ _ \ | |/ _` | '_ \ / _ \/ _ \| | | / __|
%  | |  | | \__ \ (_|  __/ | | (_| | | | |  __/ (_) | |_| \__ \
%  |_|  |_|_|___/\___\___|_|_|\__,_|_| |_|\___|\___/ \__,_|___/
%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                            
                                                             




% fcn_Path_debugPrintStringToNCharacters(input_sequence,N)
% fcn_Path_debugPrintStringToNCharacters
% Given a string and an integer N representing the number of characters to
% keep or pad, creates a new string of exactly length N by cropping the
% string or padding it (to the right) with spaces.
%
% NOTE: this function is being deprecated. See the GitHub library that
% hosts debugging tools that includes a better version:
% https://github.com/ivsg-psu/Errata_Tutorials_DebugTools/wiki



%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§

%% function fcn_INTERNAL_clearUtilitiesFromPathAndFolders
function fcn_INTERNAL_clearUtilitiesFromPathAndFolders
% Clear out the variables
clear global flag* FLAG*
clear flag*
clear path

% Clear out any path directories under Utilities
if ispc
    path_dirs = regexp(path,'[;]','split');
elseif ismac
    path_dirs = regexp(path,'[:]','split');
elseif isunix
    path_dirs = regexp(path,'[;]','split');
else
    error('Unknown operating system. Unable to continue.');
end

utilities_dir = fullfile(pwd,filesep,'Utilities');
for ith_dir = 1:length(path_dirs)
    utility_flag = strfind(path_dirs{ith_dir},utilities_dir);
    if ~isempty(utility_flag)
        rmpath(path_dirs{ith_dir})
    end
end

% Delete the Utilities folder, to be extra clean!
if  exist(utilities_dir,'dir')
    [status,message,message_ID] = rmdir(utilities_dir,'s');
    if 0==status
        error('Unable remove directory: %s \nReason message: %s \nand message_ID: %s\n',utilities_dir, message,message_ID);
    end
end

end % Ends fcn_INTERNAL_clearUtilitiesFromPathAndFolders
%% fcn_INTERNAL_initializeUtilities
function  fcn_INTERNAL_initializeUtilities(library_name,library_folders,library_url,this_project_folders)
% Reset all flags for installs to empty
clear global FLAG*

fprintf(1,'Installing utilities necessary for code ...\n');

% Dependencies and Setup of the Code
% This code depends on several other libraries of codes that contain
% commonly used functions. We check to see if these libraries are installed
% into our "Utilities" folder, and if not, we install them and then set a
% flag to not install them again.

% Set up libraries
for ith_library = 1:length(library_name)
    dependency_name = library_name{ith_library};
    dependency_subfolders = library_folders{ith_library};
    dependency_url = library_url{ith_library};

    fprintf(1,'\tAdding library: %s ...',dependency_name);
    fcn_INTERNAL_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url);
    clear dependency_name dependency_subfolders dependency_url
    fprintf(1,'Done.\n');
end

% Set dependencies for this project specifically
fcn_DebugTools_addSubdirectoriesToPath(pwd,this_project_folders);

disp('Done setting up libraries, adding each to MATLAB path, and adding current repo folders to path.');
end % Ends fcn_INTERNAL_initializeUtilities


function fcn_INTERNAL_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url, varargin)
%% FCN_DEBUGTOOLS_INSTALLDEPENDENCIES - MATLAB package installer from URL
%
% FCN_DEBUGTOOLS_INSTALLDEPENDENCIES installs code packages that are
% specified by a URL pointing to a zip file into a default local subfolder,
% "Utilities", under the root folder. It also adds either the package
% subfoder or any specified sub-subfolders to the MATLAB path.
%
% If the Utilities folder does not exist, it is created.
% 
% If the specified code package folder and all subfolders already exist,
% the package is not installed. Otherwise, the folders are created as
% needed, and the package is installed.
% 
% If one does not wish to put these codes in different directories, the
% function can be easily modified with strings specifying the
% desired install location.
% 
% For path creation, if the "DebugTools" package is being installed, the
% code installs the package, then shifts temporarily into the package to
% complete the path definitions for MATLAB. If the DebugTools is not
% already installed, an error is thrown as these tools are needed for the
% path creation.
% 
% Finally, the code sets a global flag to indicate that the folders are
% initialized so that, in this session, if the code is called again the
% folders will not be installed. This global flag can be overwritten by an
% optional flag input.
%
% FORMAT:
%
%      fcn_DebugTools_installDependencies(...
%           dependency_name, ...
%           dependency_subfolders, ...
%           dependency_url)
%
% INPUTS:
%
%      dependency_name: the name given to the subfolder in the Utilities
%      directory for the package install
%
%      dependency_subfolders: in addition to the package subfoder, a list
%      of any specified sub-subfolders to the MATLAB path. Leave blank to
%      add only the package subfolder to the path. See the example below.
%
%      dependency_url: the URL pointing to the code package.
%
%      (OPTIONAL INPUTS)
%      flag_force_creation: if any value other than zero, forces the
%      install to occur even if the global flag is set.
%
% OUTPUTS:
%
%      (none)
%
% DEPENDENCIES:
%
%      This code will automatically get dependent files from the internet,
%      but of course this requires an internet connection. If the
%      DebugTools are being installed, it does not require any other
%      functions. But for other packages, it uses the following from the
%      DebugTools library: fcn_DebugTools_addSubdirectoriesToPath
%
% EXAMPLES:
%
% % Define the name of subfolder to be created in "Utilities" subfolder
% dependency_name = 'DebugTools_v2023_01_18';
%
% % Define sub-subfolders that are in the code package that also need to be
% % added to the MATLAB path after install; the package install subfolder
% % is NOT added to path. OR: Leave empty ({}) to only add 
% % the subfolder path without any sub-subfolder path additions. 
% dependency_subfolders = {'Functions','Data'};
%
% % Define a universal resource locator (URL) pointing to the zip file to
% % install. For example, here is the zip file location to the Debugtools
% % package on GitHub:
% dependency_url = 'https://github.com/ivsg-psu/Errata_Tutorials_DebugTools/blob/main/Releases/DebugTools_v2023_01_18.zip?raw=true';
%
% % Call the function to do the install
% fcn_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url)
%
% This function was written on 2023_01_23 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2023_01_23:
% -- wrote the code originally
% 2023_04_20:
% -- improved error handling
% -- fixes nested installs automatically

% TO DO
% -- Add input argument checking

flag_do_debug = 0; % Flag to show the results for debugging
flag_do_plots = 0; % % Flag to plot the final results
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
end


%% check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag_check_inputs
    % Are there the right number of inputs?
    narginchk(3,4);
end

%% Set the global variable - need this for input checking
% Create a variable name for our flag. Stylistically, global variables are
% usually all caps.
flag_varname = upper(cat(2,'flag_',dependency_name,'_Folders_Initialized'));

% Make the variable global
eval(sprintf('global %s',flag_varname));

if nargin==4
    if varargin{1}
        eval(sprintf('clear global %s',flag_varname));
    end
end

%% Main code starts here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if ~exist(flag_varname,'var') || isempty(eval(flag_varname))
    % Save the root directory, so we can get back to it after some of the
    % operations below. We use the Print Working Directory command (pwd) to
    % do this. Note: this command is from Unix/Linux world, but is so
    % useful that MATLAB made their own!
    root_directory_name = pwd;

    % Does the directory "Utilities" exist?
    utilities_folder_name = fullfile(root_directory_name,'Utilities');
    if ~exist(utilities_folder_name,'dir')
        % If we are in here, the directory does not exist. So create it
        % using mkdir
        [success_flag,error_message,message_ID] = mkdir(root_directory_name,'Utilities');

        % Did it work?
        if ~success_flag
            error('Unable to make the Utilities directory. Reason: %s with message ID: %s\n',error_message,message_ID);
        elseif ~isempty(error_message)
            warning('The Utilities directory was created, but with a warning: %s\n and message ID: %s\n(continuing)\n',error_message, message_ID);
        end

    end

    % Does the directory for the dependency folder exist?
    dependency_folder_name = fullfile(root_directory_name,'Utilities',dependency_name);
    if ~exist(dependency_folder_name,'dir')
        % If we are in here, the directory does not exist. So create it
        % using mkdir
        [success_flag,error_message,message_ID] = mkdir(utilities_folder_name,dependency_name);

        % Did it work?
        if ~success_flag
            error('Unable to make the dependency directory: %s. Reason: %s with message ID: %s\n',dependency_name, error_message,message_ID);
        elseif ~isempty(error_message)
            warning('The %s directory was created, but with a warning: %s\n and message ID: %s\n(continuing)\n',dependency_name, error_message, message_ID);
        end

    end

    % Do the subfolders exist?
    flag_allFoldersThere = 1;
    if isempty(dependency_subfolders{1})
        flag_allFoldersThere = 0;
    else
        for ith_folder = 1:length(dependency_subfolders)
            subfolder_name = dependency_subfolders{ith_folder};
            
            % Create the entire path
            subfunction_folder = fullfile(root_directory_name, 'Utilities', dependency_name,subfolder_name);
            
            % Check if the folder and file exists that is typically created when
            % unzipping.
            if ~exist(subfunction_folder,'dir')
                flag_allFoldersThere = 0;
            end
        end
    end

    % Do we need to unzip the files?
    if flag_allFoldersThere==0
        % Files do not exist yet - try unzipping them.
        save_file_name = tempname(root_directory_name);
        zip_file_name = websave(save_file_name,dependency_url);
        % CANT GET THIS TO WORK --> unzip(zip_file_url, debugTools_folder_name);

        % Is the file there?
        if ~exist(zip_file_name,'file')
            error(['The zip file: %s for dependency: %s did not download correctly.\n' ...
                'This is usually because permissions are restricted on ' ...
                'the current directory. Check the code install ' ...
                '(see README.md) and try again.\n'],zip_file_name, dependency_name);
        end

        % Try unzipping
        unzip(zip_file_name, dependency_folder_name);

        % Did this work? If so, directory should not be empty
        directory_contents = dir(dependency_folder_name);
        if isempty(directory_contents)
            error(['The necessary dependency: %s has an error in install ' ...
                'where the zip file downloaded correctly, ' ...
                'but the unzip operation did not put any content ' ...
                'into the correct folder. ' ...
                'This suggests a bad zip file or permissions error ' ...
                'on the local computer.\n'],dependency_name);
        end

        % Check if is a nested install (for example, installing a folder
        % "Toolsets" under a folder called "Toolsets"). This can be found
        % if there's a folder whose name contains the dependency_name
        flag_is_nested_install = 0;
        for ith_entry = 1:length(directory_contents)
            if contains(directory_contents(ith_entry).name,dependency_name)
                if directory_contents(ith_entry).isdir
                    flag_is_nested_install = 1;
                    install_directory_from = fullfile(directory_contents(ith_entry).folder,directory_contents(ith_entry).name);
                    install_files_from = fullfile(directory_contents(ith_entry).folder,directory_contents(ith_entry).name,'*'); % BUG FIX - For Macs, must be *, not *.*
                    install_location_to = fullfile(directory_contents(ith_entry).folder);
                end
            end
        end

        if flag_is_nested_install
            [status,message,message_ID] = movefile(install_files_from,install_location_to);
            if 0==status
                error(['Unable to move files from directory: %s\n ' ...
                    'To: %s \n' ...
                    'Reason message: %s\n' ...
                    'And message_ID: %s\n'],install_files_from,install_location_to, message,message_ID);
            end
            [status,message,message_ID] = rmdir(install_directory_from);
            if 0==status
                error(['Unable remove directory: %s \n' ...
                    'Reason message: %s \n' ...
                    'And message_ID: %s\n'],install_directory_from,message,message_ID);
            end
        end

        % Make sure the subfolders were created
        flag_allFoldersThere = 1;
        if ~isempty(dependency_subfolders{1})
            for ith_folder = 1:length(dependency_subfolders)
                subfolder_name = dependency_subfolders{ith_folder};
                
                % Create the entire path
                subfunction_folder = fullfile(root_directory_name, 'Utilities', dependency_name,subfolder_name);
                
                % Check if the folder and file exists that is typically created when
                % unzipping.
                if ~exist(subfunction_folder,'dir')
                    flag_allFoldersThere = 0;
                end
            end
        end
         % If any are not there, then throw an error
        if flag_allFoldersThere==0
            error(['The necessary dependency: %s has an error in install, ' ...
                'or error performing an unzip operation. The subfolders ' ...
                'requested by the code were not found after the unzip ' ...
                'operation. This suggests a bad zip file, or a permissions ' ...
                'error on the local computer, or that folders are ' ...
                'specified that are not present on the remote code ' ...
                'repository.\n'],dependency_name);
        else
            % Clean up the zip file
            delete(zip_file_name);
        end

    end


    % For path creation, if the "DebugTools" package is being installed, the
    % code installs the package, then shifts temporarily into the package to
    % complete the path definitions for MATLAB. If the DebugTools is not
    % already installed, an error is thrown as these tools are needed for the
    % path creation.
    %
    % In other words: DebugTools is a special case because folders not
    % added yet, and we use DebugTools for adding the other directories
    if strcmp(dependency_name(1:10),'DebugTools')
        debugTools_function_folder = fullfile(root_directory_name, 'Utilities', dependency_name,'Functions');

        % Move into the folder, run the function, and move back
        cd(debugTools_function_folder);
        fcn_DebugTools_addSubdirectoriesToPath(dependency_folder_name,dependency_subfolders);
        cd(root_directory_name);
    else
        try
            fcn_DebugTools_addSubdirectoriesToPath(dependency_folder_name,dependency_subfolders);
        catch
            error(['Package installer requires DebugTools package to be ' ...
                'installed first. Please install that before ' ...
                'installing this package']);
        end
    end


    % Finally, the code sets a global flag to indicate that the folders are
    % initialized.  Check this using a command "exist", which takes a
    % character string (the name inside the '' marks, and a type string -
    % in this case 'var') and checks if a variable ('var') exists in matlab
    % that has the same name as the string. The ~ in front of exist says to
    % do the opposite. So the following command basically means: if the
    % variable named 'flag_CodeX_Folders_Initialized' does NOT exist in the
    % workspace, run the code in the if statement. If we look at the bottom
    % of the if statement, we fill in that variable. That way, the next
    % time the code is run - assuming the if statement ran to the end -
    % this section of code will NOT be run twice.

    eval(sprintf('%s = 1;',flag_varname));
end


%% Plot the results (for debugging)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _
%  |  __ \     | |
%  | |  | | ___| |__  _   _  __ _
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_do_plots

    % Nothing to do!



end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends function fcn_DebugTools_installDependencies