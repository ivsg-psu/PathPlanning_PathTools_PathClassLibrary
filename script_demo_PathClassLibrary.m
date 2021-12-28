% This is a demonstration script to show the primary functionality of the
% path class library.

%% Set up workspace
if ~exist('flag_paths_were_added_already','var')
    
    clc
    close all
    
    % add necessary directories for functions recursively
    addpath(genpath([pwd, filesep, 'Functions']))
    
    % % add necessary directories for data
    % addpath([pwd '/Data'])  % This is where sample maps are stored
    
    % add necessary directories for Utilities to the path
    if(exist([pwd, filesep,  'Utilities'],'dir'))
        addpath(genpath([pwd, ilesep, 'Utilities']))  % This is where GPS utilities are stored
    else
        error('No Utilities directory exists to be added to the path. Please create one (see README.md) and run again.');
    end
    
    % set a flag so we do not have to do this again
    flag_paths_were_added_already = 1;
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

%% Show how input arguments are checked, fcn_Path_checkInputsToFunctions
% TO-DO - move debug tools into utilities and remove the checkinputs capability out of Path class
path_test = [4 1; 2 1];
fcn_Path_checkInputsToFunctions(path_test, 'path');

%% Show how to load sample paths, fcn_Path_fillSamplePaths  
% Call the function to fill in an array of "path" type
paths_array = fcn_Path_fillSamplePaths;

% We can even save one of these as a single "path"
single_path = paths_array{1};

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
yaw_angles = fcn_Path_calcYawFromPathSegments(path_to_check,fig_num); %#ok<*NASGU>


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

%% Show how to snap a point onto a path, fcn_Path_snapPointOntoNearestPath
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
    fcn_Path_snapPointOntoNearestPath(point, pathXY,fignum);
 
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
% finds the location on a nearby traversal that is closest to the central
% traveral at each station point. Closest is defined via an orthogonal
% projection (or modifications of orthogonal projections) from the central
% traversal outward toward nearby traversals.


% Set up data
close all
central_path = [0 0; 1 1; 2 0];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [-1 0.5; 0.5 2; 1.5 2; 3 0.5];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);
stations = [0; 1; 2^0.5-0.1; 2^0.5; 2^0.5+.1; 2; central_traversal.Station(end)];
search_radius = 20; % Distance to search for nearby segments

% AVERAGING example 1 - default setting
fig_num = 1;
flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via prior segment (default, flag=1)');

% AVERAGING example 2 - use following segment
fig_num = 2;
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via following segment (flag=2)');

% AVERAGING example 3 - use average of both segments
fig_num = 3;
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via averaging prior and following segment at vertex (flag=3)');

% AVERAGING example 4 - use average always
fig_num = 4;
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via averaging everywhere (flag=4)');

%% Find closest distance of many traversals to one central taversal, fcn_Path_findOrthoScatterFromTraversalToTraversals
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
flag_rounding_type = 3;
search_radius = 40;
fig_num = 1;

[closestXs, closestYs, closestDistances] = ...
       fcn_Path_findOrthoScatterFromTraversalToTraversals(...
       reference_station_points, reference_traversal, all_traversals,...
       flag_rounding_type,search_radius,fig_num); %#ok<*ASGLU>
   
figure(11);
histogram([closestDistances(:,1);closestDistances(:,3)],30);
title('Histogram of all orthogonal distance projections');
