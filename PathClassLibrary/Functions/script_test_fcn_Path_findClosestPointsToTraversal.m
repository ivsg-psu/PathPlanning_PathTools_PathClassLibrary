% script_test_fcn_Path_findClosestPointsToTraversal.m
% This is a script to exercise the function: fcn_Path_findClosestPointsToTraversal.m
% This function was written on 2020_11_15 by S. Brennan
%     Modified on 2020_11_15 to prep for Path class
% Questions or comments? sbrennan@psu.edu

% Here is the function format:
% function [closestXs,closestYs,closestZs,closestYaws] = ...
% fcn_Path_findClosestPointsToTraversal(...
% reference_traversal,data,flag_yaw,flag_3D)

% Revision history
%     2021_01_09 
%     -- added more comments



%% BASIC examples
close all;
clear data;

% Create a dummy central path and convert it to a traversal
central_path = [0 0; 2 0];  
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [0 4; 2 4];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);
data.traversal{1} = nearby_traversal;
flag_yaw = 0;
flag_3D = 0;


% BASIC example 1 - parallel lines, query is in middle area
fig_num = 1;
% Calculate the closest point and distance on the nearby path
[closestXs,closestYs,closestZs] = ...
    fcn_Path_findClosestPointsToTraversal(central_traversal,data,flag_yaw,flag_3D,fig_num); %#ok<*ASGLU>


% BASIC example 2 - angled line segment adjacent to endpoint query
fig_num = 2;
nearby_path = [0 4; 2 7];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);
data.traversal{1} = nearby_traversal;
[closestXs,closestYs,closestZs] = ...
    fcn_Path_findClosestPointsToTraversal(central_traversal,data,flag_yaw,flag_3D,fig_num);


% BASIC example 3 - angled line segment adjacent to endpoint query but near-miss
fig_num = 3;
central_path = [0 0; 10 0];
central_traversal = ...
    fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [0 4; 10 7];
nearby_traversal = ...
    fcn_Path_convertPathToTraversalStructure(nearby_path);
data.traversal{1} = nearby_traversal;
[closestXs,closestYs,closestZs] = ...
    fcn_Path_findClosestPointsToTraversal(central_traversal,data,flag_yaw,flag_3D,fig_num);

% BASIC example 4 - angled line segment adjacent to startpoint query
fig_num = 4;
central_path = [0 0; 10 0];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [-1 4; 12 7];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);
data.traversal{1} = nearby_traversal;
[closestXs,closestYs,closestZs] = ...
    fcn_Path_findClosestPointsToTraversal(central_traversal,data,flag_yaw,flag_3D,fig_num);


% BASIC example 5 - parallel line segment adjacent to startpoint query but near-miss
fig_num = 5;
central_path = [0 0; 10 0];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [0 4; 10 4];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);
data.traversal{1} = nearby_traversal;
[closestXs,closestYs,closestZs] = ...
    fcn_Path_findClosestPointsToTraversal(central_traversal,data,flag_yaw,flag_3D,fig_num);


%% AVERAGING examples

% Set up data
close all
central_path = [0 0; 1 1; 2 0];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [-1 0.5; 0.5 2; 1.5 2; 3 0.5];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);
data.traversal{1} = nearby_traversal;

% AVERAGING example 1 - default setting
fig_num = 1;
[closestXs,closestYs,closestZs] = ...
    fcn_Path_findClosestPointsToTraversal(central_traversal,data,flag_yaw,flag_3D,fig_num);


%% NEGATIVE examples

% Prep the example and workspace
close all;
central_path = [-2 1; 1 4; 3 2];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [-1 0.5; 0.5 2; 1.5 2; 3 0.5];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);
data.traversal{1} = nearby_traversal;


% NEGATIVE example 1 - default setting
fig_num = 1;
[closestXs,closestYs,closestZs] = ...
    fcn_Path_findClosestPointsToTraversal(central_traversal,data,flag_yaw,flag_3D,fig_num);



%% MULTICROSS examples
close all;

% Setup
central_path = [-2 1; 1 4; 3 2; 5 2; 6 3; 7 2];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [-1 0.5; 0.5 2; 1.5 2; 3 0.5; 4 3; 5 4; 6 3; 7 1];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);
data.traversal{1} = nearby_traversal;

% MULTICROSS example 1 - default setting
fig_num = 1;
[closestXs,closestYs,closestZs] = ...
    fcn_Path_findClosestPointsToTraversal(central_traversal,data,flag_yaw,flag_3D,fig_num);



%% Real path examples
close all;

% Fill in some dummy data
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
for i_Path = 1:length(paths_array)
    traversal = fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    data.traversal{i_Path} = traversal;
end

% Call the plot command to show results in XY
fig_num = 12;
fcn_Path_plotTraversalsXY(data,fig_num);


central_traversal = data.traversal{1};
fig_num = 1;
[closestXs,closestYs,closestZs] = ...
    fcn_Path_findClosestPointsToTraversal(central_traversal,data,flag_yaw,flag_3D,fig_num);



% function print_results(stations,closest_path_point,distances)
% fprintf(1,'\n\nStation \t Location X \t Location Y \t Distance \n');
% for i_station =1:length(stations)
%     fprintf(1,'%.2f \t\t %.2f \t\t\t %.2f \t\t\t %.2f\n',stations(i_station),closest_path_point(i_station,1),closest_path_point(i_station,2),distances(i_station));
% end
% end
