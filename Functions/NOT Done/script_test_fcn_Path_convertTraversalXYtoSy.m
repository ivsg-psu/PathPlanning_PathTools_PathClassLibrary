% script_test_fcn_Path_convertTraversalXYtoSy.m
% Tests script_test_fcn_Path_convertTraversalXYtoSy

% Revision history:
%      2021_03_21
%      -- first write of the code (from
%      script_test_fcn_Path_findOrthoScatterFromTraversalToTraversals)


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
if 1==1
    %     fig_num = 12;
    %     fcn_Path_plotTraversalsYaw(all_traversals,fig_num);
    fig_num = 13;
    fcn_Path_plotTraversalsXY(all_traversals,fig_num);
    title('Original paths in XY');
end

%% Test case 1: basic call
reference_traversal = all_traversals.traversal{2};
reference_station_points = (0:30:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 40;
fig_num = 1;

[closestXs, closestYs, closestDistances] = ...
    fcn_Path_convertTraversalXYtoSy(...
    reference_station_points, reference_traversal, all_traversals,...
    flag_rounding_type,search_radius,fig_num); %#ok<*ASGLU>

%% Test case 2: basic call with different reference
reference_traversal = all_traversals.traversal{1};
reference_station_points = (0:30:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 40;
fig_num = 111;

[closestXs, closestYs, closestDistances] = ...
    fcn_Path_convertTraversalXYtoSy(...
    reference_station_points, reference_traversal, all_traversals,...
    flag_rounding_type,search_radius,fig_num); %#ok<*ASGLU>

reference_traversal = all_traversals.traversal{3};
reference_station_points = (0:10:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 40;
fig_num = 222;

[closestXs, closestYs, closestDistances] = ...
    fcn_Path_convertTraversalXYtoSy(...
    reference_station_points, reference_traversal, all_traversals,...
    flag_rounding_type,search_radius,fig_num); %#ok<*ASGLU>



%% Test case 3: basic call 1 with finer resolution
reference_traversal = all_traversals.traversal{2};
reference_station_points = (0:10:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 40;
fig_num = 2;

[closestXs, closestYs, closestDistances] = ...
    fcn_Path_convertTraversalXYtoSy(...
    reference_station_points, reference_traversal, all_traversals,...
    flag_rounding_type,search_radius,fig_num);

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

%% Another advanced test, with crossigs
% Set up data
% close all
reference_path = [0 0; 1 2; 2 0];
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