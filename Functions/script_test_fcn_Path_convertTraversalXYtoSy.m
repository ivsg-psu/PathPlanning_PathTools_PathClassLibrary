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
reference_station_points = (0:10:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 40;
fig_num = 1;

[closestXs, closestYs, closestDistances] = ...
       fcn_Path_convertTraversalXYtoSy(...
       reference_station_points, reference_traversal, all_traversals,...
       flag_rounding_type,search_radius,fig_num); %#ok<*ASGLU>
   

%% Test case 2: basic call with finer resolution
reference_traversal = all_traversals.traversal{2};
reference_station_points = (0:1:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 40;
fig_num = 2;

[closestXs, closestYs, closestDistances] = ...
       fcn_Path_convertTraversalXYtoSy(...
       reference_station_points, reference_traversal, all_traversals,...
       flag_rounding_type,search_radius,fig_num);
   