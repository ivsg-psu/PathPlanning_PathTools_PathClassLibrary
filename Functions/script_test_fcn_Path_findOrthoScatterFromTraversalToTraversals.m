% script_test_fcn_Path_findOrthoScatterFromTraversalToTraversals.m
% Tests fcn_Path_findOrthoScatterFromTraversalToTraversals
       
% Revision history:
%      2021_01_02
%      -- first write of the code
%      2021_01_07:
%      -- renamed function to clarify paths versus traversals
%     2020_01_09
%     -- added more comments during clean-up

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
    fig_num = 12;
    fcn_Path_plotTraversalsYaw(all_traversals,fig_num);
    fig_num = 13;
    fcn_Path_plotTraversalsXY(all_traversals,fig_num);
end

%% Test case 1: basic call
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


%% Test case 2: basic call with finer resolution
reference_traversal = all_traversals.traversal{2};
reference_station_points = (0:1:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 40;
fig_num = 2;

[closestXs, closestYs, closestDistances] = ...
       fcn_Path_findOrthoScatterFromTraversalToTraversals(...
       reference_station_points, reference_traversal, all_traversals,...
       flag_rounding_type,search_radius,fig_num);
   
figure(22);
histogram([closestDistances(:,1);closestDistances(:,3)],30);
title('Histogram of all orthogonal distance projections');
