% script_test_fcn_fcn_Path_calcYawFromPathSegments
% Tests fcn_Path_calcYawFromPathSegments
       
% Revision history:
%      2021_01_06
%      -- first write of the code
%      2021_01_07
%      -- updated function calls to reflect paths vs traversals
%      -- minor comment clean-ups

close all

% Clear any old variables
clear all_traversals

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Pick first path as reference_traversal structure
path_to_check = paths_array{1};

% Pick first path as reference_traversal structure
traversal_to_check = fcn_Path_convertPathToTraversalStructure(paths_array{1});
all_traversals.traversal{1} = traversal_to_check;


% Plot the results? (Note: they are plotted below as well)
if 1==1
    fig_num = 3;
    fcn_Path_plotTraversalsYaw(all_traversals,fig_num);
    fig_num = 2;
    fcn_Path_plotTraversalsXY(all_traversals,fig_num);
end

%% Test case 1: Basic call with one path
fig_num = 1;
yaw_angles = fcn_Path_calcYawFromPathSegments(path_to_check,fig_num); %#ok<*NASGU>


%% Test case 2: Multiple paths
close all;

for i_path = 1:length(paths_array)
    % Pick first path as reference_traversal structure
    path_to_check = paths_array{i_path};
    
    % Pick first path as reference_traversal structure
    traversal_to_check = fcn_Path_convertPathToTraversalStructure(paths_array{i_path});
    all_traversals.traversal{i_path} = traversal_to_check;
    
    fig_num = 2222;
    yaw_angles = fcn_Path_calcYawFromPathSegments(path_to_check,fig_num);


end

% Plot the results? (Note: they are plotted below as well)
if 1==1
    fig_num = 3;
    fcn_Path_plotTraversalsYaw(all_traversals,fig_num);
    fig_num = 2;
    fcn_Path_plotTraversalsXY(all_traversals,fig_num);
end

%% Test case 3: basic call with degenerate path
fig_num = 3333;
path_to_check = [1 1; 0 0];
yaw_angles = fcn_Path_calcYawFromPathSegments(path_to_check,fig_num);
