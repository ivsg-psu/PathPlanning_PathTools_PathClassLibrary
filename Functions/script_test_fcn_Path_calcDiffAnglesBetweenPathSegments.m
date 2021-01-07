% script_test_fcn_Path_calcDiffAnglesBetweenPathSegments.m
% Tests fcn_Path_calcDiffAnglesBetweenPathSegments
       
% Revision history:
% 2021_01_03
% -- first write of the code

close all
clc

% Clear any old variables
clear all_traversals

% Fill in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;

% Pick first path as reference_traversal structure
path_to_check = paths{1};

% Pick first path as reference_traversal structure
traversal_to_check = fcn_Path_convertPathToTraversalStructure(paths{1});
all_traversals.traversal{1} = traversal_to_check;


% Plot the results? (Note: they are plotted below as well)
if 1==1
    fig_num = 3;
    fcn_Path_plotTraversalsYaw(all_traversals,fig_num);
    fig_num = 2;
    fcn_Path_plotTraversalsXY(all_traversals,fig_num);
end

%% Test case 1: basic call
fig_num = 1;
diff_angles = fcn_Path_calcDiffAnglesBetweenPathSegments(path_to_check,fig_num);
