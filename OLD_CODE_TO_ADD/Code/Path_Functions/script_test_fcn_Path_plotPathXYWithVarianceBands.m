% script_test_fcn_Path_plotPathXYWithVarianceBands.m
% Tests fcn_Path_plotPathXYWithVarianceBands
       
% Revision history:
% 2021_01_05
% -- first write of the code

close all
clc

% Clear any old variables
clear all_traversals

% Fill in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertXYtoTraversalStructure(paths{1}(:,1),paths{1}(:,2));
all_traversals.traversal{1} = reference_traversal;


%% Test case 1: basic call for one trajectory
fcn_Path_plotPathXYWithVarianceBands(reference_traversal);


%% Test case 2: advanced call for one trajectory - specify figure
fig_num = 22;
std_deviation = [];
fcn_Path_plotPathXYWithVarianceBands(reference_traversal,...
    std_deviation,fig_num);

%% Test case 3: advanced call for one trajectory - specify std_deviation
fig_num = 31;
std_deviation = 1;
fcn_Path_plotPathXYWithVarianceBands(reference_traversal,...
    std_deviation,fig_num);
title(sprintf('Standard deviation: %.0d meters',std_deviation));

fig_num = 32;
std_deviation = 2;
fcn_Path_plotPathXYWithVarianceBands(reference_traversal,...
    std_deviation,fig_num);
title(sprintf('Standard deviation: %.0d meters',std_deviation));

fig_num = 35;
std_deviation = 5;
fcn_Path_plotPathXYWithVarianceBands(reference_traversal,...
    std_deviation,fig_num);
title(sprintf('Standard deviation: %.0d meters',std_deviation));

%% Test case 4: advanced call for multiple trajectories
fig_num = 4;
std_deviation = 2;
for i_Path = 1:length(paths)
    reference_traversal = fcn_Path_convertXYtoTraversalStructure(paths{i_Path}(:,1),paths{i_Path}(:,2));
    fcn_Path_plotPathXYWithVarianceBands(reference_traversal,...
        std_deviation,fig_num);
end