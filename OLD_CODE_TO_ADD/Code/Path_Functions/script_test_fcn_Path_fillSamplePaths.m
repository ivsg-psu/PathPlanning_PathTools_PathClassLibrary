% script_test_fcn_Path_fillSamplePaths.m
% tests fcn_Path_fillSamplePaths.m

% Revision history
%     2020_11_10
%     -- first write of the code
%     2021_01_07
%     -- cleanup to include trajectories functions

clear all
clc

% Call the function
paths = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
for i_Path = 1:length(paths)
    traversal = fcn_Path_convertPathToTraversalStructure(paths{i_Path});
    data.traversal{i_Path} = traversal;
end


% Call the plot command to show results in XY
fig_num = 12;
fcn_Path_plotTraversalsXY(data,fig_num);


