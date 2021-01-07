% script_test_fcn_Path_fillSamplePaths.m
% tests fcn_Path_fillSamplePaths.m

% Revision history
%     2020_11_10
%     -- first write of the code
%     2021_01_07
%     -- cleanup to include trajectories functions

clear all
clc

% Call the function to fill in an array of "path" type
path_array = fcn_Path_fillSamplePaths;

% We can even save one of these as a single "path"
single_path = path_array{1};

% Convert paths to traversals structures. Each traversal instance is a
% "traversal" type, and the array called "data" below is a "traversals"
% type.
for i_Path = 1:length(path_array)
    traversal = fcn_Path_convertPathToTraversalStructure(path_array{i_Path});
    data.traversal{i_Path} = traversal;
end


% Call the plot command to show results in XY
fig_num = 12;
fcn_Path_plotTraversalsXY(data,fig_num);


