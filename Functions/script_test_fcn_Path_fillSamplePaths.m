% script_test_fcn_Path_fillSamplePaths.m
% tests fcn_Path_fillSamplePaths.m

% Revision history
%     2020_11_10
%     -- first write of the code
%     2021_01_07
%     -- cleanup to include trajectories functions

close all
clc

% Call the function to fill in an array of "path" type
paths_array = fcn_Path_fillSamplePaths;

% We can even save one of these as a single "path"
single_path = paths_array{1};

% Convert paths to traversals structures. Each traversal instance is a
% "traversal" type, and the array called "data" below is a "traversals"
% type.
for i_Path = 1:length(paths_array)
    traversal = fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    data.traversal{i_Path} = traversal;
end


% Call the plot command to show results in XY
fig_num = 12;
fcn_Path_plotTraversalsXY(data,fig_num);

%% Test specific paths
fig_num = 22;

clear data
single_path = fcn_Path_fillSamplePaths(2);
traversal = fcn_Path_convertPathToTraversalStructure(single_path);
data.traversal{1} = traversal;
fcn_Path_plotTraversalsXY(data,fig_num);

%% Test specific paths
% This is a test path for testing the conversions from XY to Sy
fig_num = 33;

clear data
single_path = fcn_Path_fillSamplePaths(4);
traversal = fcn_Path_convertPathToTraversalStructure(single_path);
data.traversal{1} = traversal;
fcn_Path_plotTraversalsXY(data,fig_num);
