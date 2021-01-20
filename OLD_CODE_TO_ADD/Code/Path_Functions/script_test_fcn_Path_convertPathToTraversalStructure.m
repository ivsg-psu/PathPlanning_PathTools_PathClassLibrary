% script_test_fcn_Path_convertPathToTraversalStructure.m
% tests fcn_Path_convertPathToTraversalStructure.m

% Revision history:
%      2020_11_12 
%      -- first wrote the code
%      2021_01_06
%      -- a bit more comments, added more plotting, renamed
%      2021_01_07
%      -- changed name to reflect that the input is a Path

close all
clc

% Fill in some dummy data
paths = fcn_Path_fillSamplePaths;


% Basic call
fig_num = 1;
for i_traveral = 1:length(paths)
    traversal = fcn_Path_convertPathToTraversalStructure(paths{i_traveral},fig_num);
    data.traversal{i_traveral} = traversal;
end


% Call the plot command to show ALL results in XY
fig_num = 12;
fcn_Path_plotTraversalsXY(data,fig_num);


