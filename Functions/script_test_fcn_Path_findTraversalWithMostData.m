% script_test_fcn_Path_findTraversalWithMostData.m
% tests function: fcn_Path_findTraversalWithMostData.m

% Revision history
%      2020_11_10
%      -- first write of the code
%      2021_01_07
%      -- renamed function to change paths to traversals
%     2020_01_09
%     -- added more comments during clean-up


%% Create some dummy test paths

% Call the function that fills in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
for i_Path = 1:length(paths)
    traversal = fcn_Path_convertPathToTraversalStructure(paths{i_Path});
    data.traversal{i_Path} = traversal;
end

% Plot the results?
if 1==1
    fig_num = 12;
    fcn_Path_plotTraversalsYaw(data,fig_num);
    fig_num = 13;
    fcn_Path_plotTraversalsXY(data,fig_num);
end

%% Call the function

index_of_longest = fcn_Path_findTraversalWithMostData(data);
fprintf(1,'The longest path of the %.0d paths was path %.0d with %.0d elements\n',length(data.traversal),index_of_longest,length(data.traversal{index_of_longest}.X));



