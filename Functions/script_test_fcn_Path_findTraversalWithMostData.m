% script_test_fcn_Path_findTraversalWithMostData.m


%% Create some dummy test paths

% Call the function that fills in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
for i_Path = 1:length(paths)
    traversal = fcn_Path_convertXYtoTraversalStructure(paths{i_Path}(:,1),paths{i_Path}(:,2));
    data.traversal{i_Path} = traversal;
end

% Plot the results?
if 1==1
    fig_num = 12;
    fcn_Path_plotPathYaw(data,fig_num);
    fig_num = 13;
    fcn_Path_plotPathXY(data,fig_num);
end

%% Call the function

index_of_longest = fcn_Path_findTraversalWithMostData(data);
fprintf(1,'The longest path of the %.0d paths was path %.0d with %.0d elements\n',length(data.traversal),index_of_longest,length(data.traversal{index_of_longest}.X));



