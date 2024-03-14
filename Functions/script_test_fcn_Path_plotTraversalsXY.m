% script_test_fcn_Path_plotTraversalsXY.m
% Tests fcn_Path_plotTraversalsXY
       
% Revision history:
%      2020_11_10
%      -- first write of the code
%      2021_01_06
%      -- code clean-up to remove unneeded yaw function
%     2021_01_07
%     -- renamed function to show that traversals being used, not path
%     -- fixed small typo in function name 

close all





%% Call the plot command to show how it works. First, put it into our figure
% to show that it will auto-label the axes and create a new figure (NOT
% figure 11 here) to plot the data.

fig_num = 11;

% Fill in some dummy data
paths = fcn_Path_fillSamplePaths;
 

% Convert paths into traversals
for i_traveral = 1:length(paths)
    traversal = fcn_Path_convertPathToTraversalStructure(paths{i_traveral});
    data.traversal{i_traveral} = traversal;
end

figure(fig_num);
clf;
temp = figure(fig_num);

fcn_Path_plotTraversalsXY(data);

assert(isempty(temp.Children))

%% Next, specify the figure number to show that it will NOT auto-label the
% axes if figure is already given and it puts the plots into this figure.

fig_num = 11;

% Fill in some dummy data
paths = fcn_Path_fillSamplePaths;
 

% Convert paths into traversals
for i_traveral = 1:length(paths)
    traversal = fcn_Path_convertPathToTraversalStructure(paths{i_traveral});
    data.traversal{i_traveral} = traversal;
end

figure(fig_num);
clf;
temp = figure(fig_num);

fcn_Path_plotTraversalsXY(data,fig_num);

assert(~isempty(temp.Children))