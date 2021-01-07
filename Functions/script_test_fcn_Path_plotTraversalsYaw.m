% script_test_fcn_Path_plotTraversalsYaw.m
% tests fcn_Path_plotTraversalsYaw

% Revision history:
%      2020_11_12 
%      -- first wrote the code
%      2021_01_06
%      -- a bit more comments, renamed function

clc
close all

% Fill in some dummy data
paths = fcn_Path_fillSamplePaths;
 

% Convert paths into traversals
for i_traveral = 1:length(paths)
    traversal = fcn_Path_convertXYtoTraversalStructure(paths{i_traveral}(:,1),paths{i_traveral}(:,2));
    data.traversal{i_traveral} = traversal;
end


%% Call the plot command to show how it works. First, put it into our figure
% to show that it will auto-label the axes and create a new figure (NOT
% figure 11 here) to plot the data.
figure(11);
fcn_Path_plotTraversalsYaw(data);


%% Next, specify the figure number to show that it will NOT auto-label the
% axes if figure is already given and it puts the plots into this figure.
fig_num = 12;
fcn_Path_plotTraversalsYaw(data,fig_num);

