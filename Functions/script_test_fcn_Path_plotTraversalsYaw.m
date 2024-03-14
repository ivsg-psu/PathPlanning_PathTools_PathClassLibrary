% script_test_fcn_Path_plotTraversalsYaw.m
% tests fcn_Path_plotTraversalsYaw

% Revision history:
% 2020_11_12 - S. Brennan
% -- first wrote the code
% 
% 2021_01_06 - S. Brennan
% -- a bit more comments, renamed function to be traversals, not path
% 
% 2021_01_07 - S. Brennan
% -- fixed typos in arguments for function calls, based on recent edits 
% 
% 2024_03_14 - S. Brennan
% -- added clear assertions, moved data loading into test sections

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

fcn_Path_plotTraversalsYaw(data);
assert(isempty(temp.Children))

%% Next, specify the figure number to show that it will NOT auto-label the
% axes if figure is already given and it puts the plots into this figure.
fig_num = 12;

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

fcn_Path_plotTraversalsYaw(data, fig_num);
assert(~isempty(temp.Children))
