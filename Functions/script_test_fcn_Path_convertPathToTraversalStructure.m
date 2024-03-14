% script_test_fcn_Path_convertPathToTraversalStructure.m
% tests fcn_Path_convertPathToTraversalStructure.m

% Revision history:
%      2020_11_12 
%      -- first wrote the code
%      2021_01_06
%      -- a bit more comments, added more plotting, renamed
%      2021_01_07
%      -- changed name to reflect that the input is a Path, and that we are
%      using an array of paths
%      -- more descriptive comments
%      2021_01_09
%      -- added another example to illustrate station errors
%      2021_03_21
%      -- added examples to illustrate 2D and 3D basic plotting

close all

%% Example 1.1 - show how it works for 2D
% Basic call 
fig_num = 11;
path_simple_2D = [1 1; 1 2; 3 4; 4 5; 7 7];
traversal = fcn_Path_convertPathToTraversalStructure(path_simple_2D,fig_num);
simple_example.traversal{1} = traversal;

fcn_Path_plotTraversalsXY(simple_example,fig_num);
xlabel('X [m]');
ylabel('Y [m]');

%% Example 1.2 - show how it works for 3D
% Basic call 
fig_num = 12;
path_simple_3D = [1 1 1; 1 2 2; 3 4 3; 4 5 3; 7 7 4];
traversal = fcn_Path_convertPathToTraversalStructure(path_simple_3D,fig_num);
simple_example.traversal{1} = traversal;

fcn_Path_plotTraversalsXY(simple_example,fig_num);
xlabel('X [m]');
ylabel('Y [m]');

%% Example 1.3 - show how it works
% Fill in some dummy data
paths_array = fcn_Path_fillSamplePaths;

% Basic call with a figure option to plot output repeatedly onto figure
fig_num = 13;
for i_traveral = 1:length(paths_array)
    traversal = fcn_Path_convertPathToTraversalStructure(paths_array{i_traveral},fig_num);
    data.traversal{i_traveral} = traversal;
end


% Call the plot command to show ALL results in XY, which gives same result
fig_num = 12;
fcn_Path_plotTraversalsXY(data,fig_num);
xlabel('X [m]');
ylabel('Y [m]');

%% Example 2: Show station discrepancies
% Plot station markers

fig_num = 2;

paths_array = fcn_Path_fillSamplePaths;

% Basic call with a figure option to plot output repeatedly onto figure
fig_num = 13;
for i_traveral = 1:length(paths_array)
    traversal = fcn_Path_convertPathToTraversalStructure(paths_array{i_traveral},fig_num);
    data.traversal{i_traveral} = traversal;
end

figure(fig_num);
clf;
hold on;
grid on;
grid minor;

Station_step = 40;

fcn_Path_plotTraversalsXY(data,fig_num);
xlabel('X [m]');
ylabel('Y [m]');

for i_traveral = 1:length(data.traversal)
    traversal_stations = data.traversal{i_traveral}.Station;
    for i_station = Station_step:Station_step:traversal_stations(end)
        index = find(traversal_stations >= i_station,1);
        plot(data.traversal{i_traveral}.X(index),...
            data.traversal{i_traveral}.Y(index),...
            'k.','Markersize',15);
        text(data.traversal{i_traveral}.X(index),...
            data.traversal{i_traveral}.Y(index),...
            sprintf('Station: %.2f',...
            data.traversal{i_traveral}.Station(index)));
    end
end
