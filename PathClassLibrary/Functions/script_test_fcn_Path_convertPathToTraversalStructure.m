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

close all
clc

% Fill in some dummy data
paths_array = fcn_Path_fillSamplePaths;

%% Example 1 - show how it works
% Basic call with a figure option to plot output repeatedly onto figure
fig_num = 1;
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

figure(fig_num);
clf;
hold on;
grid on;
grid minor;

Station_step = 40;

fcn_Path_plotTraversalsXY(data,fig_num);
xlabel('X [m]');
ylabel('Y [m]');

for i_traveral = 1:length(paths)
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
