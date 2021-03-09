% script_test_fcn_Path_findAverageTraversalViaStationAlignment
% tests the function: script_test_fcn_Path_findAverageTraversalViaStationAlignment

% Revision history
%      2020_11_10
%      -- first write of the code
%      2021_01_07
%      -- lots of bug fixes as we demo for the team (lol)

close all;
clc;

% Fill in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
for i_Path = 1:length(paths)
    traversal = fcn_Path_convertPathToTraversalStructure(paths{i_Path});
    data.traversal{i_Path} = traversal;
end

% Plot the results?
if 1==0
    fig_num = 12;
    fcn_Path_plotTraversalsYaw(data,fig_num);
    fig_num = 13;
    fcn_Path_plotTraversalsXY(data,fig_num);
end

%% Call the station-averaging function

[aligned_Data_ByStation,mean_Data] = ...
    fcn_Path_findAverageTraversalViaStationAlignment(data);

% Plot the final XY result of mean station
path_points_fig = 11111;
fcn_Path_plotTraversalsXY(data,path_points_fig);
hold on
plot(mean_Data.mean_xEast,mean_Data.mean_yNorth,'Linewidth',4);
title('Original paths and final average path via station averaging')
xlabel('X [m]')
ylabel('Y [m]')


