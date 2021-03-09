% script_test_fcn_Path_findAverageTraversalViaOrthoProjection
% Tests the function: fcn_Path_findAverageTraversalViaOrthoProjection

% Revision history
%     2020_11_10
%     -- first write of the code
%     2021_01_07
%     -- lots of bug fixes as we demo for the team (lol)
%     2020_01_09
%     -- added more comments during clean-up

close all;
clc;

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
for i_Path = 1:length(paths_array)
    traversal = fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    data.traversal{i_Path} = traversal;
end

% Plot the results?
if 1==0
    fig_num = 12;
    fcn_Path_plotTraversalsYaw(data,fig_num);
    fig_num = 13;
    fcn_Path_plotTraversalsXY(data,fig_num);
end


%% Example 1: Basic call 
path_average_final3 = fcn_Path_findAverageTraversalViaOrthoProjection(data);

% Plot the final XY result of orthogonal
path_points_fig = 33333;
fcn_Path_plotTraversalsXY(data,path_points_fig);
hold on
plot(path_average_final3.X,path_average_final3.Y,'Linewidth',4);
title('Original paths and final average path via orthogonal projections')
xlabel('X [m]')
ylabel('Y [m]')


%% Show some extra features of the orthogonal projection averaging function
path_points_fig = 444444;
figure(path_points_fig);
clf;
hold on;

[path_average_final4, closestXs, closestYs, closestDistances]  = ...
    fcn_Path_findAverageTraversalViaOrthoProjection(data);

% Plot the original inputs? 
if 1==0
    fcn_Path_plotTraversalsXY(data,path_points_fig);
    hold on
end


% Show the "hit" points
for i_path = 1:length(closestXs(1,:))
    plot(closestXs(:,i_path),closestYs(:,i_path),'o-');
end

% Calculate the standard deviation
std_deviation = std(closestDistances,0,'all');

% Overlay the average with 3 times the standard deviations
fcn_Path_plotTraversalXYWithVarianceBands(path_average_final4,...
    3*std_deviation,path_points_fig);

num_outliers = length(find(abs(closestDistances)>3*std_deviation));
num_points = length(closestXs(:,1))*length(closestXs(1,:));

title(sprintf('There were %.0d points out of %.0d points that were more than 3 standard deviations from the average path, or %.2f percent',num_outliers,num_points,num_outliers/num_points));
xlabel('X [m]')
ylabel('Y [m]')

figure(44);
histogram(closestDistances,50);
title('Histogram of distances of paths from the final average path');
