% script_test_fcn_Path_findAverageTraversalViaOrthoProjection
% Tests the function: fcn_Path_findAverageTraversalViaOrthoProjection

% Revision history
%     2020_11_10
%     -- first write of the code
%     2021_01_07
%     -- lots of bug fixes as we demo for the team (lol)
%     2020_01_09
%     -- added more comments during clean-up
%     2021_12_27:
%     -- corrected dependencies in comments
%     2022_01_03:
%     -- corrected typos in comments
%     -- fixed a bug where the Z value is not defined in loop
%     2022_01_06:
%     -- refactored code, added weighted averaging to prevent iteration
%     bouncing

%      [traversal_average, closestXs, closestYs, closestDistances] = ...
%      fcn_Path_findAverageTraversalViaOrthoProjection(...
%            data,
%            (reference_traversal),...
%            (num_iterations),
%            (weight_for_averaging),
%            (fig_num));

close all;




%% Example 1: Basic call 

clear data

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

path_average_final3 = fcn_Path_findAverageTraversalViaOrthoProjection(data);

% Plot the final XY result of orthogonal
path_points_fig = 33333;
fcn_Path_plotTraversalsXY(data,path_points_fig);
hold on
plot(path_average_final3.X,path_average_final3.Y,'Linewidth',4);
title('Original paths and final average path via orthogonal projections')
xlabel('X [m]')
ylabel('Y [m]')


%% BASIC CALL - Show some extra features of the orthogonal projection averaging function
path_points_fig = 444444;
figure(path_points_fig);
clf;
hold on;


clear data

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

[path_average_final4, closestXs, closestYs, closestDistances]  = ...
    fcn_Path_findAverageTraversalViaOrthoProjection(data);



% Show the "hit" points
for i_path = 1:length(closestXs(1,:))
    plot(closestXs(:,i_path),closestYs(:,i_path),'o-');
end

% Convert closestDistances into a single column
closestDistancesColumn = reshape(closestDistances,[length(closestDistances(:,1))*length(closestDistances(1,:)) 1]);

% Calculate the standard deviation
std_deviation = std(closestDistancesColumn,'omitnan');

% Overlay the average with 3 times the standard deviations
fcn_Path_plotTraversalXYWithVarianceBands(path_average_final4,...
    3*std_deviation,path_points_fig);

num_outliers = length(find(abs(closestDistancesColumn)>3*std_deviation));
num_points = length(closestXs(:,1))*length(closestXs(1,:));

title(sprintf('There were %.0d points out of %.0d points that were more than 3 standard deviations from the average path, or %.2f percent',num_outliers,num_points,num_outliers/num_points));
xlabel('X [m]')
ylabel('Y [m]')

figure(44);
histogram(closestDistancesColumn,50);
title('Histogram of distances of paths from the final average path');


%% Example 2: Basic call showing effect of reference traversal


clear data

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

% choose first one for traversal
reference_traversal = data.traversal{1};
path_average_final3 = fcn_Path_findAverageTraversalViaOrthoProjection(data,reference_traversal);

path_points_fig = 3331;
fcn_Path_plotTraversalsXY(data,path_points_fig);
hold on
plot(path_average_final3.X,path_average_final3.Y,'Linewidth',4);
title('Final average using first trajectory for reference')
xlabel('X [m]')
ylabel('Y [m]')

% choose second one for traversal
reference_traversal = data.traversal{2};
path_average_final3 = fcn_Path_findAverageTraversalViaOrthoProjection(data,reference_traversal);

path_points_fig = 3332;
fcn_Path_plotTraversalsXY(data,path_points_fig);
hold on
plot(path_average_final3.X,path_average_final3.Y,'Linewidth',4);
title('Final average using second trajectory for reference')
xlabel('X [m]')
ylabel('Y [m]')

% choose second one for traversal
reference_traversal = data.traversal{3};
path_average_final3 = fcn_Path_findAverageTraversalViaOrthoProjection(data,reference_traversal);

path_points_fig = 3333;
fcn_Path_plotTraversalsXY(data,path_points_fig);
hold on
plot(path_average_final3.X,path_average_final3.Y,'Linewidth',4);
title('Final average using third trajectory for reference')
xlabel('X [m]')
ylabel('Y [m]')


% choose a "dumb" one for traversal
%dumb_traversal = fcn_Path_convertPathToTraversalStructure([10 10; 15 65; 10 90]);
dumb_traversal = fcn_Path_convertPathToTraversalStructure([0 0; 40 40; 80 80]);
path_average_final3 = fcn_Path_findAverageTraversalViaOrthoProjection(data,dumb_traversal);

path_points_fig = 3334;
fcn_Path_plotTraversalsXY(data,path_points_fig);
hold on
plot(dumb_traversal.X,dumb_traversal.Y,'Linewidth',4);
plot(path_average_final3.X,path_average_final3.Y,'Linewidth',4);
title('Final average using dumb trajectory for reference')
legend('trajectory 1','trajectory 2','trajectory 3','dumb traversal','average')
xlabel('X [m]')
ylabel('Y [m]')


%% Example 2: Basic call showing effect of num_iterations


clear data

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

reference_traversal = data.traversal{3};

% Show the effect of one iteration
num_iterations = 1;
path_average_final3 = fcn_Path_findAverageTraversalViaOrthoProjection(data,reference_traversal, num_iterations);
path_points_fig = 3331;
fcn_Path_plotTraversalsXY(data,path_points_fig);
hold on
plot(path_average_final3.X,path_average_final3.Y,'Linewidth',4);
title('Final average path via orthogonal projections, N = 1')
xlabel('X [m]')
ylabel('Y [m]')

% Show the effect of 3 iterations
num_iterations = 3;
path_average_final3 = fcn_Path_findAverageTraversalViaOrthoProjection(data,reference_traversal, num_iterations);
path_points_fig = 3332;
fcn_Path_plotTraversalsXY(data,path_points_fig);
hold on
plot(path_average_final3.X,path_average_final3.Y,'Linewidth',4);
title('Final average path via orthogonal projections, N = 3')
xlabel('X [m]')
ylabel('Y [m]')

% Show the effect of 10 iterations
num_iterations = 10;
path_average_final3 = fcn_Path_findAverageTraversalViaOrthoProjection(data,reference_traversal, num_iterations);
path_points_fig = 3333;
fcn_Path_plotTraversalsXY(data,path_points_fig);
hold on
plot(path_average_final3.X,path_average_final3.Y,'Linewidth',4);
title('Final average path via orthogonal projections, N = 10')
xlabel('X [m]')
ylabel('Y [m]')

% Show the effect of 40 iterations
num_iterations = 40;
path_average_final3 = fcn_Path_findAverageTraversalViaOrthoProjection(data,reference_traversal, num_iterations);
path_points_fig = 3334;
fcn_Path_plotTraversalsXY(data,path_points_fig);
hold on
plot(path_average_final3.X,path_average_final3.Y,'Linewidth',4);
title('Final average path via orthogonal projections, N = 40')
xlabel('X [m]')
ylabel('Y [m]')

%% Example 3: Basic call showing effect of weight_for_averaging


clear data

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

reference_traversal = data.traversal{3};
num_iterations = 40;

% Show the effect of no averaging
path_points_fig = 3331;
weight_for_averaging = 0;
fcn_Path_findAverageTraversalViaOrthoProjection(data,reference_traversal, num_iterations,weight_for_averaging,path_points_fig);

% Show the effect of small averaging
path_points_fig = 3332;
weight_for_averaging = 0.5;
fcn_Path_findAverageTraversalViaOrthoProjection(data,reference_traversal, num_iterations,weight_for_averaging,path_points_fig);

% Show the effect of medium averaging
path_points_fig = 3333;
weight_for_averaging = 0.8;
fcn_Path_findAverageTraversalViaOrthoProjection(data,reference_traversal, num_iterations,weight_for_averaging,path_points_fig);

num_iterations = 80;
% Show the effect of heavy averaging
path_points_fig = 3334;
weight_for_averaging = 0.90;
fcn_Path_findAverageTraversalViaOrthoProjection(data,reference_traversal, num_iterations,weight_for_averaging,path_points_fig);

