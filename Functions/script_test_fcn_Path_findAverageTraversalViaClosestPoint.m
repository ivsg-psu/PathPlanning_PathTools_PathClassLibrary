% script_test_fcn_Path_findAverageTraversalViaClosestPoint
% Tests the function: fcn_Path_findAverageTraversalViaClosestPoint

% Revision history
%      2021_01_07
%      -- first write of the code

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

%% Example 1: Basic call
path_average_final_1 = fcn_Path_findAverageTraversalViaClosestPoint(data);

% Plot the final XY result of closest point
path_points_fig = 11111;
fcn_Path_plotTraversalsXY(data,path_points_fig);
hold on
plot(path_average_final_1.X,path_average_final_1.Y,'Linewidth',4);
title('Original paths and final average path via closest point averaging')
xlabel('X [m]')
ylabel('Y [m]')


%% Example 2: Set the reference traversal 
% and show they give same results if lots of iterations

clc;

reference_traversal_1 = data.traversal{1};
reference_traversal_2 = data.traversal{2};
reference_traversal_3 = data.traversal{3};

path_average_final_1 = fcn_Path_findAverageTraversalViaClosestPoint(...
    data,reference_traversal_1);
path_average_final_2 = fcn_Path_findAverageTraversalViaClosestPoint(...
    data,reference_traversal_2);
path_average_final_3 = fcn_Path_findAverageTraversalViaClosestPoint(...
    data,reference_traversal_3);

% Plot the final XY results
path_points_fig = 22222;
fcn_Path_plotTraversalsXY(data,path_points_fig);
hold on
plot(path_average_final_1.X,path_average_final_1.Y,'Linewidth',4);
plot(path_average_final_2.X,path_average_final_2.Y,'Linewidth',3);
plot(path_average_final_3.X,path_average_final_3.Y,'Linewidth',2);
title('Original paths and final averages for 3 different initial references')
xlabel('X [m]')
ylabel('Y [m]')

%% Example 3: Set the reference traversal and number of iterations
% and show they give different results if few iterations (Note: convergence
% usually is strong after 2 iterations)

num_iterations = 1;

close all;
clc;

reference_traversal_1 = data.traversal{1};
reference_traversal_2 = data.traversal{2};
reference_traversal_3 = data.traversal{3};

path_average_final_1 = fcn_Path_findAverageTraversalViaClosestPoint(...
    data,reference_traversal_1,num_iterations);
path_average_final_2 = fcn_Path_findAverageTraversalViaClosestPoint(...
    data,reference_traversal_2,num_iterations);
path_average_final_3 = fcn_Path_findAverageTraversalViaClosestPoint(...
    data,reference_traversal_3,num_iterations);

% Plot the final XY results
path_points_fig = 33333;
fcn_Path_plotTraversalsXY(data,path_points_fig);
hold on
plot(path_average_final_1.X,path_average_final_1.Y,'Linewidth',4);
plot(path_average_final_2.X,path_average_final_2.Y,'Linewidth',3);
plot(path_average_final_3.X,path_average_final_3.Y,'Linewidth',2);
title('Original paths and final averages for 3 different initial references')
xlabel('X [m]')
ylabel('Y [m]')

%% Example 4: Basic call with a figure number
% This turns on quite a bit of debug plotting
close all;
clc;

path_points_fig = 44444;
path_average_final_4 = fcn_Path_findAverageTraversalViaClosestPoint(...
    data,[],[],path_points_fig);

% Plot the final XY result of closest point
fcn_Path_plotTraversalsXY(data,path_points_fig);
hold on
plot(path_average_final_4.X,path_average_final_4.Y,'Linewidth',4);
title('Original paths and final average path via closest point averaging')
xlabel('X [m]')
ylabel('Y [m]')