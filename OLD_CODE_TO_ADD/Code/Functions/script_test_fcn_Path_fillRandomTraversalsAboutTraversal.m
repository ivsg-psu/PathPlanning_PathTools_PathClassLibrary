% script_test_fcn_Path_fillRandomTraversalsAboutTraversal.m
% Tests fcn_Path_fillRandomTraversalsAboutTraversal
       
% Revision history:
%      2021_01_03
%      -- first write of the code
%      2021_01_07
%      -- cleared up function calls for traversals vs paths

close all
clc

% Clear any old variables
clear all_traversals

% Fill in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths{1});
all_traversals.traversal{1} = reference_traversal;


% Plot the results? (Note: they are plotted below as well)
if 1==0
    fig_num = 12;
    fcn_Path_plotTraversalsYaw(all_traversals,fig_num);
    fig_num = 13;
    fcn_Path_plotTraversalsXY(all_traversals,fig_num);
end

%% Test case 1: basic call for one trajectory
random_traversals = fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal);

figure(1);
plot(random_traversals.traversal{1}.X,random_traversals.traversal{1}.Y,'r.-','Linewidth',3);

%% Test case 2: advanced call for one trajectory - specify figure
%      random_traversals = ...
%      fcn_Path_fillRandomTraversalsAboutTraversal(...
%            reference_traversal,...
%            (std_deviation),...
%            (num_trajectories),...
%            (num_points),...
%            (flag_generate_random_stations),...
%            (spatial_smoothness),...
%            (fig_num));
emtpy_value = [];
fig_num = 2;
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    emtpy_value,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    emtpy_value,... % (flag_generate_random_stations),...
    emtpy_value,... % (spatial_smoothness),...
    fig_num);

%% Test case 3: use same station points
%      random_traversals = ...
%      fcn_Path_fillRandomTraversalsAboutTraversal(...
%            reference_traversal,...
%            (std_deviation),...
%            (num_trajectories),...
%            (num_points),...
%            (flag_generate_random_stations),...
%            (spatial_smoothness),...
%            (fig_num));
emtpy_value = [];
fig_num = 3;
flag_generate_random_stations = 0;
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    emtpy_value,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    emtpy_value,... % (spatial_smoothness),...
    fig_num);

%% Test case 4: show effects of spatial smoothness with many trajectories
%      random_traversals = ...
%      fcn_Path_fillRandomTraversalsAboutTraversal(...
%            reference_traversal,...
%            (std_deviation),...
%            (num_trajectories),...
%            (num_points),...
%            (flag_generate_random_stations),...
%            (spatial_smoothness),...
%            (fig_num));
emtpy_value = [];
flag_generate_random_stations = 0;
num_trajectories = 5;

fig_num = 41;
spatial_smoothness = 2;  % Units are meters
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title('Spatial smoothness: 2 meters (generates warning)');

fig_num = 42;
spatial_smoothness = 5;  % Units are meters
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('Spatial smoothness: %.0d meters',spatial_smoothness));


fig_num = 43;
spatial_smoothness = 10;  % Units are meters
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('Spatial smoothness: %.0d meters',spatial_smoothness));

fig_num = 44;
spatial_smoothness = 15;  % Units are meters
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('Spatial smoothness: %.0d meters',spatial_smoothness));

fig_num = 45;
spatial_smoothness = 25;  % Units are meters
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('Spatial smoothness: %.0d meters',spatial_smoothness));


%% Test case 5: show effects of standard deviation
%      random_traversals = ...
%      fcn_Path_fillRandomTraversalsAboutTraversal(...
%            reference_traversal,...
%            (std_deviation),...
%            (num_trajectories),...
%            (num_points),...
%            (flag_generate_random_stations),...
%            (spatial_smoothness),...
%            (fig_num));
emtpy_value = [];
flag_generate_random_stations = 0;
num_trajectories = 5;
spatial_smoothness = 7;  % Units are meters

fig_num = 51;
std_deviation = 1;  % Units are meters
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    std_deviation,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('Standard deviation: %.0d meters',std_deviation));

fig_num = 52;
std_deviation = 2;  % Units are meters
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    std_deviation,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('Standard deviation: %.0d meters',std_deviation));

fig_num = 53;
std_deviation = 5;  % Units are meters
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    std_deviation,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('Standard deviation: %.0d meters',std_deviation));

%% Test case 6: show effects of num_points
%      random_traversals = ...
%      fcn_Path_fillRandomTraversalsAboutTraversal(...
%            reference_traversal,...
%            (std_deviation),...
%            (num_trajectories),...
%            (num_points),...
%            (flag_generate_random_stations),...
%            (spatial_smoothness),...
%            (fig_num));
emtpy_value = [];
flag_generate_random_stations = 1;
num_trajectories = 5;
spatial_smoothness = 10;  % Units are meters

fig_num = 61;
num_points = 10;
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    num_points,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('num points: %.0d',num_points));

fig_num = 62;
num_points = 50;
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    num_points,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('num points: %.0d',num_points));

fig_num = 63;
num_points = 200;
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    num_points,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('num points: %.0d',num_points));