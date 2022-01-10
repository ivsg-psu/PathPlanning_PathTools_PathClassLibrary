% script_test_fcn_Path_fillOffsetTraversalsAboutTraversal
% Tests fcn_Path_fillOffsetTraversalsAboutTraversal
       
% Revision history:
%      2021_01_24
%      -- first write of the code


close all
clc

% Clear any old variables
clear all_traversals

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});
all_traversals.traversal{1} = reference_traversal;


% Plot the results? (Note: they are plotted below as well)
if 1==0
    fig_num = 12;
    fcn_Path_plotTraversalsYaw(all_traversals,fig_num);
    fig_num = 13;
    fcn_Path_plotTraversalsXY(all_traversals,fig_num);
end

%% Test case 1: basic call for one trajectory
offsets = 2; 
offset_traversals = fcn_Path_fillOffsetTraversalsAboutTraversal(reference_traversal, offsets);

figure(1);
plot(offset_traversals.traversal{1}.X,offset_traversals.traversal{1}.Y,'r.-','Linewidth',3);

%% Test case 2: basic call for one trajectory - specify figure
fig_num = 2;
offsets = 2; 
offset_traversals = fcn_Path_fillOffsetTraversalsAboutTraversal(reference_traversal, offsets,fig_num); %#ok<*NASGU>

%% Test case 3: basic call for two trajectories 
fig_num = 3;
offsets = [2; -2]; 
offset_traversals = fcn_Path_fillOffsetTraversalsAboutTraversal(reference_traversal, offsets,fig_num);

%% Test case 4: show how "pinching" can happen
fig_num = 4;
% Grab the "curve" of the path
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1}(13:20,:));
offsets = (-10:1:10)'; 
offset_traversals = fcn_Path_fillOffsetTraversalsAboutTraversal(reference_traversal, offsets,fig_num);
axis equal;

%% Test case 5: Show how offsets can link lane markers
fig_num = 5;
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1}(30:end,:));
path2 = [30 10; 25 44];
second_traversal = fcn_Path_convertPathToTraversalStructure(path2);

offsets = [2]; 
offset_traversal_1 = fcn_Path_fillOffsetTraversalsAboutTraversal(reference_traversal, offsets,fig_num);

offsets = [2; -2]; 
offset_traversal_2 = fcn_Path_fillOffsetTraversalsAboutTraversal(second_traversal, offsets,fig_num);

% Find intersections
[right_intersection_points,...
    ~,...
    ~] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    offset_traversal_1.traversal{1},...
    offset_traversal_2.traversal{1});
plot(right_intersection_points(:,1),right_intersection_points(:,2),'ro','Markersize',10);

[left_intersection_points,...
    ~,...
    ~] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    offset_traversal_1.traversal{1},...
    offset_traversal_2.traversal{2});
plot(left_intersection_points(:,1),left_intersection_points(:,2),'bo','Markersize',10);
title('Illustration of how to use offsets to link lane edge designations');







