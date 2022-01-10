% script_test_fcn_Path_removePinchPointInTraversal.m
% Tests fcn_Path_removePinchPointInTraversal
       
% Revision history:
%      2021_01_24
%      -- first write of the code


close all
clc


%% Test case 1: simple test case
fig_num = 1;
path = [0 0; 0 4; 2 4; -1 1];
traversal_with_pinch_point = fcn_Path_convertPathToTraversalStructure(path);
[traversal_no_pinch_point] = ...
    fcn_Path_removePinchPointInTraversal(...
    traversal_with_pinch_point,...
    fig_num); %#ok<*NASGU>


%% Test case 2: simple test case producing warning
fig_num = 2;
path = [0 0; 0 4; 4 4; -1 1];
traversal_with_pinch_point = fcn_Path_convertPathToTraversalStructure(path);
[traversal_no_pinch_point] = ...
    fcn_Path_removePinchPointInTraversal(...
    traversal_with_pinch_point,...
    fig_num);

%% Test case 3: simple test case - no intersections
fig_num = 3;
path = [0 0; 0 4; 4 5; 5 5];
traversal_with_pinch_point = fcn_Path_convertPathToTraversalStructure(path);
[traversal_no_pinch_point] = ...
    fcn_Path_removePinchPointInTraversal(...
    traversal_with_pinch_point,...
    fig_num);

%% Test case 4: simple test case with multiple crossings
fig_num = 4;
path = [0 0; 0 4; 2 4; -1 1; -1 3; -0.5 1];
traversal_with_pinch_point = fcn_Path_convertPathToTraversalStructure(path);
[traversal_no_pinch_point] = ...
    fcn_Path_removePinchPointInTraversal(...
    traversal_with_pinch_point,...
    fig_num);
axis([-2 4 -2 6]);



%% Test case 5: simple test case with multiple crossings
fig_num = 5;
path = [0 0; 0 4; 2 4; -1 1; -1 3; -0.5 1; -0.5 2.5];
traversal_with_pinch_point = fcn_Path_convertPathToTraversalStructure(path);
[traversal_no_pinch_point] = ...
    fcn_Path_removePinchPointInTraversal(...
    traversal_with_pinch_point,...
    fig_num);

%% Test case 6: crossing back toward itself
fig_num = 6;
path = [0 0; 0 4; 0 2];
traversal_with_pinch_point = fcn_Path_convertPathToTraversalStructure(path);
[traversal_no_pinch_point] = ...
    fcn_Path_removePinchPointInTraversal(...
    traversal_with_pinch_point,...
    fig_num);

%% Test case 7: crossing back toward itself and continuing
fig_num = 7;
path = [0 0; 0 4; 0 2; 0 5];
traversal_with_pinch_point = fcn_Path_convertPathToTraversalStructure(path);
[traversal_no_pinch_point] = ...
    fcn_Path_removePinchPointInTraversal(...
    traversal_with_pinch_point,...
    fig_num);

%% Test case 8: crossing back toward itself and continuing elsewhere
fig_num = 8;
path = [0 0; 0 4; 0 2; 1 1];
traversal_with_pinch_point = fcn_Path_convertPathToTraversalStructure(path);
[traversal_no_pinch_point] = ...
    fcn_Path_removePinchPointInTraversal(...
    traversal_with_pinch_point,...
    fig_num);

%% Test case 9: crossing back toward itself without being in sequence
fig_num = 9;
path = [0 0; 0 4; 4 4; 0 0; -1 0];
traversal_with_pinch_point = fcn_Path_convertPathToTraversalStructure(path);
[traversal_no_pinch_point] = ...
    fcn_Path_removePinchPointInTraversal(...
    traversal_with_pinch_point,...
    fig_num);

%% Test case 10: crossing back toward itself without being in sequence
fig_num = 10;
path = [0 0; 0 4; 4 4; 0 0; 0 2];
traversal_with_pinch_point = fcn_Path_convertPathToTraversalStructure(path);
[traversal_no_pinch_point] = ...
    fcn_Path_removePinchPointInTraversal(...
    traversal_with_pinch_point,...
    fig_num);


%% Advanced test case 3: show a pinch point in practice
fig_num = 1111;
figure(fig_num);
clf;
axis equal

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

% Grab the "curve" of the path
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1}(13:20,:));
offsets = (0:1:10)'; 
offset_traversals = fcn_Path_fillOffsetTraversalsAboutTraversal(reference_traversal, offsets,fig_num);

% Fill in an array of "fixed" traversals
clear fixed_traversals
for ith_traversal = 1:length(offset_traversals.traversal)
    traversal_with_pinch_point = offset_traversals.traversal{ith_traversal};
    [traversal_no_pinch_point] = ...
        fcn_Path_removePinchPointInTraversal(...
        traversal_with_pinch_point);
    fixed_traversals.traversal{ith_traversal} = traversal_no_pinch_point; 
end

% Plot the results
fixed_fig_num = 2222;
figure(fixed_fig_num);
clf;
axis equal
hold on;
plot(reference_traversal.X,reference_traversal.Y,'b','Linewidth',3);
fcn_Path_plotTraversalsXY(fixed_traversals,fixed_fig_num)

