% script_test_fcn_Path_findCenterPathBetweenTwoPaths
% Tests the following:
%    center_path = ...
%     fcn_Path_findCenterPathBetweenTwoPaths(...
%     first_path,second_path,(flag_rounding_type),(search_radius),(fig_num))

% Revision history:
% 2023_09_04 by S. Brennan
% -- first write of the code

close all;



%% Basic demo 1 - Demonstration of fcn_Path_findCenterPathBetweenTwoPaths
% This function finds the center projected from one traversal toward
% another
first_path = [0 0; 1 1; 2 1; 3 2];
second_path   = [0.5 1.5; 1.5 2.1; 4 6];

flag_rounding_type = 1;
search_radius = 10;
core_fig_num = 1000;
fig_num = core_fig_num+1;

figure(fig_num);
clf;
hold on;
grid on;
axis equal;

center_path = ...
    fcn_Path_findCenterPathBetweenTwoPaths(...
    first_path,second_path,(flag_rounding_type),(search_radius),(fig_num));
assert(length(center_path(:,1)) == (length(first_path(:,1)) + length(second_path(:,1))));

%% Demonstration of fcn_Path_findCenterPathBetweenTwoPaths
% This function finds the center projected from one traversal toward
% another
first_path = [0 0; 1 1; 2 1; 3 4];
second_path   = first_path + ones(length(first_path(:,1)),1)*[0 1];

flag_rounding_type = 1;
search_radius = 10;
core_fig_num = 1000;
fig_num = core_fig_num+2;

figure(fig_num);
clf;
hold on;
grid on;
axis equal;

center_path = ...
    fcn_Path_findCenterPathBetweenTwoPaths(...
    first_path,second_path,(flag_rounding_type),(search_radius),(fig_num));
assert(length(center_path(:,1)) == (length(first_path(:,1)) + length(second_path(:,1))));

%% Demonstration of effect of flag_rounding_type
% This function finds the center projected from one traversal toward
% another
first_path = [0 0; 1 1; 2 1; 3 2];
second_path   = [0.5 1.5; 1.5 2.1; 4 6];

core_fig_num = 1000;
fig_num = core_fig_num+3;
search_radius = 10;


figure(fig_num);
clf;


for flag_rounding_type = 1:3
    subplot(2,2,flag_rounding_type);
    hold on;
    grid on;
    axis equal;
    center_path = ...
        fcn_Path_findCenterPathBetweenTwoPaths(...
        first_path,second_path,(flag_rounding_type),(search_radius),(fig_num));
    title(sprintf('flag_rounding_type: %.0d',flag_rounding_type),'Interpreter','none');
end
sgtitle('Effect of flag_rounding_type','Interpreter','none');
assert(length(center_path(:,1)) == (length(first_path(:,1)) + length(second_path(:,1))));

%% Test case where first path is not intersecting second path with any projections

first_path = [0 0; 10 0];
second_path   = [1 1; 9 1];

core_fig_num = 1000;
fig_num = core_fig_num+4;
flag_rounding_type = 1;
search_radius = 10;


figure(fig_num);
clf;
hold on;
grid on;
axis equal;

center_path = ...
    fcn_Path_findCenterPathBetweenTwoPaths(...
    first_path,second_path,(flag_rounding_type),(search_radius),(fig_num));
assert(length(center_path(:,1)) == (length(first_path(:,1)) + length(second_path(:,1))));

%% Test case where second path is not intersecting first path with any projections

second_path = [0 0; 10 0];
first_path   = [1 1; 9 1];

core_fig_num = 1000;
fig_num = core_fig_num+5;
flag_rounding_type = 1;
search_radius = 10;


figure(fig_num);
clf;
hold on;
grid on;
axis equal;

center_path = ...
    fcn_Path_findCenterPathBetweenTwoPaths(...
    first_path,second_path,(flag_rounding_type),(search_radius),(fig_num));
assert(length(center_path(:,1)) == (length(first_path(:,1)) + length(second_path(:,1))));

%% Test case where first path and second path splay
first_path = [-5 0; 0 0; 10 0; 30 0];
second_path   = [-5 5; 1 5; 11 11; 30 11];

core_fig_num = 1000;
fig_num = core_fig_num+6;
flag_rounding_type = 1;
search_radius = 20;


figure(fig_num);
clf;
hold on;
grid on;
axis equal;

center_path = ...
    fcn_Path_findCenterPathBetweenTwoPaths(...
    first_path,second_path,(flag_rounding_type),(search_radius),(fig_num));
assert(length(center_path(:,1)) == (length(first_path(:,1)) + length(second_path(:,1))));
%% Call the center calculation function with a real-world path
paths = fcn_Path_fillSamplePaths;
first_path = paths{1};
second_path = paths{2};

flag_rounding_type = 1;
search_radius = 10;
core_fig_num = 1000;
fig_num = core_fig_num+7;


figure(fig_num);
clf;
hold on;
grid on;
axis equal;

center_path = ...
    fcn_Path_findCenterPathBetweenTwoPaths(...
    first_path,second_path,(flag_rounding_type),(search_radius),(fig_num));
assert(length(center_path(:,1)) == (length(first_path(:,1)) + length(second_path(:,1))));



%% FAIL CASES
if 1==0
    %% Show failure when neither path has projections that intersect each other

    first_path = [0 0; 10 0];
    second_path   = [-1 0; -10 10];

    core_fig_num = 1000;
    fig_num = core_fig_num+8;
    flag_rounding_type = 1;
    search_radius = 10;


    figure(fig_num);
    clf;
    hold on;
    grid on;
    axis equal;

    center_path = ...
        fcn_Path_findCenterPathBetweenTwoPaths(...
        first_path,second_path,(flag_rounding_type),(search_radius),(fig_num));
    
end