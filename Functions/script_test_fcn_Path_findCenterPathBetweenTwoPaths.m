% script_test_fcn_Path_findCenterPathBetweenTwoPaths
% Tests the following:
%    center_path = ...
%     fcn_Path_findCenterPathBetweenTwoPaths(...
%     first_path,second_path,(flag_rounding_type),(search_radius),(fig_num))

% Revision history:
% 2023_09_04 by S. Brennan
% -- first write of the code

close all;
clc;

% Fill in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;


% Plot the results?
if 1==0
    fig_num = 12;
    fcn_Path_plotTraversalsYaw(data,fig_num);
    fig_num = 13;
    fcn_Path_plotTraversalsXY(data,fig_num);
end

%% Initialize the figure numbers
fig_num = 1000;


%% Demonstration of fcn_Path_findCenterPathBetweenTwoPaths
% This function finds the center projected from one traversal toward
% another
first_path = [0 0; 1 1; 2 1; 3 2];
second_path   = [0.5 1.5; 1.5 2.1; 4 6];

flag_rounding_type = 1;
search_radius = 10;
fig_num = fig_num+1;

figure(fig_num);
clf;
hold on;
grid on;
axis equal;

center_path = ...
    fcn_Path_findCenterPathBetweenTwoPaths(...
    first_path,second_path,(flag_rounding_type),(search_radius),(fig_num));

%% Demonstration of fcn_Path_findCenterPathBetweenTwoPaths
% This function finds the center projected from one traversal toward
% another
first_path = [0 0; 1 1; 2 1; 3 4];
second_path   = first_path + ones(length(first_path(:,1)),1)*[0 1];

flag_rounding_type = 1;
search_radius = 10;
fig_num = fig_num+1;

figure(fig_num);
clf;
hold on;
grid on;
axis equal;

center_path = ...
    fcn_Path_findCenterPathBetweenTwoPaths(...
    first_path,second_path,(flag_rounding_type),(search_radius),(fig_num));


%% Demonstration of effect of flag_rounding_type
% This function finds the center projected from one traversal toward
% another
first_path = [0 0; 1 1; 2 1; 3 2];
second_path   = [0.5 1.5; 1.5 2.1; 4 6];

fig_num = fig_num+1;
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

%% Call the center calculation function
first_path = paths{1};
second_path = paths{2};

flag_rounding_type = 1;
search_radius = 10;
fig_num = fig_num+1;


figure(fig_num);
clf;
hold on;
grid on;
axis equal;

center_path = ...
    fcn_Path_findCenterPathBetweenTwoPaths(...
    first_path,second_path,(flag_rounding_type),(search_radius),(fig_num));

