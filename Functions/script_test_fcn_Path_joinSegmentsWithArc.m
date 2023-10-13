% script_test_fcn_Path_joinSegmentsWithArc
% Tests fcn_Path_joinSegmentsWithArc
       
% Revision history:
%  2023_10_09
%  -- first write of the code by S. Brennan (sbrennan@psu.edu)


%% clear the workspace
close all
clc

% Clear any old variables
clear all_traversals

%% Basic Example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   ____            _        ______                           _      
%  |  _ \          (_)      |  ____|                         | |     
%  | |_) | __ _ ___ _  ___  | |__  __  ____ _ _ __ ___  _ __ | | ___ 
%  |  _ < / _` / __| |/ __| |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \
%  | |_) | (_| \__ \ | (__  | |____ >  < (_| | | | | | | |_) | |  __/
%  |____/ \__,_|___/_|\___| |______/_/\_\__,_|_| |_| |_| .__/|_|\___|
%                                                      | |           
%                                                      |_|          
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Basic%20Example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§

%% Basic example 1
% simple line segment with straight lines and positive radius, 90 degree
% turn, adjoint to the 2nd segment

fig_num = 4501;
figure(fig_num);
clf;
grid on;
hold on;

segment_1 = [0 0; 1 1];
segment_2 = [2 3; 1 4];
number_of_points_on_curve = 50;

[arc_points, line_segment_1, line_segment_2,...
    radius_of_arc, circle_center,...
    start_angle_radians, end_angle_radians] = ...
    fcn_Path_joinSegmentsWithArc(...
    segment_1, ...
    segment_2,...
    number_of_points_on_curve,...
    fig_num);

% Check the size of the vectors
assert(length(arc_points(:,1))==number_of_points_on_curve);
assert(length(line_segment_1(1,:))==2);
assert(length(line_segment_2(1,:))==2);
assert(length(radius_of_arc(1,:))==1);
assert(length(circle_center(1,:))==2);
assert(length(start_angle_radians(1,:))==1);
assert(length(end_angle_radians(1,:))==1);

% Check the angles
assert((abs(end_angle_radians - start_angle_radians) - pi/2)<1E-8);

%% Basic example 2
% simple line segment with straight lines and positive radius, 90 degree
% turn, adjoint to the 2nd segment

fig_num = 4502;
segment_1 = [0 0; 2 2];
segment_2 = [1 4; 0 5];
number_of_points_on_curve = 5;

[arc_points, line_segment_1, line_segment_2,...
    radius_of_arc, circle_center,...
    start_angle_radians, end_angle_radians] = ...
    fcn_Path_joinSegmentsWithArc(...
    segment_1, ...
    segment_2,...
    number_of_points_on_curve,...
    fig_num);

% Check the size of the vectors
assert(length(arc_points(:,1))==number_of_points_on_curve);
assert(length(line_segment_1(1,:))==2);
assert(length(line_segment_2(1,:))==2);
assert(length(radius_of_arc(1,:))==1);
assert(length(circle_center(1,:))==2);
assert(length(start_angle_radians(1,:))==1);
assert(length(end_angle_radians(1,:))==1);

% Check the angles
assert((abs(end_angle_radians - start_angle_radians) - pi/2)<1E-8);

%% Basic example 3
% simple line segments that overlap

fig_num = 4503;
segment_1 = [0 0; 4 4];
segment_2 = [5 0; 0 5];
number_of_points_on_curve = 5;

[arc_points, line_segment_1, line_segment_2,...
    radius_of_arc, circle_center,...
    start_angle_radians, end_angle_radians] = ...
    fcn_Path_joinSegmentsWithArc(...
    segment_1, ...
    segment_2,...
    number_of_points_on_curve,...
    fig_num);

% Check the size of the vectors
assert(all(isnan(arc_points)));
assert(all(isnan(line_segment_1)));
assert(all(isnan(line_segment_2)));
assert(isequal(radius_of_arc,0)); 
assert(all(isnan(line_segment_2)));
assert(isequal(circle_center,[2.5 2.5])); 
assert(all(isnan(start_angle_radians)));
assert(all(isnan(end_angle_radians)));


%% Basic example 4
% simple line segment with straight lines and negative radius, 90 degree
% turn, adjoint to the 2nd segment

fig_num = 4504;
segment_1 = [0 0; 1 1];
segment_2 = [3 2; 4 1];
number_of_points_on_curve = 50;

[arc_points, line_segment_1, line_segment_2,...
    radius_of_arc, circle_center,...
    start_angle_radians, end_angle_radians] = ...
    fcn_Path_joinSegmentsWithArc(...
    segment_1, ...
    segment_2,...
    number_of_points_on_curve,...
    fig_num);

% Check the size of the vectors
assert(length(arc_points(:,1))==number_of_points_on_curve);
assert(length(line_segment_1(1,:))==2);
assert(length(line_segment_2(1,:))==2);
assert(length(radius_of_arc(1,:))==1);
assert(length(circle_center(1,:))==2);
assert(length(start_angle_radians(1,:))==1);
assert(length(end_angle_radians(1,:))==1);

% Check the angles
assert((abs(end_angle_radians - start_angle_radians) - pi/2)<1E-8);

%% Basic example 5
% simple line segment with straight lines and negative radius, 90 degree
% turn, adjoint to the 2nd segment

fig_num = 4505;
segment_1 = [0 0; 2 2];
segment_2 = [4 1; 5 0];
number_of_points_on_curve = 5;

[arc_points, line_segment_1, line_segment_2,...
    radius_of_arc, circle_center,...
    start_angle_radians, end_angle_radians] = ...
    fcn_Path_joinSegmentsWithArc(...
    segment_1, ...
    segment_2,...
    number_of_points_on_curve,...
    fig_num);

% Check the size of the vectors
assert(length(arc_points(:,1))==number_of_points_on_curve);
assert(length(line_segment_1(1,:))==2);
assert(length(line_segment_2(1,:))==2);
assert(length(radius_of_arc(1,:))==1);
assert(length(circle_center(1,:))==2);
assert(length(start_angle_radians(1,:))==1);
assert(length(end_angle_radians(1,:))==1);

% Check the angles
assert((abs(end_angle_radians - start_angle_radians) - pi/2)<1E-8);


%% Advanced case

fig_num = 4599;
figure(fig_num);
clf;
axis equal
hold on;
grid on;


segment_1 = [0 0; 2 2];
center = [2.5 2.5];
original_segment = [1 0; 3 0];


for ith_angle = -120:10:210
    angle_in_radians = ith_angle*pi/180;
    segment_2 = original_segment*[cos(angle_in_radians) sin(angle_in_radians); -sin(angle_in_radians) cos(angle_in_radians)] + ones(2,1)*center;
    number_of_points_on_curve = 5;

    [arc_points, line_segment_1, line_segment_2,...
        radius_of_arc, circle_center,...
        start_angle_radians, end_angle_radians] = ...
        fcn_Path_joinSegmentsWithArc(...
        segment_1, ...
        segment_2,...
        number_of_points_on_curve,...
        fig_num);

    axis([-1 5 -1 5]);
    pause(0.01);
end


%% Advanced case 2

fig_num = 45992;
figure(fig_num);
clf;
axis equal
hold on;
grid on;


segment_1 = [0 0; -2 0];
center = [-2.5 0];
original_segment = [-4 0; -7 0];


for ith_angle = -200:2:200
    angle_in_radians = ith_angle*pi/180;
    segment_2 = original_segment*[cos(angle_in_radians) sin(angle_in_radians); -sin(angle_in_radians) cos(angle_in_radians)] + ones(2,1)*center;
    number_of_points_on_curve = 5;

    [arc_points, line_segment_1, line_segment_2,...
        radius_of_arc, circle_center,...
        start_angle_radians, end_angle_radians] = ...
        fcn_Path_joinSegmentsWithArc(...
        segment_1, ...
        segment_2,...
        number_of_points_on_curve,...
        fig_num);

    axis([-10 5 -5 5]);
    pause(0.01);
end



%% FAIL CASES
if 1==0
    %% Case where paths completely overlap
    fig_num = 4500;
    segment_1 = [0 0; 4 4];
    segment_2 = segment_1;
    radius_of_curve = 1;
    plot_color = [];
    plot_line_width = [];
    plot_text = '';
    number_of_points_on_curve = 50;
    [ TransitionCurves, ...
        closest_path_point1, closest_path_point2,...
        angle_point1_radians,angle_point2_radians] = ...
        fcn_Path_joinSegmentsWithArc(...
        segment_1, ...
        segment_2,...
        radius_of_curve, ...
        number_of_points_on_curve,...
        plot_color,...
        plot_line_width,...
        plot_text,...
        fig_num);
end
