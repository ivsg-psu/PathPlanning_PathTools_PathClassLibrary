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

%% Basic example 1: simple line segment with straight lines and positive radius

fig_num = 4501;
segment_1 = [0 0; 1 1];
segment_2 = [2 3; 1 4];
plot_color = [];
plot_line_width = [];
plot_text = '';
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
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(length(closest_path_point1(:,1))==1);
assert(length(closest_path_point2(:,1))==1);

% Check the angles
assert(isequal(angle_point1_radians,-pi/4));
assert(abs(angle_point2_radians - pi/4)<eps*1000);

%% Basic example 1: simple line segment with straight lines and positive radius

fig_num = 4501;
segment_1 = [0 0; 4 4];
segment_2 = [4 1; 1 4];
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
    number_of_points_on_curve,...
    plot_color,...
    plot_line_width,...
    plot_text,...
    fig_num);

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(length(closest_path_point1(:,1))==1);
assert(length(closest_path_point2(:,1))==1);

% Check the angles
assert(isequal(angle_point1_radians,-pi/4));
assert(abs(angle_point2_radians - pi/4)<eps*1000);

%% Basic example 2: simple line segment with straight lines and negative radius

fig_num = 4502;
segment_1 = [0 0; 4 4];
segment_2 = flipud([4 1; 1 4]);  % <--- note direction change
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

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(length(closest_path_point1(:,1))==1);
assert(length(closest_path_point2(:,1))==1);
assert(length(angle_point1_radians(:,1))==1);
assert(length(angle_point2_radians(:,1))==1);

%% Basic example 3: simple line segment with line changing direction 

fig_num = 4503;
segment_1 = flipud([0 0; 4 4]); % <--- note direction change
segment_2 = flipud([4 1; 1 4]);  % <--- note direction change
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

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(length(closest_path_point1(:,1))==1);
assert(length(closest_path_point2(:,1))==1);
assert(length(angle_point1_radians(:,1))==1);
assert(length(angle_point2_radians(:,1))==1);

%% Basic example 4: simple line segment with line changing direction 

fig_num = 4504;
segment_1 = flipud([0 0; 4 4]); % <--- note direction change
segment_2 = [4 1; 1 4]; 
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

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(length(closest_path_point1(:,1))==1);
assert(length(closest_path_point2(:,1))==1);
assert(length(angle_point1_radians(:,1))==1);
assert(length(angle_point2_radians(:,1))==1);

%% Basic example 5a: simple curve, outside fit

fig_num = 45051;
th = linspace( pi/2, -pi/2, 100);
Radius = 10;  %or whatever radius you want
semi_circle_x = Radius*cos(th) + 5;
semi_circle_y = Radius*sin(th) + 4;
segment_1 = [semi_circle_x', semi_circle_y'];
segment_2 = [5 0; 30 0];
radius_of_curve = 6;
number_of_points_on_curve = 50;
plot_color = [];
plot_line_width = 3;
plot_text = 'Transition_Curve';

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

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(length(closest_path_point1(:,1))==1);
assert(length(closest_path_point2(:,1))==1);
assert(length(angle_point1_radians(:,1))==1);
assert(length(angle_point2_radians(:,1))==1);

%% Basic example 5b: simple curve, inside fit

fig_num = 45051;
figure(fig_num);
clf;
hold on;
grid on;


th = linspace( pi/2, -pi/2, 100);
Radius = 10;  %or whatever radius you want
semi_circle_x = Radius*cos(th) + 5;
semi_circle_y = Radius*sin(th) + 4;
segment_1 = [semi_circle_x', semi_circle_y'];
segment_2 = [5 0; 30 0];
radius_of_curve = 2;
number_of_points_on_curve = 50;
plot_color = [];
plot_line_width = 3;
plot_text = 'Transition_Curve';

[ TransitionCurves, ...
    closest_path_point1, closest_path_point2,...
    angle_point1_radians,angle_point2_radians] = ...
    fcn_Path_joinSegmentsWithArc(...
    segment_2, ...
    segment_1,...
    radius_of_curve, ...
    number_of_points_on_curve,...
    plot_color,...
    plot_line_width,...
    plot_text,...
    fig_num);

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(length(closest_path_point1(:,1))==1);
assert(length(closest_path_point2(:,1))==1);
assert(length(angle_point1_radians(:,1))==1);
assert(length(angle_point2_radians(:,1))==1);

%% Basic example 5c: simple curve, inside fit, other direction

fig_num = 45051;
figure(fig_num);
clf;
hold on;
grid on;


th = linspace( pi/2, -pi/2, 100);
Radius = 10;  %or whatever radius you want
semi_circle_x = Radius*cos(th) + 5;
semi_circle_y = Radius*sin(th) + 4;
segment_1 = [semi_circle_x', semi_circle_y'];
segment_2 = [5 0; 30 0];
radius_of_curve = 2;
number_of_points_on_curve = 50;
plot_color = [];
plot_line_width = 3;
plot_text = 'Transition_Curve';

[ TransitionCurves, ...
    closest_path_point1, closest_path_point2,...
    angle_point1_radians,angle_point2_radians] = ...
    fcn_Path_joinSegmentsWithArc(...
    segment_2, ...
    flipud(segment_1),...
    radius_of_curve, ...
    number_of_points_on_curve,...
    plot_color,...
    plot_line_width,...
    plot_text,...
    fig_num);

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(length(closest_path_point1(:,1))==1);
assert(length(closest_path_point2(:,1))==1);
assert(length(angle_point1_radians(:,1))==1);
assert(length(angle_point2_radians(:,1))==1);


%% Basic example 5d: simple curve, inside fit, other direction, coarse

fig_num = 45051;
figure(fig_num);
clf;
hold on;
grid on;


th = linspace( pi/2, -pi/2, 5);
Radius = 10;  %or whatever radius you want
semi_circle_x = Radius*cos(th) + 5;
semi_circle_y = Radius*sin(th) + 4;
segment_1 = [semi_circle_x', semi_circle_y'];
segment_2 = [5 0; 30 0];
radius_of_curve = 5;
number_of_points_on_curve = 50;
plot_color = [];
plot_line_width = 3;
plot_text = 'Transition_Curve';

[ TransitionCurves, ...
    closest_path_point1, closest_path_point2,...
    angle_point1_radians,angle_point2_radians] = ...
    fcn_Path_joinSegmentsWithArc(...
    segment_2, ...
    flipud(segment_1),...
    radius_of_curve, ...
    number_of_points_on_curve,...
    plot_color,...
    plot_line_width,...
    plot_text,...
    fig_num);

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(length(closest_path_point1(:,1))==1);
assert(length(closest_path_point2(:,1))==1);
assert(length(angle_point1_radians(:,1))==1);
assert(length(angle_point2_radians(:,1))==1);


%%  Basic example 6a: multiple intersections between the line segments 

fig_num = 45061;
segment_1 = [2 2; 5 8];
segment_2 = [3 2; 3 6; 5 6];
radius_of_curve = 2;

plot_color = [];
plot_line_width = [];
plot_text = [];

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

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(length(closest_path_point1(:,1))==1);
assert(length(closest_path_point2(:,1))==1);
assert(length(angle_point1_radians(:,1))==1);
assert(length(angle_point2_radians(:,1))==1);

%%  Basic example 6b: multiple intersections between the line segments 

fig_num = 45062;
segment_1 = [2 2; 5 8];
segment_2 = [5 2; 2 6; 5 6];
radius_of_curve = 2;

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

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(length(closest_path_point1(:,1))==1);
assert(length(closest_path_point2(:,1))==1);
assert(length(angle_point1_radians(:,1))==1);
assert(length(angle_point2_radians(:,1))==1);

%%  Basic example 7: multiple intersections between the line segments 

fig_num = 4507;
segment_1 = [2 2; 5 8];
segment_2 = [5 2; 2 6; 5 6];
radius_of_curve = 2;
number_of_points_on_curve = 50;
plot_color = [];
plot_line_width = [];
plot_text = '';

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

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(length(closest_path_point1(:,1))==1);
assert(length(closest_path_point2(:,1))==1);
assert(length(angle_point1_radians(:,1))==1);
assert(length(angle_point2_radians(:,1))==1);

%% HARD example 8: line segments have infinite intersections for a period

fig_num = 4508;
segment_1 = [0 0; 4 4; 6 6; 10 5];
segment_2 = [4 1; 4 3; 5 5; 7 7; 5 10];
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

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(length(closest_path_point1(:,1))==1);
assert(length(closest_path_point2(:,1))==1);
assert(length(angle_point1_radians(:,1))==1);
assert(length(angle_point2_radians(:,1))==1);


%% Example for which code breaks
% should get error: Lines are collinear, curve is not possible between them

fig_num = 4506;
segment_1 = [4 9; 4 5; 2 2];
segment_2 = [4 9; 4 2];
radius_of_curve = 2;
number_of_points_on_curve = 50;
plot_color = [];
plot_line_width = [];
plot_text = '';

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

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(length(closest_path_point1(:,1))==1);
assert(length(closest_path_point2(:,1))==1);
assert(length(angle_point1_radians(:,1))==1);
assert(length(angle_point2_radians(:,1))==1);


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
%% function to remove NaNs from ENU data
function data_no_nans = fcn_INTERNAL_remove_NaNs(data_with_NaNs)

% INPUT
% data_with_NaNs : ENU data that has NaNs

% OUTPUT
% data_no_nans : ENU data without NaNs

data_no_nans = data_with_NaNs(~isnan(data_with_NaNs(:,1)),:); % taking all data that is not NaNs and putting it in the variable called data_no_nans
data_no_nans = [data_no_nans(:,1), data_no_nans(:,2)]; % take only first 2 columns 

end