% script_test_fcn_Path_snapPointToPathViaVectors.m
% This is a script to exercise the function: fcn_Path_snapPointToPathViaVectors.m
% This function was written on 2020_10_10 by S. Brennan, sbrennan@psu.edu

% Revision history:
% 2023_09_28 - S. Brennan
% -- documented bug via test case wherein closest path not correctly found
% in prior snap function, fcn_Path_snapPointOntoNearestPath, and copied the
% test script from fcn_Path_snapPointOntoNearestPath into here.


close all;

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


%% BASIC example
% A simple line segment, a simple query, zero distance in rear segments
fig_num = 10001;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

points = [0.5 0.5];
pathXY = [0 0; 2 2];
flag_snap_type = 1;

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length,distance_real,distance_imaginary] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);


% Make sure function worked
true_closest_path_point = [0.5 0.5];
true_s_coordinate = 0.5*(2^0.5);
true_first_path_point_index = 1;
true_second_path_point_index = 2;
true_percent_along_length = 0.25;
true_distance_real = 0;
true_distance_imaginary = 0;

assert(isequal(round(closest_path_point,4),round(true_closest_path_point,4)));
assert(isequal(round(s_coordinate,4),round(true_s_coordinate,4)));
assert(isequal(round(first_path_point_index,4),round(true_first_path_point_index,4)));
assert(isequal(round(second_path_point_index,4),round(true_second_path_point_index,4)));
assert(isequal(round(percent_along_length,4),round(true_percent_along_length,4)));
assert(isequal(round(distance_real,4),round(true_distance_real,4)));
assert(isequal(round(distance_imaginary,4),round(true_distance_imaginary,4)));


fprintf(1,['Figure: %d,\n\t\t Closest point is: %.2f %.2f \n' ...
    '\t\t Matched to the path segment given by indices %d and %d, \n' ...
    '\t\t S-coordinate is: %.2f, \n' ...
    '\t\t percent_along_length is: %.2f\n' ...
    '\t\t real distance is: %.2f\n, ' ...
    '\t\t imag distance is %.2f\n, '],...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length,...
    distance_real,distance_imaginary);


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% HARD one - was not working on 2023-09-29, now fixed
fig_num = 10002;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

flag_snap_type = 1;
points = [-2 -1]; 
pathXY = [-1 0; 1 0];


% Snap the point onto the path
[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length,distance_real,distance_imaginary] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);


% Make sure function worked
true_closest_path_point = [-1 0];
true_s_coordinate = -1;
true_first_path_point_index = 1;
true_second_path_point_index = 1;
true_percent_along_length = -0.50;
true_distance_real = -1;
true_distance_imaginary = 1;

assert(isequal(round(closest_path_point,4),round(true_closest_path_point,4)));
assert(isequal(round(s_coordinate,4),round(true_s_coordinate,4)));
assert(isequal(round(first_path_point_index,4),round(true_first_path_point_index,4)));
assert(isequal(round(second_path_point_index,4),round(true_second_path_point_index,4)));
assert(isequal(round(percent_along_length,4),round(true_percent_along_length,4)));
assert(isequal(round(distance_real,4),round(true_distance_real,4)));
assert(isequal(round(distance_imaginary,4),round(true_distance_imaginary,4)));

% Print results to the workspace
fprintf(1,'Figure: %d\n',fig_num);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Testing percentage after end
fig_num = 10003;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

flag_snap_type = 1;
points = [2 1]; 
pathXY = [-1 0; 1 0];

% Snap the point onto the path
[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length,distance_real,distance_imaginary] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);


% Make sure function worked
true_closest_path_point = [1 0];
true_s_coordinate = 3;
true_first_path_point_index = 2;
true_second_path_point_index = 2;
true_percent_along_length = 1.50;
true_distance_real = 1;
true_distance_imaginary = -1;

assert(isequal(round(closest_path_point,4),round(true_closest_path_point,4)));
assert(isequal(round(s_coordinate,4),round(true_s_coordinate,4)));
assert(isequal(round(first_path_point_index,4),round(true_first_path_point_index,4)));
assert(isequal(round(second_path_point_index,4),round(true_second_path_point_index,4)));
assert(isequal(round(percent_along_length,4),round(true_percent_along_length,4)));
assert(isequal(round(distance_real,4),round(true_distance_real,4)));
assert(isequal(round(distance_imaginary,4),round(true_distance_imaginary,4)));

% Print results to the workspace
fprintf(1,'Figure: %d\n',fig_num);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);

%% HARD one 
% A 90-degree line segment with multiple surrounding queries
%XY_points = XY_points(6,:);
fig_num = 10004;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

flag_snap_type = 1;
points = [2 0]; 
pathXY = [-1 0; 1 0; 1 -1];

% Snap the point onto the path
[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length,distance_real,distance_imaginary] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);


% Make sure function worked
true_closest_path_point = [1 0];
true_s_coordinate = 2;
true_first_path_point_index = 2;
true_second_path_point_index = 2;
true_percent_along_length = 0;
true_distance_real = 0;
true_distance_imaginary = -1;

assert(isequal(round(closest_path_point,4),round(true_closest_path_point,4)));
assert(isequal(round(s_coordinate,4),round(true_s_coordinate,4)));
assert(isequal(round(first_path_point_index,4),round(true_first_path_point_index,4)));
assert(isequal(round(second_path_point_index,4),round(true_second_path_point_index,4)));
assert(isequal(round(percent_along_length,4),round(true_percent_along_length,4)));
assert(isequal(round(distance_real,4),round(true_distance_real,4)));
assert(isequal(round(distance_imaginary,4),round(true_distance_imaginary,4)));

% Print results to the workspace
fprintf(1,'Figure: %d\n',fig_num);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% HARD one - was not working on 2023-08-27, now fixed
fig_num = 10005;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

flag_snap_type = 1;
points = [2 0]; 
pathXY = [-1 0; 1 0; 1 -1];

% Snap the point onto the path
[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length,distance_real,distance_imaginary] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);


% Make sure function worked
true_closest_path_point = [1 0];
true_s_coordinate = 2;
true_first_path_point_index = 2;
true_second_path_point_index = 2;
true_percent_along_length = 0;
true_distance_real = 0;
true_distance_imaginary = -1;

assert(isequal(round(closest_path_point,4),round(true_closest_path_point,4)));
assert(isequal(round(s_coordinate,4),round(true_s_coordinate,4)));
assert(isequal(round(first_path_point_index,4),round(true_first_path_point_index,4)));
assert(isequal(round(second_path_point_index,4),round(true_second_path_point_index,4)));
assert(isequal(round(percent_along_length,4),round(true_percent_along_length,4)));
assert(isequal(round(distance_real,4),round(true_distance_real,4)));
assert(isequal(round(distance_imaginary,4),round(true_distance_imaginary,4)));

% Print results to the workspace
fprintf(1,'Figure: %d\n',fig_num);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% HARD one - was not working on 2023-08-27, now fixed
fig_num = 10006;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

flag_snap_type = 2;
points = [1 1]; 
pathXY = [-1 0; 1 0; 1 -1];


% Snap the point onto the path
[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length,distance_real,distance_imaginary] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);


% Make sure function worked
true_closest_path_point = [1 0];
true_s_coordinate = 2;
true_first_path_point_index = 2;
true_second_path_point_index = 2;
true_percent_along_length = 0;
true_distance_real = 0;
true_distance_imaginary = 1;

assert(isequal(round(closest_path_point,4),round(true_closest_path_point,4)));
assert(isequal(round(s_coordinate,4),round(true_s_coordinate,4)));
assert(isequal(round(first_path_point_index,4),round(true_first_path_point_index,4)));
assert(isequal(round(second_path_point_index,4),round(true_second_path_point_index,4)));
assert(isequal(round(percent_along_length,4),round(true_percent_along_length,4)));
assert(isequal(round(distance_real,4),round(true_distance_real,4)));
assert(isequal(round(distance_imaginary,4),round(true_distance_imaginary,4)));

% Print results to the workspace
fprintf(1,'Figure: %d\n',fig_num);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Known fail case - requires different function approach
% There is a point ON the path, and yet it is snapping to a different
% portion of the path. This is because the snap function FIRST finds the
% closest vertx point on the path to the query point, and if this closet
% point is on a path segment nearby but NOT the closest, the closest path
% will not be checked.
fig_num = 10007;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

flag_snap_type = 1;
points = [1 1]; 
pathXY = [-5 1; 5 1; 5 0; 2 0; 2 2];

% Snap the point onto the path
[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length,distance_real,distance_imaginary] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);


% Make sure function worked
true_closest_path_point = [1 1];
true_s_coordinate = 6;
true_first_path_point_index = 1;
true_second_path_point_index = 2;
true_percent_along_length = 0.6;
true_distance_real = 0;
true_distance_imaginary = 0;

assert(isequal(round(closest_path_point,4),round(true_closest_path_point,4)));
assert(isequal(round(s_coordinate,4),round(true_s_coordinate,4)));
assert(isequal(round(first_path_point_index,4),round(true_first_path_point_index,4)));
assert(isequal(round(second_path_point_index,4),round(true_second_path_point_index,4)));
assert(isequal(round(percent_along_length,4),round(true_percent_along_length,4)));
assert(isequal(round(distance_real,4),round(true_distance_real,4)));
assert(isequal(round(distance_imaginary,4),round(true_distance_imaginary,4)));

% Print results to the workspace
fprintf(1,'Figure: %d\n',fig_num);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example 0.1
% A simple line segment, a simple query, zero distance in rear segments,
% negative percent length
fig_num = 10008;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

points = [-0.5 -0.5];
pathXY = [0 0;2 2];
flag_snap_type = 1;

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length,distance_real,distance_imaginary] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);
fprintf(1,['Figure: %d,\n\t\t Closest point is: %.2f %.2f \n' ...
    '\t\t Matched to the path segment given by indices %d and %d, \n' ...
    '\t\t S-coordinate is: %.2f, \n' ...
    '\t\t percent_along_length is: %.2f\n' ...
    '\t\t real distance is: %.2f\n, ' ...
    '\t\t imag distance is %.2f\n, '],...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length,...
    distance_real,distance_imaginary);


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example 0.2
% A simple line segment, a simple query, zero distance in rear segments,
% positive percent length over 100%
fig_num = 10009;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

points = [3 3];
pathXY = [0 0;2 2];
flag_snap_type = 1;

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length,distance_real,distance_imaginary] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);
fprintf(1,['Figure: %d,\n\t\t Closest point is: %.2f %.2f \n' ...
    '\t\t Matched to the path segment given by indices %d and %d, \n' ...
    '\t\t S-coordinate is: %.2f, \n' ...
    '\t\t percent_along_length is: %.2f\n' ...
    '\t\t real distance is: %.2f\n, ' ...
    '\t\t imag distance is %.2f\n, '],...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length,...
    distance_real,distance_imaginary);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example 1
% A simple line segment, a simple query, positive distance in rear segments
fig_num = 10010;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

points = [0.5 1.5];
pathXY = [0 0;3 0; 5 2; 8 2];
flag_snap_type = 1;

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length,distance_real,distance_imaginary] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);
fprintf(1,['Figure: %d,\n\t\t Closest point is: %.2f %.2f \n' ...
    '\t\t Matched to the path segment given by indices %d and %d, \n' ...
    '\t\t S-coordinate is: %.2f, \n' ...
    '\t\t percent_along_length is: %.2f\n' ...
    '\t\t real distance is: %.2f\n, ' ...
    '\t\t imag distance is %.2f\n, '],...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length,...
    distance_real,distance_imaginary);


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example 1.1
% A simple line segment, a simple query, negative distance in both front
% and rear segments
fig_num = 10011;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;


points = [4.5 1];
pathXY = [0 0;3 0; 5 2; 8 2];

flag_snap_type = 1;

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length,distance_real,distance_imaginary] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);
fprintf(1,['Figure: %d,\n\t\t Closest point is: %.2f %.2f \n' ...
    '\t\t Matched to the path segment given by indices %d and %d, \n' ...
    '\t\t S-coordinate is: %.2f, \n' ...
    '\t\t percent_along_length is: %.2f\n' ...
    '\t\t real distance is: %.2f\n, ' ...
    '\t\t imag distance is %.2f\n, '],...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length,...
    distance_real,distance_imaginary);


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example 1.2
% A simple line segment, a simple query, negative distance in both front
% and rear segments
points = [5.4 1.5];
pathXY = [0 0;3 0; 5 2; 8 2];

flag_snap_type = 1;

fig_num = 111;
[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length,distance_real,distance_imaginary] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);
fprintf(1,['Figure: %d,\n\t\t Closest point is: %.2f %.2f \n' ...
    '\t\t Matched to the path segment given by indices %d and %d, \n' ...
    '\t\t S-coordinate is: %.2f, \n' ...
    '\t\t percent_along_length is: %.2f\n' ...
    '\t\t real distance is: %.2f\n, ' ...
    '\t\t imag distance is %.2f\n, '],...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length,...
    distance_real,distance_imaginary);


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% BASIC example 1.3
% A simple line segment, a simple query, postive distance in both front
% and rear segments
fig_num = 10012;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

points = [2 3];
pathXY = [0 0;3 0; 5 2; 8 2];

flag_snap_type = 1;

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length,distance_real,distance_imaginary] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);
fprintf(1,['Figure: %d,\n\t\t Closest point is: %.2f %.2f \n' ...
    '\t\t Matched to the path segment given by indices %d and %d, \n' ...
    '\t\t S-coordinate is: %.2f, \n' ...
    '\t\t percent_along_length is: %.2f\n' ...
    '\t\t real distance is: %.2f\n, ' ...
    '\t\t imag distance is %.2f\n, '],...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length,...
    distance_real,distance_imaginary);


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example 2
% A simple line segment, a simple query, negative distance purely in rear
% segment
fig_num = 10013;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

points = [0.5 -1];
pathXY = [0 0;3 0; 5 2; 8 2];

flag_snap_type = 1;

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example 3
% A simple line segment, a pre-start query, positive distance
fig_num = 10014;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

points = [-0.5 1];
pathXY = [0 0;3 0; 5 2; 8 2];

flag_snap_type = 1;

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example 4
% A simple line segment, a pre-start query, negative distance
fig_num = 10015;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

points = [-0.5 -0.2];
pathXY = [0 0;3 0; 5 2; 8 2];

flag_snap_type = 1;

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example 5
% A simple line segment, a post-end query, positive distance
fig_num = 10016;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

points = [9 3];
pathXY = [0 0;3 0; 5 2; 8 2];

flag_snap_type = 1;

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example 6
% A simple line segment, a post-end query, positive distance
fig_num = 10017;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

points = [9 1];
pathXY = [0 0;3 0; 5 2; 8 2];

flag_snap_type = 1;

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);



% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example 7.1
% A right angle, at midpoint in angle
fig_num = 10018;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

points = [1 1];
pathXY = [1 0;0 0; 0 1];

flag_snap_type = 1;

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example 7.2
% A right angle, at midpoint in angle
fig_num = 10019;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

points = [0.5 0.5];
pathXY = [1 0;0 0; 0 1];

flag_snap_type = 1;

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example 7.2
% A right angle, at nudge from midpoint in angle
fig_num = 10020;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

points = [0.40 0.5];
pathXY = [1 0;0 0; 0 1];

flag_snap_type = 1;

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example 7.3
% A right angle, at nudge from midpoint in angle
fig_num = 10021;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

points = [0.40 0.5];
pathXY = [1 0;0 0; 0 1];

flag_snap_type = 1;

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Flag Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   ______ _               _______        _       
%  |  ____| |             |__   __|      | |      
%  | |__  | | __ _  __ _     | | ___  ___| |_ ___ 
%  |  __| | |/ _` |/ _` |    | |/ _ \/ __| __/ __|
%  | |    | | (_| | (_| |    | |  __/\__ \ |_\__ \
%  |_|    |_|\__,_|\__, |    |_|\___||___/\__|___/
%                   __/ |                         
%                  |___/                          
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Flag%20Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§
close all

%% BASIC example 20001
fig_num = 20001;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Tests the flag_snap_type cases with an outside corner point
points = [2.5 1.3];
pathXY = [1 1; 2 1; 2 0]; % Define an XY path
flags_snap_type = [1; 2; 3];

for ith_flag = 1:length(flags_snap_type(:,1))
    flag_snap_type = flags_snap_type(ith_flag,:);
    temp_fig_num = fig_num+ith_flag-1;
    figure(temp_fig_num); clf;

    [closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
        fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,temp_fig_num);
    title(sprintf('Flag type: %.0d',flag_snap_type));

    axis([-1 4 -1 3]);
    fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
        temp_fig_num, closest_path_point(1,1),closest_path_point(1,2),...
        first_path_point_index,second_path_point_index, ...
        s_coordinate, percent_along_length);
end


%% BASIC example 
% Tests the flag_snap_type = 1 case
fig_num = 30001;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

points = [0.5 1.5; 1 1.5; 1.5 1.5; 2 1.5; 2.5 2; 2.5 1.5; 2.5 1; 2.5 0.5; 2.5 0; 2.5 -0.5; 2 -0.5; 1.5 -0.5; 1.5 0; 1.5 0.5; 1 0.5; 0.5 0.5; 0.5 1];
pathXY = [1 1; 2 1; 2 0]; % Define an XY path
flag_snap_type = 1;


for ith_point = 1:length(points(:,1))
    point = points(ith_point,:);
    figure(fig_num); clf;
    [closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
        fcn_Path_snapPointToPathViaVectors(point, pathXY,flag_snap_type,fig_num);
    axis([-1 4 -1 3]);
    fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
        fig_num, closest_path_point(1,1),closest_path_point(1,2),...
        first_path_point_index,second_path_point_index, ...
        s_coordinate, percent_along_length);
end

% Now, call all the points at once to show the function is vectorized
fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num+10);

%% BASIC example 1.02
% Tests the flag_snap_type = 2 case
fig_num = 30002;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

points = [0.5 1.5; 1 1.5; 1.5 1.5; 2 1.5; 2.5 2; 2.5 1.5; 2.5 1; 2.5 0.5; 2.5 0; 2.5 -0.5; 2 -0.5; 1.5 -0.5; 1.5 0; 1.5 0.5; 1 0.5; 0.5 0.5; 0.5 1];
pathXY = [1 1; 2 1; 2 0]; % Define an XY path
flag_snap_type = 2;


for ith_point = 1:length(points(:,1))
    point = points(ith_point,:);
    [closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
        fcn_Path_snapPointToPathViaVectors(point, pathXY,flag_snap_type,fig_num);
    axis([-1 4 -1 3]);
    fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
        fig_num, closest_path_point(1,1),closest_path_point(1,2),...
        first_path_point_index,second_path_point_index, ...
        s_coordinate, percent_along_length);
end


%% BASIC example 1.03
fig_num = 30003;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Tests the flag_snap_type = 3 case
points = [0.5 1.5; 1 1.5; 1.5 1.5; 2 1.5; 2.5 2; 2.5 1.5; 2.5 1; 2.5 0.5; 2.5 0; 2.5 -0.5; 2 -0.5; 1.5 -0.5; 1.5 0; 1.5 0.5; 1 0.5; 0.5 0.5; 0.5 1];
pathXY = [1 1; 2 1; 2 0]; % Define an XY path
flag_snap_type = 3;


for ith_point = 1:length(points(:,1))
    point = points(ith_point,:);
    [closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
        fcn_Path_snapPointToPathViaVectors(point, pathXY,flag_snap_type,fig_num);
    axis([-1 4 -1 3]);
    fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
        fig_num, closest_path_point(1,1),closest_path_point(1,2),...
        first_path_point_index,second_path_point_index, ...
        s_coordinate, percent_along_length);
end


% %% BASIC example 1.04
% % Tests the flag_snap_type = 4 case
% points = [0.5 1.5; 1 1.5; 1.5 1.5; 2 1.5; 2.5 2; 2.5 1.5; 2.5 1; 2.5 0.5; 2.5 0; 2.5 -0.5; 2 -0.5; 1.5 -0.5; 1.5 0; 1.5 0.5; 1 0.5; 0.5 0.5; 0.5 1];
% pathXY = [1 1; 2 1; 2 0]; % Define an XY path
% flag_snap_type = 4;
% 
% 
% for ith_point = 1:length(points(:,1))
%     point = points(ith_point,:);
%     fignum = 11104;
%     [closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
%         fcn_Path_snapPointToPathViaVectors(point, pathXY,flag_snap_type,fignum);
%     axis([-1 4 -1 3]);
%     fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
%         fignum, closest_path_point(1,1),closest_path_point(1,2),...
%         first_path_point_index,second_path_point_index, ...
%         s_coordinate, percent_along_length);
% end

%%

close all;

%% BASIC example 1.01
fig_num = 40001;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Tests the before start case
flag_snap_type = 1;
points = [0.4 0.2];
pathXY = [0.5 0.2; 0.9 0.9; 1 0.4]; % Define an XY path

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

%% BASIC example 1.02
% Tests the after end case
fig_num = 40002;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

flag_snap_type = 1;
points = [1.1 0.4];
pathXY = [0.5 0.2; 0.9 0.9; 1 0.4]; % Define an XY path

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);


%% BASIC example 1.03
% Tests the after back and before front case
fig_num = 40003;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;


flag_snap_type = 1;
points = [0.9 1];
pathXY = [0.5 0.2; 0.9 0.9; 1 0.4]; % Define an XY path

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);


%% BASIC example 1.04
% Tests the after back and in front case
fig_num = 40004;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

flag_snap_type = 1;
points = [1.2 0.9];
pathXY = [0.5 0.2; 0.9 0.9; 1 0.4]; % Define an XY path

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

%% BASIC example 1.05
% Tests the in back and before front case
fig_num = 40005;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

flag_snap_type = 1;
points = [0.6 0.9];
pathXY = [0.5 0.2; 0.9 0.9; 1 0.4]; % Define an XY path

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);


%% BASIC example 1.06
% Tests the in back and in front case, where back wins
fig_num = 40006;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

flag_snap_type = 1;
points = [0.8 0.7];
pathXY = [0.5 0.2; 0.9 0.9; 1 0.4]; % Define an XY path

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);


%% BASIC example 1.07
% Tests the in back and in front case, where front wins
fig_num = 40007;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

flag_snap_type = 1;
points = [0.9 0.7];
pathXY = [0.5 0.2; 0.9 0.9; 1 0.4]; % Define an XY path

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

%% BASIC example 1.01
fig_num = 40008;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

flag_snap_type = 1;
points = [0.5 0.5];
pathXY = [0.5 0.2; 0.9 0.9; 1.5 0.6]; % Define an XY path

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

%% BASIC example 1.2 - works
fig_num = 40009;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

flag_snap_type = 1;
points = [1.4 1.3]; % Define the query point
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 1.5 0.6; 3 0]; % Define an XY path

% Snap the point onto the path
[closest_path_point,s_coordinate,...
    first_path_point_index,second_path_point_index,...
    percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);

% Print results to the workspace
fprintf(1,'Figure: %d\n',fig_num);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);

%% BASIC example 1.21 - works
fig_num = 40010;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

flag_snap_type = 1;
points = [0.5 0.75]; % Define the query point
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 1.5 0.6; 3 0]; % Define an XY path

% Snap the point onto the path
[closest_path_point,s_coordinate,...
    first_path_point_index,second_path_point_index,...
    percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);

% Print results to the workspace
fprintf(1,'Figure: %d\n',fig_num);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);


%% BASIC example 1.3 - works
fig_num = 40011;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

flag_snap_type = 1;
points = [1.5 1];
% pathXY = [0 0; 0.5 0.2; 0.9 0.9; 3 0];
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 3 0];

% Snap the point onto the path
[closest_path_point,s_coordinate,...
    first_path_point_index,second_path_point_index,...
    percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);

% Print results to the workspace
fprintf(1,'Figure: %d\n',fig_num);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);


%% BASIC example 1.4 - works, but shows that it is on neither segement
fig_num = 40012;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

flag_snap_type = 1;
points = [0.9 1.4]; 
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 1.5 0.6; 3 0]; % Define an XY path

% Snap the point onto the path
[closest_path_point,s_coordinate,...
    first_path_point_index,second_path_point_index,...
    percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);

% Print results to the workspace
fprintf(1,'Figure: %d\n',fig_num);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);

%% BASIC example 1.5 - works, but is on BOTH segments
fig_num = 40013;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

flag_snap_type = 1;
points = [1 0.5]; 
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 1.5 0.6; 3 0]; % Define an XY path

% Snap the point onto the path
[closest_path_point,s_coordinate,...
    first_path_point_index,second_path_point_index,...
    percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);

% Print results to the workspace
fprintf(1,'Figure: %d\n',fig_num);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);

%% BASIC example 1.5 - works, but is on BOTH segments
fig_num = 40014;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

flag_snap_type = 1;
points = [1 0.5]; 
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 1.5 0.6; 3 0]; % Define an XY path

% Snap the point onto the path
[closest_path_point,s_coordinate,...
    first_path_point_index,second_path_point_index,...
    percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);

% Print results to the workspace
fprintf(1,'Figure: %d\n',fig_num);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);

%% BASIC example 1.6 - all points
fig_num = 40015;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

flag_snap_type = 1;
%pathXY = [0 0; 0.5 0.2; 0.9 0.9; 1.5 0.6; 3 0]; % Define an XY path
pathXY = [0.5 0.2; 0.9 0.9; 1.5 0.6]; % Define an XY path
% pathXY = [0 0; 1 1]; 

x_range = linspace(-4,4,100);
y_range = linspace(-0.5,1.5,100);

[X,Y] = meshgrid(x_range,y_range);
Z = zeros(size(X));
S = Z;
first_point = Z;
second_point = Z;
Percent = Z;

for ith_row = 1:length(X(:,1))
    for jth_col = 1:length(X(1,:))
        query_point = [X(ith_row,jth_col) Y(ith_row,jth_col)];
        [closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
            fcn_Path_snapPointToPathViaVectors(query_point, pathXY);
        
        distance = sum((query_point-closest_path_point).^2,2).^0.5;
        Z(ith_row,jth_col) = distance;
        S(ith_row,jth_col) = s_coordinate;
        first_point(ith_row,jth_col) = first_path_point_index;
        second_point(ith_row,jth_col) = second_path_point_index;
        Percent(ith_row,jth_col) = percent_along_length;
    end
end
figure(fig_num); 
clf;
hold on;
surfc(X,Y,Z-0.2);
plot(pathXY(:,1),pathXY(:,2),'r-','LineWidth',3);
view(30,50)
title('Distance');


figure(fig_num*10); 
clf;
hold on;
surfc(X,Y,S);
plot(pathXY(:,1),pathXY(:,2),'r-','LineWidth',3);
view(30,50)
title('Station');


figure(fig_num*10+1); 
clf;
hold on;
surfc(X,Y,first_point);
plot(pathXY(:,1),pathXY(:,2),'r-','LineWidth',3);
view(30,50)
title('First point');


figure(fig_num*10+2); 
clf;
hold on;
surfc(X,Y,second_point);
plot(pathXY(:,1),pathXY(:,2),'r-','LineWidth',3);
view(30,50)
title('Second point');

figure(fig_num*10+3); 
clf;
hold on;
surfc(X,Y,Percent);
plot(pathXY(:,1),pathXY(:,2),'r-','LineWidth',3);
view(30,50)
title('Percent along path');

%%

close all

%% BASIC example 3 - positive s-coords (after path ends)
fig_num = 50001;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;


flag_snap_type = 1;

points = [4 0.2];
pathXY = [0 0; 1 0; 2 0];

% [closest_path_point,s_coordinate] = ...
%     fcn_Path_snapPointToPathViaVectors(point, path,flag_snap_type,fignum);
% fprintf(1,'Figure: %d, Closest point is: %.2f %.2f, S-coordinate is: %.2f \n',fignum, closest_path_point(1,1),closest_path_point(1,2), s_coordinate);

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

%% BASIC example 4 - an example of percentage along segment greater than 100% even though "inside" path
fig_num = 50002;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

points = [0.8 1.3];
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 3 0];
flag_snap_type = 1;

% [closest_path_point,s_coordinate] = ...
%     fcn_Path_snapPointToPathViaVectors(point, path,flag_snap_type,fignum);
% fprintf(1,'Figure: %d, Closest point is: %.2f %.2f, S-coordinate is: %.2f \n',fignum, closest_path_point(1,1),closest_path_point(1,2), s_coordinate);

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,fig_num);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);



%% ADVANCED example
fig_num = 60001;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Create some paths:
path1 = [
    8.1797    8.6006
    8.4101   15.8892
    10.0230   25.5102
    10.4839   39.7959
    10.7143   50.8746
    10.9447   61.9534
    13.0184   73.6152
    16.7051   82.9446
    24.7696   89.3586
    35.3687   91.1079
    42.5115   91.3994
    52.8802   89.3586
    58.1797   85.2770
    60.9447   78.5714
    57.9493   72.1574
    57.2581   63.7026
    59.1014   58.1633
    63.0184   57.2886
    67.3963   56.9971
    69.9309   56.7055
    74.3088   56.1224
    78.6866   54.0816
    80.9908   51.4577
    82.3733   49.1254
    84.6774   40.6706
    84.6774   34.5481
    82.8341   28.7172
    80.0691   26.9679
    76.3825   25.2187
    69.2396   20.2624
    65.5530   18.2216
    60.9447   18.8047
    57.4885   22.0117
    50.5760   28.1341
    47.1198   30.7580
    43.8940   34.8397
    39.7465   37.7551
    37.6728   40.9621
    32.6037   42.7114
    30.0691   43.0029
    28.2258   43.2945
    26.3825   43.2945
    24.5392   44.4606
    20.8525   47.3761
    19.2396   49.7085
    16.7051   53.2070
    14.1705   58.1633
    13.0184   62.2449
    10.0230   70.1166
    8.6406   74.4898
    7.7189   79.7376
    6.5668   82.9446
    5.1843   86.7347
    4.2627   88.4840
    3.8018   89.0671];

figure(111); plot(path1(:,1),path1(:,2),'r-o');
text(path1(1,1),path1(1,2),'Start');

% Create a query
points = [75 45];
flag_snap_type = 1;
[closest_path_point,s_coordinate] = ...
    fcn_Path_snapPointToPathViaVectors(points, path1,flag_snap_type,fig_num);
fprintf(1,'Figure: %d, Closest point is: %.2f %.2f, S-coordinate is: %.2f \n',fig_num, closest_path_point(1,1),closest_path_point(1,2), s_coordinate);


%% ADVANCED example - more complicated path
fig_num = 60002;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

npoints = 20;
Ntests = 10;
seeds{Ntests} = 0;
flag_snap_type = 1;
for i_test = 1:Ntests
    seeds{i_test} = rng;    
    rng(seeds{i_test});
    rand_x = rand(npoints,1);
    rand_y = rand(npoints,1);
    
    pathXY = [cumsum(rand_x),cumsum(rand_y)];
    points = mean(pathXY,1);
    
    temp_fig_num = fig_num+i_test-1;
    [closest_path_point,s_coordinate] = ...
        fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,temp_fig_num);
    fprintf(1,'Figure: %d, Closest point is: %.2f %.2f, S-coordinate is: %.2f \n',temp_fig_num, closest_path_point(1,1),closest_path_point(1,2), s_coordinate);
end


%% Fast Mode Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ______        _     __  __           _        _______        _
% |  ____|      | |   |  \/  |         | |      |__   __|      | |
% | |__ __ _ ___| |_  | \  / | ___   __| | ___     | | ___  ___| |_ ___
% |  __/ _` / __| __| | |\/| |/ _ \ / _` |/ _ \    | |/ _ \/ __| __/ __|
% | | | (_| \__ \ |_  | |  | | (_) | (_| |  __/    | |  __/\__ \ |_\__ \
% |_|  \__,_|___/\__| |_|  |_|\___/ \__,_|\___|    |_|\___||___/\__|___/
%
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Fast%20Mode%20Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures start with 8

close all;
fprintf(1,'Figure: 8XXXXXX: Demo of fast mode cases\n');

%% Basic example - NO FIGURE
fig_num = 80001;
fprintf(1,'Figure: %.0f: Demo of fast mode, empty fig_num\n',fig_num);
figure(fig_num); close(fig_num);

points = [0.5 0.5];
pathXY = [0 0; 2 2];
flag_snap_type = 1;

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length,distance_real,distance_imaginary] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,[]);


% Make sure function worked
true_closest_path_point = [0.5 0.5];
true_s_coordinate = 0.5*(2^0.5);
true_first_path_point_index = 1;
true_second_path_point_index = 2;
true_percent_along_length = 0.25;
true_distance_real = 0;
true_distance_imaginary = 0;

assert(isequal(round(closest_path_point,4),round(true_closest_path_point,4)));
assert(isequal(round(s_coordinate,4),round(true_s_coordinate,4)));
assert(isequal(round(first_path_point_index,4),round(true_first_path_point_index,4)));
assert(isequal(round(second_path_point_index,4),round(true_second_path_point_index,4)));
assert(isequal(round(percent_along_length,4),round(true_percent_along_length,4)));
assert(isequal(round(distance_real,4),round(true_distance_real,4)));
assert(isequal(round(distance_imaginary,4),round(true_distance_imaginary,4)));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

points = [0.5 0.5];
pathXY = [0 0; 2 2];
flag_snap_type = 1;

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length,distance_real,distance_imaginary] = ...
    fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,-1);


% Make sure function worked
true_closest_path_point = [0.5 0.5];
true_s_coordinate = 0.5*(2^0.5);
true_first_path_point_index = 1;
true_second_path_point_index = 2;
true_percent_along_length = 0.25;
true_distance_real = 0;
true_distance_imaginary = 0;

assert(isequal(round(closest_path_point,4),round(true_closest_path_point,4)));
assert(isequal(round(s_coordinate,4),round(true_s_coordinate,4)));
assert(isequal(round(first_path_point_index,4),round(true_first_path_point_index,4)));
assert(isequal(round(second_path_point_index,4),round(true_second_path_point_index,4)));
assert(isequal(round(percent_along_length,4),round(true_percent_along_length,4)));
assert(isequal(round(distance_real,4),round(true_distance_real,4)));
assert(isequal(round(distance_imaginary,4),round(true_distance_imaginary,4)));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

% Fill in sample paths (as a starter)
points = [0.5 0.5];
pathXY = [0 0; 2 2];
flag_snap_type = 1;

Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length,distance_real,distance_imaginary] = ...
        fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,[]);
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length,distance_real,distance_imaginary] = ...
        fcn_Path_snapPointToPathViaVectors(points, pathXY,flag_snap_type,-1);
end
fast_method = toc;

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

% Plot results as bar chart
figure(373737);
clf;
hold on;

X = categorical({'Normal mode','Fast mode'});
X = reordercats(X,{'Normal mode','Fast mode'}); % Forces bars to appear in this exact order, not alphabetized
Y = [slow_method fast_method ]*1000/Niterations;
bar(X,Y)
ylabel('Execution time (Milliseconds)')


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% BUG cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ____  _    _  _____
% |  _ \| |  | |/ ____|
% | |_) | |  | | |  __    ___ __ _ ___  ___  ___
% |  _ <| |  | | | |_ |  / __/ _` / __|/ _ \/ __|
% | |_) | |__| | |__| | | (_| (_| \__ \  __/\__ \
% |____/ \____/ \_____|  \___\__,_|___/\___||___/
%
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=BUG%20cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All bug case figures start with the number 9

% close all;

%% BUG



%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§

