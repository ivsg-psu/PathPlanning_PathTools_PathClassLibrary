% script_test_fcn_Path_snapPointOntoNearestPath.m
% This is a script to exercise the function: fcn_Path_snapPointOntoNearestPath.m
% This function was written on 2020_10_10 by S. Brennan, sbrennan@psu.edu

% Revision history:
%     2020_10_10 
%     -- first write of the code
%     2020_11_10 
%     -- changed function names in prep for DataClean class
%     2020_12_03 
%     -- updated some of the plotting/debug details to improve
%     2020_11_12 
%     -- modified to prep for Path class
%     2021_01_08 
%     -- cleaned up comments
%     -- added argument checking
%     -- cleaned up notation to show path vs pathSXY
%     2021_01_09
%     -- completely converted code to path form, not pathSXY
%     2021_03_21
%     -- modified to allow 3D snapping
%     -- changed input checks to include 3D paths
%     2023_06_03
%     -- modified to allow vectorized snapping
%     2023_08_08
%     -- bug fixes for points before start or after end


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

%% BASIC example 1
% A simple line segment, a simple query, positive distance in both front
% and rear segments
point = [0.5 2];
pathXY = [0 0;3 0; 5 2; 8 2];
flag_snap_type = 1;

fignum = 111;
[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length,distance_real,distance_imaginary] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);
fprintf(1,['Figure: %d,\n\t\t Closest point is: %.2f %.2f \n' ...
    '\t\t Matched to the path segment given by indices %d and %d, \n' ...
    '\t\t S-coordinate is: %.2f, \n' ...
    '\t\t percent_along_length is: %.2f\n' ...
    '\t\t real distance is: %.2f\n, ' ...
    '\t\t imag distance is %.2f\n, '],...
    fignum, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length,...
    distance_real,distance_imaginary);

%% BASIC example 1.1
% A simple line segment, a simple query, negative distance in both front
% and rear segments
point = [5 0.5];
pathXY = [0 0;3 0; 5 2; 8 2];

flag_snap_type = 1;

fignum = 111;
[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length,distance_real,distance_imaginary] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);
fprintf(1,['Figure: %d,\n\t\t Closest point is: %.2f %.2f \n' ...
    '\t\t Matched to the path segment given by indices %d and %d, \n' ...
    '\t\t S-coordinate is: %.2f, \n' ...
    '\t\t percent_along_length is: %.2f\n' ...
    '\t\t real distance is: %.2f\n, ' ...
    '\t\t imag distance is %.2f\n, '],...
    fignum, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length,...
    distance_real,distance_imaginary);

%% BASIC example 1.2
% A simple line segment, a simple query, negative distance in both front
% and rear segments
point = [4 1];
pathXY = [0 0;3 0; 5 2; 8 2];

flag_snap_type = 1;

fignum = 111;
[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length,distance_real,distance_imaginary] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);
fprintf(1,['Figure: %d,\n\t\t Closest point is: %.2f %.2f \n' ...
    '\t\t Matched to the path segment given by indices %d and %d, \n' ...
    '\t\t S-coordinate is: %.2f, \n' ...
    '\t\t percent_along_length is: %.2f\n' ...
    '\t\t real distance is: %.2f\n, ' ...
    '\t\t imag distance is %.2f\n, '],...
    fignum, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length,...
    distance_real,distance_imaginary);


%% BASIC example 2
% A simple line segment, a simple query, negative distance purely in rear
% segment
point = [0.5 -1];
pathXY = [0 0;3 0; 5 2; 8 2];

fignum = 222;
[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fignum, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

%% BASIC example 3
% A simple line segment, a pre-start query, positive distance
point = [-0.5 1];
pathXY = [0 0;3 0; 5 2; 8 2];

fignum = 333;
[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fignum, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

%% BASIC example 4
% A simple line segment, a pre-start query, negative distance
point = [-0.5 -0.2];
pathXY = [0 0;3 0; 5 2; 8 2];


fignum = 444;
[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fignum, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

%% BASIC example 5
% A simple line segment, a post-end query, positive distance
point = [1 2];
pathXY = [0 0;3 0; 5 2; 8 2];



fignum = 555;
[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fignum, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);
xlim([-1 3])

%% BASIC example 6
% A simple line segment, a post-end query, positive distance
point = [2 1];
pathXY = [0 0;3 0; 5 2; 8 2];


fignum = 666;
[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fignum, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);
xlim([-1 3])


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

%% BASIC example 1.01
% Tests the flag_snap_type = 1 case
points = [0.5 1.5; 1 1.5; 1.5 1.5; 2 1.5; 2.5 2; 2.5 1.5; 2.5 1; 2.5 0.5; 2.5 0; 2.5 -0.5; 2 -0.5; 1.5 -0.5; 1.5 0; 1.5 0.5; 1 0.5; 0.5 0.5; 0.5 1];
pathXY = [1 1; 2 1; 2 0]; % Define an XY path
flag_snap_type = 1;


for ith_point = 1:length(points(:,1))
    point = points(ith_point,:);
    fignum = 11100+ith_point;
    [closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
        fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);
    axis([-1 4 -1 3]);
    fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
        fignum, closest_path_point(1,1),closest_path_point(1,2),...
        first_path_point_index,second_path_point_index, ...
        s_coordinate, percent_along_length);
end

% Now, call all the points at once to show the function is vectorized
fig_num = 111111;
fcn_Path_snapPointOntoNearestPath(points, pathXY,flag_snap_type,fignum);

%% BASIC example 1.02
% Tests the flag_snap_type = 2 case
points = [0.5 1.5; 1 1.5; 1.5 1.5; 2 1.5; 2.5 2; 2.5 1.5; 2.5 1; 2.5 0.5; 2.5 0; 2.5 -0.5; 2 -0.5; 1.5 -0.5; 1.5 0; 1.5 0.5; 1 0.5; 0.5 0.5; 0.5 1];
pathXY = [1 1; 2 1; 2 0]; % Define an XY path
flag_snap_type = 2;


for ith_point = 1:length(points(:,1))
    point = points(ith_point,:);
    fignum = 11100+ith_point;
    [closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
        fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);
    axis([-1 4 -1 3]);
    fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
        fignum, closest_path_point(1,1),closest_path_point(1,2),...
        first_path_point_index,second_path_point_index, ...
        s_coordinate, percent_along_length);
end


%% BASIC example 1.03
% Tests the flag_snap_type = 3 case
points = [0.5 1.5; 1 1.5; 1.5 1.5; 2 1.5; 2.5 2; 2.5 1.5; 2.5 1; 2.5 0.5; 2.5 0; 2.5 -0.5; 2 -0.5; 1.5 -0.5; 1.5 0; 1.5 0.5; 1 0.5; 0.5 0.5; 0.5 1];
pathXY = [1 1; 2 1; 2 0]; % Define an XY path
flag_snap_type = 3;


for ith_point = 1:length(points(:,1))
    point = points(ith_point,:);
    fignum = 11100+ith_point;
    [closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
        fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);
    axis([-1 4 -1 3]);
    fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
        fignum, closest_path_point(1,1),closest_path_point(1,2),...
        first_path_point_index,second_path_point_index, ...
        s_coordinate, percent_along_length);
end


%% BASIC example 1.04
% Tests the flag_snap_type = 4 case
points = [0.5 1.5; 1 1.5; 1.5 1.5; 2 1.5; 2.5 2; 2.5 1.5; 2.5 1; 2.5 0.5; 2.5 0; 2.5 -0.5; 2 -0.5; 1.5 -0.5; 1.5 0; 1.5 0.5; 1 0.5; 0.5 0.5; 0.5 1];
pathXY = [1 1; 2 1; 2 0]; % Define an XY path
flag_snap_type = 4;


for ith_point = 1:length(points(:,1))
    point = points(ith_point,:);
    fignum = 11100+ith_point;
    [closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
        fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);
    axis([-1 4 -1 3]);
    fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
        fignum, closest_path_point(1,1),closest_path_point(1,2),...
        first_path_point_index,second_path_point_index, ...
        s_coordinate, percent_along_length);
end


%% BASIC example 1.01
% Tests the before start case
point = [0.4 0.2];
pathXY = [0.5 0.2; 0.9 0.9; 1 0.4]; % Define an XY path

fignum = 11101;
[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fignum, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

%% BASIC example 1.02
% Tests the after end case
point = [1.1 0.4];
pathXY = [0.5 0.2; 0.9 0.9; 1 0.4]; % Define an XY path

fignum = 11102;
[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fignum, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);


%% BASIC example 1.03
% Tests the after back and before front case
point = [0.9 1];
pathXY = [0.5 0.2; 0.9 0.9; 1 0.4]; % Define an XY path

fignum = 11103;
[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fignum, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);


%% BASIC example 1.04
% Tests the after back and in front case
point = [1.2 0.9];
pathXY = [0.5 0.2; 0.9 0.9; 1 0.4]; % Define an XY path

fignum = 11104;
[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fignum, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

%% BASIC example 1.05
% Tests the in back and before front case
point = [0.6 0.9];
pathXY = [0.5 0.2; 0.9 0.9; 1 0.4]; % Define an XY path

fignum = 11105;
[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fignum, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);


%% BASIC example 1.06
% Tests the in back and in front case, where back wins
point = [0.8 0.7];
pathXY = [0.5 0.2; 0.9 0.9; 1 0.4]; % Define an XY path

fignum = 11106;
[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fignum, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);


%% BASIC example 1.07
% Tests the in back and in front case, where front wins
point = [0.9 0.7];
pathXY = [0.5 0.2; 0.9 0.9; 1 0.4]; % Define an XY path

fignum = 11106;
[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fignum, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);











%% BASIC example 1.01
point = [0.5 0.5];
pathXY = [0.5 0.2; 0.9 0.9; 1.5 0.6]; % Define an XY path

fignum = 11101;
[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fignum, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

%% BASIC example 1.2 - works
point = [1.4 1.3]; % Define the query point
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 1.5 0.6; 3 0]; % Define an XY path
fignum = 112; % Define the figure number

% Snap the point onto the path
[closest_path_point,s_coordinate,...
    first_path_point_index,second_path_point_index,...
    percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);

% Print results to the workspace
fprintf(1,'Figure: %d\n',fignum);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);

%% BASIC example 1.21 - works
point = [0.5 0.75]; % Define the query point
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 1.5 0.6; 3 0]; % Define an XY path
fignum = 112; % Define the figure number

% Snap the point onto the path
[closest_path_point,s_coordinate,...
    first_path_point_index,second_path_point_index,...
    percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);

% Print results to the workspace
fprintf(1,'Figure: %d\n',fignum);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);


%% BASIC example 1.3 - works
point = [1.5 1];
% pathXY = [0 0; 0.5 0.2; 0.9 0.9; 3 0];
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 3 0];

fignum = 113;

% Snap the point onto the path
[closest_path_point,s_coordinate,...
    first_path_point_index,second_path_point_index,...
    percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);

% Print results to the workspace
fprintf(1,'Figure: %d\n',fignum);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);


%% BASIC example 1.4 - works, but shows that it is on neither segement
point = [0.9 1.4]; 
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 1.5 0.6; 3 0]; % Define an XY path

fignum = 114;

% Snap the point onto the path
[closest_path_point,s_coordinate,...
    first_path_point_index,second_path_point_index,...
    percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);

% Print results to the workspace
fprintf(1,'Figure: %d\n',fignum);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);

%% BASIC example 1.5 - works, but is on BOTH segments
point = [1 0.5]; 
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 1.5 0.6; 3 0]; % Define an XY path

fignum = 115;

% Snap the point onto the path
[closest_path_point,s_coordinate,...
    first_path_point_index,second_path_point_index,...
    percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);

% Print results to the workspace
fprintf(1,'Figure: %d\n',fignum);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);

%% BASIC example 1.5 - works, but is on BOTH segments
point = [1 0.5]; 
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 1.5 0.6; 3 0]; % Define an XY path

fignum = 115;

% Snap the point onto the path
[closest_path_point,s_coordinate,...
    first_path_point_index,second_path_point_index,...
    percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);

% Print results to the workspace
fprintf(1,'Figure: %d\n',fignum);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);

%% BASIC example 1.6 - all points
%pathXY = [0 0; 0.5 0.2; 0.9 0.9; 1.5 0.6; 3 0]; % Define an XY path
pathXY = [0.5 0.2; 0.9 0.9; 1.5 0.6]; % Define an XY path


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
            fcn_Path_snapPointOntoNearestPath(query_point, pathXY);
        
        distance = sum((query_point-closest_path_point).^2,2).^0.5;
        Z(ith_row,jth_col) = distance;
        S(ith_row,jth_col) = s_coordinate;
        first_point(ith_row,jth_col) = first_path_point_index;
        second_point(ith_row,jth_col) = second_path_point_index;
        Percent(ith_row,jth_col) = percent_along_length;
    end
end
figure(1161); 
clf;
hold on;
surfc(X,Y,Z-0.2);
plot(pathXY(:,1),pathXY(:,2),'r-','LineWidth',3);
view(30,50)
title('Distance');


figure(1162); 
clf;
hold on;
surfc(X,Y,S);
plot(pathXY(:,1),pathXY(:,2),'r-','LineWidth',3);
view(30,50)
title('Station');


figure(1163); 
clf;
hold on;
surfc(X,Y,first_point);
plot(pathXY(:,1),pathXY(:,2),'r-','LineWidth',3);
view(30,50)
title('First point');


figure(1164); 
clf;
hold on;
surfc(X,Y,second_point);
plot(pathXY(:,1),pathXY(:,2),'r-','LineWidth',3);
view(30,50)
title('Second point');

figure(1165); 
clf;
hold on;
surfc(X,Y,Percent);
plot(pathXY(:,1),pathXY(:,2),'r-','LineWidth',3);
view(30,50)
title('Percent along path');

%% BASIC example 3 - positive s-coords (after path ends)
point = [4 0.2];
pathXY = [0 0; 1 0; 2 0];

fignum = 333;
% [closest_path_point,s_coordinate] = ...
%     fcn_Path_snapPointOntoNearestPath(point, path,flag_snap_type,fignum);
% fprintf(1,'Figure: %d, Closest point is: %.2f %.2f, S-coordinate is: %.2f \n',fignum, closest_path_point(1,1),closest_path_point(1,2), s_coordinate);

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fignum, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

%% BASIC example 4 - an example of percentage along segment greater than 100% even though "inside" path
point = [0.8 1.3];
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 3 0];

fignum = 444;
% [closest_path_point,s_coordinate] = ...
%     fcn_Path_snapPointOntoNearestPath(point, path,flag_snap_type,fignum);
% fprintf(1,'Figure: %d, Closest point is: %.2f %.2f, S-coordinate is: %.2f \n',fignum, closest_path_point(1,1),closest_path_point(1,2), s_coordinate);

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fignum, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ____  _____       _____                     _______        _       
%  |___ \|  __ \     / ____|                   |__   __|      | |      
%    __) | |  | |   | (___  _ __   __ _ _ __      | | ___  ___| |_ ___ 
%   |__ <| |  | |    \___ \| '_ \ / _` | '_ \     | |/ _ \/ __| __/ __|
%   ___) | |__| |    ____) | | | | (_| | |_) |    | |  __/\__ \ |_\__ \
%  |____/|_____/    |_____/|_| |_|\__,_| .__/     |_|\___||___/\__|___/
%                                      | |                             
%                                      |_|                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% BASIC example 3D - simple 3D snapping onto a vertex
point = [0.8 1.3 2.1];
pathXYZ = [0 0 0; 0.5 0.2 0.4; 0.9 0.9 0.8; 3 0 1];

fignum = 3331;

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXYZ,flag_snap_type,fignum);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fignum, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

%% BASIC example 3D - simple 3D snapping onto a vertex
point = [2 1.3 2.1];
pathXYZ = [0 0 0; 0.5 0.2 0.4; 0.9 0.9 0.8; 3 0 1];

fignum = 3332;

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXYZ,flag_snap_type,fignum);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fignum, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

%% ADVANCED example
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
fignum = 2222;
point = [75 45];
[closest_path_point,s_coordinate] = ...
    fcn_Path_snapPointOntoNearestPath(point, path1,flag_snap_type,fignum);
fprintf(1,'Figure: %d, Closest point is: %.2f %.2f, S-coordinate is: %.2f \n',fignum, closest_path_point(1,1),closest_path_point(1,2), s_coordinate);


%% ADVANCED example - more complicated path
npoints = 20;
Ntests = 10;
seeds{Ntests} = 0;
for i_test = 1:Ntests
    seeds{i_test} = rng;    
    rng(seeds{i_test});
    rand_x = rand(npoints,1);
    rand_y = rand(npoints,1);
    
    pathXY = [cumsum(rand_x),cumsum(rand_y)];
    point = mean(pathXY,1);
    
    fignum = i_test;
    [closest_path_point,s_coordinate] = ...
        fcn_Path_snapPointOntoNearestPath(point, pathXY,flag_snap_type,fignum);
    fprintf(1,'Figure: %d, Closest point is: %.2f %.2f, S-coordinate is: %.2f \n',fignum, closest_path_point(1,1),closest_path_point(1,2), s_coordinate);
end

