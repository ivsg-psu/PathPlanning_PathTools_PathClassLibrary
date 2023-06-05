% script_test_fcn_Path_findNearestPathPoints
% This is a script to exercise the function: fcn_Path_findNearestPathPoints.m
% This function was written on 2023_06_02 by S. Brennan in support of the
% fcn_Path_snapPointOntoNearestPath function expansion.
%
% Questions or comments? sbrennan@psu.edu 

% Revision history:   
% 2023_06_02 by sbrennan@psu.edu
% -- first write of the code


close all;

%% BASIC example 1
% Small path with single query point
query_points = [0.8 0.4];
pathXY = [0 0; 1 0; 2 0; 2 1];

fignum = 111;
closest_path_point_indicies = ...
    fcn_Path_findNearestPathPoints(query_points, pathXY,fignum);
fprintf(1,'Figure: %d,\n\t\t Closest point index is: %.0d %.2f\n',...
    fignum, closest_path_point_indicies(1,1));
title('fcn_Path_findNearestPathPoints: tested with a single point query in 2D','Interpreter','none');

%% BASIC example 2
% Small path with two query points
query_points = [0.8 0.4; 1.6 0.8];
pathXY = [0 0; 1 0; 2 0; 2 1];

fignum = 222;
closest_path_point_indicies = ...
    fcn_Path_findNearestPathPoints(query_points, pathXY,fignum);
title('fcn_Path_findNearestPathPoints: tested with two point queries in 2D','Interpreter','none');
fprintf(1,'Figure: %d,\n\t\t Closest point index is: %.0d %.2f\n',...
    fignum, closest_path_point_indicies(1,1));

%% BASIC example 3
% Small path with many query points
query_points = rand(20,2);
query_points(:,1)=query_points(:,1)*4-1;
query_points(:,2)=query_points(:,2)*3-1.5;

pathXY = [0 0; 1 0; 2 0; 2 1];

fignum = 333;
closest_path_point_indicies = ...
    fcn_Path_findNearestPathPoints(query_points, pathXY,fignum);
fprintf(1,'Figure: %d,\n\t\t Closest point index is: %.0d %.2f\n',...
    fignum, closest_path_point_indicies(1,1));
title('fcn_Path_findNearestPathPoints: tested with many point queries in 2D','Interpreter','none');



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
query_points = [0.8 1.3 2.1];
pathXYZ = [0 0 0; 0.5 0.2 0.4; 0.9 0.9 0.8; 3 0 1];

fignum = 3331;

closest_path_point_indicies = ...
    fcn_Path_findNearestPathPoints(query_points, pathXYZ,fignum);
title('fcn_Path_findNearestPathPoints: tested with single point queries in 3D','Interpreter','none');

fprintf(1,'Figure: %d, Closest point index is: %.0d \n',...
    fignum, closest_path_point_indicies(1,1));
view(3);


%% BASIC example 3D - simple 3D snapping onto a vertex
query_points = [2 1.3 2.1];
pathXYZ = [0 0 0; 0.5 0.2 0.4; 0.9 0.9 0.8; 3 0 1];

fignum = 3332;

closest_path_point_indicies = ...
    fcn_Path_findNearestPathPoints(query_points, pathXYZ,fignum);
title('fcn_Path_findNearestPathPoints: tested with single point queries in 3D','Interpreter','none');
xlabel('X'); ylabel('Y'); zlabel('Z');

fprintf(1,'Figure: %d, Closest point index is: %.0d \n',...
    fignum, closest_path_point_indicies(1,1));
view(3);

%% BASIC example 3D - simple 3D snapping with two query points
query_points = [0.8 1.3 2.1; 2 1.3 2.1];
pathXYZ = [0 0 0; 0.5 0.2 0.4; 0.9 0.9 0.8; 3 0 1];

fignum = 3333;

closest_path_point_indicies = ...
    fcn_Path_findNearestPathPoints(query_points, pathXYZ,fignum);
title('fcn_Path_findNearestPathPoints: tested with single point queries in 3D','Interpreter','none');
xlabel('X'); ylabel('Y'); zlabel('Z');

fprintf(1,'Figure: %d, Closest point index is: %.0d \n',...
    fignum, closest_path_point_indicies(1,1));
view(3);


%% BASIC example 3
% Small path with many query points
query_points = rand(20,3);
query_points(:,1)=query_points(:,1)*5-2.5;
query_points(:,2)=query_points(:,2)*3-1.5;
query_points(:,3)=query_points(:,2)*3-1.5;

pathXYZ = [0 0 0; 0.5 0.2 0.4; 0.9 0.9 0.8; 3 0 1];

fignum = 3333;

closest_path_point_indicies = ...
    fcn_Path_findNearestPathPoints(query_points, pathXYZ,fignum);
title('fcn_Path_findNearestPathPoints: tested with single point queries in 3D','Interpreter','none');
xlabel('X'); ylabel('Y'); zlabel('Z');


fprintf(1,'Figure: %d, Closest point index is: %.0d \n',...
    fignum, closest_path_point_indicies(1,1));
view(3);


