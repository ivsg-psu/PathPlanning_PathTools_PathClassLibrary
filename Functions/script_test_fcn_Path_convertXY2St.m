% script_test_fcn_Path_convertXY2St.m
% This is a script to exercise the function: fcn_Path_convertXY2St.m
% This function was written on 2023_08_26 by S. Brennan, sbrennan@psu.edu

% Revision history:
% 2023_08_26 by S. Brennan
% -- first write of the code



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
% A simple line segment, a simple query, zero distance in rear segments
XY_points = [0 1];
referencePath = [-1 0; 1 0];
flag_snap_type = 1;

fig_num = 111;
St_points = fcn_Path_convertXY2St(referencePath,XY_points, flag_snap_type,fig_num);
expected_solution = [1 1];
assert(abs(sum((St_points - expected_solution).^2,2))<1E-10);

%% BASIC example 2
% A simple line segment, a simple query, zero distance in rear segments
XY_points = [0 1];
referencePath = [0 0; 2 2];
flag_snap_type = 1;

fig_num = 222;
St_points = fcn_Path_convertXY2St(referencePath,XY_points, flag_snap_type,fig_num);
expected_solution = [1 1]./(2^0.5);
assert(abs(sum((St_points - expected_solution).^2,2))<1E-10);


%% BASIC example 3
% A simple line segment, a simple query, zero distance in rear segments
XY_points = [0 -1];
referencePath = [-1 0; 1 0];
flag_snap_type = 1;

fig_num = 333;
St_points = fcn_Path_convertXY2St(referencePath,XY_points, flag_snap_type,fig_num);
expected_solution = [1 -1];
assert(abs(sum((St_points - expected_solution).^2,2))<1E-10);

%% BASIC example 4
% A simple line segment, a complex number in rear segments
XY_points = [-2 -1];
referencePath = [-1 0; 1 0];
flag_snap_type = 1;

fig_num = 444;
St_points = fcn_Path_convertXY2St(referencePath,XY_points, flag_snap_type,fig_num);
expected_solution = [-1 -1+1i];
assert(abs(sum((St_points - expected_solution).^2,2))<1E-10);

%% BASIC example 5 - many points
% A 90-degree line segment with multiple surrounding queries
XY_points = [-2 1; -1 1; 0 1; 1 1; 2 1; 2 0; 2 -1; 2 -2; 1 -2; 0 -2; 0 -1; -1 -1; -2 -1; -2 0];
referencePath = [-1 0; 1 0; 1 -1];
flag_snap_type = 1;

fig_num = 555;
St_points = fcn_Path_convertXY2St(referencePath,XY_points, flag_snap_type,fig_num);
assert(length(St_points(:,1))==length(XY_points(:,1)));

%% BASIC example 6 - FLAG 1, use the prior segment
% A 90-degree line segment, a simple query, zero distance in rear segments
XY_points = [ 2 1];
referencePath = [-1 0; 1 0; 1 -1];
flag_snap_type = 1;

fig_num = 666;
St_points = fcn_Path_convertXY2St(referencePath,XY_points, flag_snap_type,fig_num);
expected_solution = [2 1-1i];
assert(abs(sum((St_points - expected_solution).^2,2))<1E-10);

%% BASIC example 7 - FLAG 2
% A 90-degree line segment, a simple query, zero distance in rear segments
XY_points = [ 2 1];
referencePath = [-1 0; 1 0; 1 -1];
flag_snap_type = 2;

fig_num = 777;
St_points = fcn_Path_convertXY2St(referencePath,XY_points, flag_snap_type,fig_num);
expected_solution = [2 1+1i];
assert(abs(sum((St_points - expected_solution).^2,2))<1E-10);

%% BASIC example 8 - FLAG 3
% A 90-degree line segment, a simple query, zero distance in rear segments
XY_points = [ 2 1];
referencePath = [-1 0; 1 0; 1 -1];
flag_snap_type = 3;

fig_num = 888;
St_points = fcn_Path_convertXY2St(referencePath,XY_points, flag_snap_type,fig_num);
expected_solution = [2 2^0.5];
assert(abs(sum((St_points - expected_solution).^2,2))<1E-10);

%% Illustrative example of fcn_Path_convertXY2St
XY_points = [-2 -1; -1 0; -0.5 0.4; 0 0; 0.5 -0.5; 1 -0.4];
referencePath = [-3 -3; -1 -0.5; 0.5 0; 3 3];
flag_snap_type = 3;

St_points_XY = fcn_Path_convertXY2St(referencePath,XY_points, flag_snap_type);
St_points_ref = fcn_Path_convertXY2St(referencePath,referencePath, flag_snap_type);


fig_num = 999;
figure(fig_num);
clf;

subplot(1,2,1);
hold on;
grid on;
axis equal;

plot(XY_points(:,1),XY_points(:,2),'b.-','LineWidth',3,'MarkerSize',20)
plot(referencePath(:,1),referencePath(:,2),'r.-','LineWidth',3,'MarkerSize',20)
title('XY coordinates');

subplot(1,2,2);
hold on;
grid on;
axis equal;
plot(St_points_XY(:,1),St_points_XY(:,2),'b.-','LineWidth',3,'MarkerSize',20)
plot(St_points_ref(:,1),St_points_ref(:,2),'r.-','LineWidth',3,'MarkerSize',20)
title('St coordinates');









