% script_test_fcn_Path_convertSt2XY.m
% This is a script to exercise the function: fcn_Path_convertSt2XY.m
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
XY_expected_solution = [0 1];
St_points_input = [1 1];
referencePath = [-1 0; 1 0];
flag_snap_type = 1;

fig_num = 111;
XY_points_calculated = fcn_Path_convertSt2XY(referencePath,St_points_input, flag_snap_type,fig_num);

assert(abs(sum((XY_points_calculated - XY_expected_solution).^2,2))<1E-10);

%% BASIC example 2
% A simple line segment, a simple query, zero distance in rear segments
XY_expected_solution = [0 1];
St_points_input = [1 1]./(2^0.5);
referencePath = [0 0; 2 2];
flag_snap_type = 1;

fig_num = 222;
XY_points_calculated = fcn_Path_convertSt2XY(referencePath,St_points_input, flag_snap_type,fig_num);

assert(abs(sum((XY_points_calculated - XY_expected_solution).^2,2))<1E-10);


%% BASIC example 3
% A simple line segment, a simple query, zero distance in rear segments
XY_expected_solution = [0 -1];
St_points_input = [1 -1];
referencePath = [-1 0; 1 0];
flag_snap_type = 1;

fig_num = 333;
XY_points_calculated = fcn_Path_convertSt2XY(referencePath,St_points_input, flag_snap_type,fig_num);

assert(abs(sum((XY_points_calculated - XY_expected_solution).^2,2))<1E-10);

%% BASIC example 4
% A simple line segment, a complex number in rear segments
XY_expected_solution = [-2 -1];
St_points_input = [-1 -1+1i];
referencePath = [-1 0; 1 0];
flag_snap_type = 1;

fig_num = 444;
XY_points_calculated = fcn_Path_convertSt2XY(referencePath,St_points_input, flag_snap_type,fig_num);

assert(abs(sum((XY_points_calculated - XY_expected_solution).^2,2))<1E-10);



%% BASIC example 6 - FLAG 1, use the prior segment
% A 90-degree line segment, a simple query, zero distance in rear segments
XY_expected_solution = [ 2 1];
St_points_input = [2 1-1i];
referencePath = [-1 0; 1 0; 1 -1];
flag_snap_type = 1;

fig_num = 666;
XY_points_calculated = fcn_Path_convertSt2XY(referencePath,St_points_input, flag_snap_type,fig_num);

assert(abs(sum((XY_points_calculated - XY_expected_solution).^2,2))<1E-10);

%% BASIC example 7 - FLAG 2
% A 90-degree line segment, a simple query, zero distance in rear segments
XY_expected_solution = [ 2 1];
St_points_input = [2 1+1i];
referencePath = [-1 0; 1 0; 1 -1];
flag_snap_type = 2;

fig_num = 777;
XY_points_calculated = fcn_Path_convertSt2XY(referencePath,St_points_input, flag_snap_type,fig_num);

assert(abs(sum((XY_points_calculated - XY_expected_solution).^2,2))<1E-10);

%% BASIC example 8 - FLAG 3
% A 90-degree line segment, a simple query, zero distance in rear segments
XY_expected_solution = [ 2 1];
St_points_input = [2 2^0.5];
referencePath = [-1 0; 1 0; 1 -1];
flag_snap_type = 3;

fig_num = 888;
XY_points_calculated = fcn_Path_convertSt2XY(referencePath,St_points_input, flag_snap_type,fig_num);

assert(abs(sum((XY_points_calculated - XY_expected_solution).^2,2))<1E-10);


%% BASIC example 5 - many points

% A 90-degree line segment with multiple surrounding queries
XY_points = [-2 1; -1 1; 0 1; 1 1; 2 1; 2 0; 2 -1; 2 -2; 1 -2; 0 -2; 0 -1; -1 -1; -2 -1; -2 0];
XY_points = XY_points(6,:);
referencePath = [-1 0; 1 0; 1 -1];
flag_snap_type = 1;

fig_num = 555;
St_points = fcn_Path_convertXY2St(referencePath,XY_points, flag_snap_type,fig_num);
assert(length(St_points(:,1))==length(XY_points(:,1)));


% Now, convert them back
fig_num = 666;
XY_points_calculated = fcn_Path_convertSt2XY(referencePath,St_points, flag_snap_type,fig_num);
figure(555);
plot(XY_points_calculated(:,1),XY_points_calculated(:,2),'bo','MarkerSize',20);

assert(abs(sum(sum((XY_points_calculated - XY_points).^2,2)))<1E-10);
