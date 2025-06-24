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

%% BASIC example
% A simple line segment, a simple query, zero distance in rear segments
fig_num = 10001;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

XY_points = [0 1];
referencePath = [-1 0; 1 0];
flag_snap_type = 1;

St_points = fcn_Path_convertXY2St(referencePath,XY_points, flag_snap_type,fig_num);
expected_solution = [1 1];
assert(abs(sum((St_points - expected_solution).^2,2))<1E-10);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example 
% A simple line segment, a simple query, zero distance in rear segments
fig_num = 10002;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

XY_points = [0 1];
referencePath = [0 0; 2 2];
flag_snap_type = 1;

St_points = fcn_Path_convertXY2St(referencePath,XY_points, flag_snap_type,fig_num);
expected_solution = [1 1]./(2^0.5);
assert(abs(sum((St_points - expected_solution).^2,2))<1E-10);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example
% A simple line segment, a simple query, zero distance in rear segments
fig_num = 10003;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

XY_points = [0 -1];
referencePath = [-1 0; 1 0];
flag_snap_type = 1;

St_points = fcn_Path_convertXY2St(referencePath,XY_points, flag_snap_type,fig_num);
expected_solution = [1 -1];
assert(abs(sum((St_points - expected_solution).^2,2))<1E-10);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example
% A simple line segment, a complex number in rear segments
fig_num = 10004;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

XY_points = [-2 -1];
referencePath = [-1 0; 1 0];
flag_snap_type = 1;

St_points = fcn_Path_convertXY2St(referencePath,XY_points, flag_snap_type,fig_num);
expected_solution = [-1 -1+1i];
assert(abs(sum((St_points - expected_solution).^2,2))<1E-10);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example 
% A simple line egment, a complex number in front segments
fig_num = 10005;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

XY_points = [2 1];
referencePath = [-1 0; 1 0];
flag_snap_type = 1;

St_points = fcn_Path_convertXY2St(referencePath,XY_points, flag_snap_type,fig_num);
expected_solution = [3 1-1i];
assert(abs(sum((St_points - expected_solution).^2,2))<1E-10);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example - many points
% A 90-degree line segment with multiple surrounding queries
fig_num = 10006;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

XY_points = [-2 1; -1 1; 0 1; 1 1; 2 1; 2 0; 2 -1; 2 -2; 1 -2; 0 -2; 0 -1; -1 -1; -2 -1; -2 0];
referencePath = [-1 0; 1 0; 1 -1];
flag_snap_type = 1;

St_points = fcn_Path_convertXY2St(referencePath,XY_points, flag_snap_type,fig_num);
assert(length(St_points(:,1))==length(XY_points(:,1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example 6 - FLAG 1, use the prior segment
% A 90-degree line segment, a simple query, zero distance in rear segments
fig_num = 10007;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

XY_points = [ 2 1];
referencePath = [-1 0; 1 0; 1 -1];
flag_snap_type = 1;

St_points = fcn_Path_convertXY2St(referencePath,XY_points, flag_snap_type,fig_num);
expected_solution = [2 1-1i];
assert(abs(sum((St_points - expected_solution).^2,2))<1E-10);

%% BASIC example - FLAG 2
% A 90-degree line segment, a simple query, zero distance in rear segments
fig_num = 10008;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

XY_points = [ 2 1];
referencePath = [-1 0; 1 0; 1 -1];
flag_snap_type = 2;

St_points = fcn_Path_convertXY2St(referencePath,XY_points, flag_snap_type,fig_num);
expected_solution = [2 1+1i];
assert(abs(sum((St_points - expected_solution).^2,2))<1E-10);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example - FLAG 3
% A 90-degree line segment, a simple query, zero distance in rear segments
fig_num = 10009;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

XY_points = [ 2 1];
referencePath = [-1 0; 1 0; 1 -1];
flag_snap_type = 3;

St_points = fcn_Path_convertXY2St(referencePath,XY_points, flag_snap_type,fig_num);
expected_solution = [2 2^0.5];
assert(abs(sum((St_points - expected_solution).^2,2))<1E-10);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Illustrative example of fcn_Path_convertXY2St
fig_num = 10010;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;


XY_points = [-2 -1; -1 0; -0.5 0.4; 0 0; 0.5 -0.5; 1 -0.4];
referencePath = [-3 -3; -1 -0.5; 0.5 0; 3 3];
flag_snap_type = 3;

St_points_XY = fcn_Path_convertXY2St(referencePath,XY_points, flag_snap_type);
St_points_ref = fcn_Path_convertXY2St(referencePath,referencePath, flag_snap_type);

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

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


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

XY_points = [0 1];
referencePath = [0 0; 2 2];
flag_snap_type = 1;

St_points = fcn_Path_convertXY2St(referencePath,XY_points, flag_snap_type,[]);
expected_solution = [1 1]./(2^0.5);
assert(abs(sum((St_points - expected_solution).^2,2))<1E-10);

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

XY_points = [0 1];
referencePath = [0 0; 2 2];
flag_snap_type = 1;

St_points = fcn_Path_convertXY2St(referencePath,XY_points, flag_snap_type,(-1));
expected_solution = [1 1]./(2^0.5);
assert(abs(sum((St_points - expected_solution).^2,2))<1E-10);

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);


XY_points = [0 1];
referencePath = [0 0; 2 2];
flag_snap_type = 1;

Niterations = 20;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    St_points = fcn_Path_convertXY2St(referencePath,XY_points, flag_snap_type,[]);
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    St_points = fcn_Path_convertXY2St(referencePath,XY_points, flag_snap_type,-1);

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

