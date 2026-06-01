% script_test_fcn_Path_calculateTransitionCurves
% Tests fcn_Path_calculateTransitionCurves
       
% Revision history:
%  2023_07_17
%  -- first write of the code by V. Wagh (vbw5054@psu.edu)
%  2023_07_31 by S. Brennan (sbrennan@psu.edu)
%  -- added assertions
%  -- slight reformatting, including renaming segments to paths
% 2023_08_09 by S. Brennan
% - Converted this function over to Path library from LoadWZ librar

%% clear the workspace
close all


% Clear any old variables
clear test_traversals

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
figNum = 10001;
fprintf(1,'Figure %.0f: basic demo 1\n',figNum);
figure(figNum); clf;

path_1 = [0 0; 4 4];
path_2 = [4 1; 1 4];
radius_of_curve = 1;
plot_color = [];
plot_line_width = [];
plot_text = '';
number_of_points_on_curve = 50;

[ TransitionCurves, ...
    closest_path_point1, closest_path_point2,...
    angle_point1_radians,angle_point2_radians] = ...
    fcn_Path_calculateTransitionCurves(...
    path_1, ...
    path_2,...
    radius_of_curve, ...
    number_of_points_on_curve,...
    plot_color,...
    plot_line_width,...
    plot_text,...
    figNum);

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(isscalar(closest_path_point1(:,1)));
assert(isscalar(closest_path_point2(:,1)));

% Check the angles
assert(isequal(round(angle_point1_radians,4),round(-pi/4,4)));
assert(abs(angle_point2_radians - pi/4)<eps*1000);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% Basic example 2: simple line segment with straight lines and negative radius
figNum = 10002;
fprintf(1,'Figure %.0f: basic demo 2\n',figNum);
figure(figNum); clf;


path_1 = [0 0; 4 4];
path_2 = flipud([4 1; 1 4]);  % <--- note direction change
radius_of_curve = 1;
plot_color = [];
plot_line_width = [];
plot_text = '';
number_of_points_on_curve = 50;

[ TransitionCurves, ...
    closest_path_point1, closest_path_point2,...
    angle_point1_radians,angle_point2_radians] = ...
    fcn_Path_calculateTransitionCurves(...
    path_1, ...
    path_2,...
    radius_of_curve, ...
    number_of_points_on_curve,...
    plot_color,...
    plot_line_width,...
    plot_text,...
    figNum);

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(isscalar(closest_path_point1(:,1)));
assert(isscalar(closest_path_point2(:,1)));
assert(isscalar(angle_point1_radians(:,1)));
assert(isscalar(angle_point2_radians(:,1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% Basic example 3: simple line segment with line changing direction 
figNum = 10003;
fprintf(1,'Figure %.0f: basic demo 1\n',figNum);
figure(figNum); clf;

path_1 = flipud([0 0; 4 4]); % <--- note direction change
path_2 = flipud([4 1; 1 4]);  % <--- note direction change
radius_of_curve = 1;
plot_color = [];
plot_line_width = [];
plot_text = '';
number_of_points_on_curve = 50;

[ TransitionCurves, ...
    closest_path_point1, closest_path_point2,...
    angle_point1_radians,angle_point2_radians] = ...
    fcn_Path_calculateTransitionCurves(...
    path_1, ...
    path_2,...
    radius_of_curve, ...
    number_of_points_on_curve,...
    plot_color,...
    plot_line_width,...
    plot_text,...
    figNum);

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(isscalar(closest_path_point1(:,1)));
assert(isscalar(closest_path_point2(:,1)));
assert(isscalar(angle_point1_radians(:,1)));
assert(isscalar(angle_point2_radians(:,1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% Basic example 4: simple line segment with line changing direction 
figNum = 10004;
fprintf(1,'Figure %.0f: basic demo 1\n',figNum);
figure(figNum); clf;

path_1 = flipud([0 0; 4 4]); % <--- note direction change
path_2 = [4 1; 1 4]; 
radius_of_curve = 1;
plot_color = [];
plot_line_width = [];
plot_text = '';
number_of_points_on_curve = 50;

[ TransitionCurves, ...
    closest_path_point1, closest_path_point2,...
    angle_point1_radians,angle_point2_radians] = ...
    fcn_Path_calculateTransitionCurves(...
    path_1, ...
    path_2,...
    radius_of_curve, ...
    number_of_points_on_curve,...
    plot_color,...
    plot_line_width,...
    plot_text,...
    figNum);

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(isscalar(closest_path_point1(:,1)));
assert(isscalar(closest_path_point2(:,1)));
assert(isscalar(angle_point1_radians(:,1)));
assert(isscalar(angle_point2_radians(:,1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% Basic example 5a: simple curve, outside fit
figNum = 100051;
fprintf(1,'Figure %.0f: basic demo 1\n',figNum);
figure(figNum); clf;


th = linspace( pi/2, -pi/2, 100);
Radius = 10;  %or whatever radius you want
semi_circle_x = Radius*cos(th) + 5;
semi_circle_y = Radius*sin(th) + 4;
path_1 = [semi_circle_x', semi_circle_y'];
path_2 = [5 0; 30 0];
radius_of_curve = 6;
number_of_points_on_curve = 50;
plot_color = [];
plot_line_width = 3;
plot_text = 'Transition_Curve';

[ TransitionCurves, ...
    closest_path_point1, closest_path_point2,...
    angle_point1_radians,angle_point2_radians] = ...
    fcn_Path_calculateTransitionCurves(...
    path_1, ...
    path_2,...
    radius_of_curve, ...
    number_of_points_on_curve,...
    plot_color,...
    plot_line_width,...
    plot_text,...
    figNum);

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(isscalar(closest_path_point1(:,1)));
assert(isscalar(closest_path_point2(:,1)));
assert(isscalar(angle_point1_radians(:,1)));
assert(isscalar(angle_point2_radians(:,1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% Basic example 5b: simple curve, inside fit
figNum = 100052;
fprintf(1,'Figure %.0f: basic demo 1\n',figNum);
figure(figNum); clf;


figure(figNum);
clf;
hold on;
grid on;


th = linspace( pi/2, -pi/2, 100);
Radius = 10;  %or whatever radius you want
semi_circle_x = Radius*cos(th) + 5;
semi_circle_y = Radius*sin(th) + 4;
path_1 = [semi_circle_x', semi_circle_y'];
path_2 = [5 0; 30 0];
radius_of_curve = 2;
number_of_points_on_curve = 50;
plot_color = [];
plot_line_width = 3;
plot_text = 'Transition_Curve';

[ TransitionCurves, ...
    closest_path_point1, closest_path_point2,...
    angle_point1_radians,angle_point2_radians] = ...
    fcn_Path_calculateTransitionCurves(...
    path_2, ...
    path_1,...
    radius_of_curve, ...
    number_of_points_on_curve,...
    plot_color,...
    plot_line_width,...
    plot_text,...
    figNum);

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(isscalar(closest_path_point1(:,1)));
assert(isscalar(closest_path_point2(:,1)));
assert(isscalar(angle_point1_radians(:,1)));
assert(isscalar(angle_point2_radians(:,1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% Basic example 5c: simple curve, inside fit, other direction
figNum = 100053;
fprintf(1,'Figure %.0f: basic demo 1\n',figNum);
figure(figNum); clf;


figure(figNum);
clf;
hold on;
grid on;


th = linspace( pi/2, -pi/2, 100);
Radius = 10;  %or whatever radius you want
semi_circle_x = Radius*cos(th) + 5;
semi_circle_y = Radius*sin(th) + 4;
path_1 = [semi_circle_x', semi_circle_y'];
path_2 = [5 0; 30 0];
radius_of_curve = 2;
number_of_points_on_curve = 50;
plot_color = [];
plot_line_width = 3;
plot_text = 'Transition_Curve';

[ TransitionCurves, ...
    closest_path_point1, closest_path_point2,...
    angle_point1_radians,angle_point2_radians] = ...
    fcn_Path_calculateTransitionCurves(...
    path_2, ...
    flipud(path_1),...
    radius_of_curve, ...
    number_of_points_on_curve,...
    plot_color,...
    plot_line_width,...
    plot_text,...
    figNum);

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(isscalar(closest_path_point1(:,1)));
assert(isscalar(closest_path_point2(:,1)));
assert(isscalar(angle_point1_radians(:,1)));
assert(isscalar(angle_point2_radians(:,1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% Basic example 5d: simple curve, inside fit, other direction, coarse
figNum = 100054;
fprintf(1,'Figure %.0f: basic demo 1\n',figNum);
figure(figNum); clf;

th = linspace( pi/2, -pi/2, 5);
Radius = 10;  %or whatever radius you want
semi_circle_x = Radius*cos(th) + 5;
semi_circle_y = Radius*sin(th) + 4;
path_1 = [semi_circle_x', semi_circle_y'];
path_2 = [5 0; 30 0];
radius_of_curve = 5;
number_of_points_on_curve = 50;
plot_color = [];
plot_line_width = 3;
plot_text = 'Transition_Curve';

[ TransitionCurves, ...
    closest_path_point1, closest_path_point2,...
    angle_point1_radians,angle_point2_radians] = ...
    fcn_Path_calculateTransitionCurves(...
    path_2, ...
    flipud(path_1),...
    radius_of_curve, ...
    number_of_points_on_curve,...
    plot_color,...
    plot_line_width,...
    plot_text,...
    figNum);

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(isscalar(closest_path_point1(:,1)));
assert(isscalar(closest_path_point2(:,1)));
assert(isscalar(angle_point1_radians(:,1)));
assert(isscalar(angle_point2_radians(:,1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%%  Basic example 6a: multiple intersections between the line segments 
figNum = 100061;
fprintf(1,'Figure %.0f: basic demo 1\n',figNum);
figure(figNum); clf;


path_1 = [2 2; 5 8];
path_2 = [3 2; 3 6; 5 6];
radius_of_curve = 2;

plot_color = [];
plot_line_width = [];
plot_text = [];

number_of_points_on_curve = 50;
[ TransitionCurves, ...
    closest_path_point1, closest_path_point2,...
    angle_point1_radians,angle_point2_radians] = ...
    fcn_Path_calculateTransitionCurves(...
    path_1, ...
    path_2,...
    radius_of_curve, ...
    number_of_points_on_curve,...
    plot_color,...
    plot_line_width,...
    plot_text,...
    figNum);

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(isscalar(closest_path_point1(:,1)));
assert(isscalar(closest_path_point2(:,1)));
assert(isscalar(angle_point1_radians(:,1)));
assert(isscalar(angle_point2_radians(:,1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%%  Basic example 6b: multiple intersections between the line segments 
figNum = 100062;
fprintf(1,'Figure %.0f: basic demo 1\n',figNum);
figure(figNum); clf;


path_1 = [2 2; 5 8];
path_2 = [5 2; 2 6; 5 6];
radius_of_curve = 2;

plot_color = [];
plot_line_width = [];
plot_text = [];

number_of_points_on_curve = 50;
[ TransitionCurves, ...
    closest_path_point1, closest_path_point2,...
    angle_point1_radians,angle_point2_radians] = ...
    fcn_Path_calculateTransitionCurves(...
    path_1, ...
    path_2,...
    radius_of_curve, ...
    number_of_points_on_curve,...
    plot_color,...
    plot_line_width,...
    plot_text,...
    figNum);

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(isscalar(closest_path_point1(:,1)));
assert(isscalar(closest_path_point2(:,1)));
assert(isscalar(angle_point1_radians(:,1)));
assert(isscalar(angle_point2_radians(:,1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%%  Basic example 7: multiple intersections between the line segments 
figNum = 10007;
fprintf(1,'Figure %.0f: basic demo 1\n',figNum);
figure(figNum); clf;


path_1 = [2 2; 5 8];
path_2 = [5 2; 2 6; 5 6];
radius_of_curve = 2;
number_of_points_on_curve = 50;

plot_color = [];
plot_line_width = [];
plot_text = '';

[ TransitionCurves, ...
    closest_path_point1, closest_path_point2,...
    angle_point1_radians,angle_point2_radians] = ...
    fcn_Path_calculateTransitionCurves(...
    path_1, ...
    path_2,...
    radius_of_curve, ...
    number_of_points_on_curve,...
    plot_color,...
    plot_line_width,...
    plot_text,...
    figNum);

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(isscalar(closest_path_point1(:,1)));
assert(isscalar(closest_path_point2(:,1)));
assert(isscalar(angle_point1_radians(:,1)));
assert(isscalar(angle_point2_radians(:,1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% HARD example 8: line segments have infinite intersections for a period
% figNum = 10008;
% fprintf(1,'Figure %.0f: basic demo 1\n',figNum);
% figure(figNum); clf;
% 
% 
% path_1 = [0 0; 4 4; 6 6; 10 5];
% path_2 = [4 1; 4 3; 5 5; 7 7; 5 10];
% radius_of_curve = 1;
% 
% plot_color = [];
% plot_line_width = [];
% plot_text = '';
% 
% number_of_points_on_curve = 50;
% 
% [ TransitionCurves, ...
%     closest_path_point1, closest_path_point2,...
%     angle_point1_radians,angle_point2_radians] = ...
%     fcn_Path_calculateTransitionCurves(...
%     path_1, ...
%     path_2,...
%     radius_of_curve, ...
%     number_of_points_on_curve,...
%     plot_color,...
%     plot_line_width,...
%     plot_text,...
%     figNum);
% 
% % Check the size of the vectors
% assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
% assert(length(closest_path_point1(:,1))==1);
% assert(length(closest_path_point2(:,1))==1);
% assert(length(angle_point1_radians(:,1))==1);
% assert(length(angle_point2_radians(:,1))==1);
% 
% % Make sure plot opened up
% assert(isequal(get(gcf,'Number'),figNum));


%% Example for which code breaks
% should get error: Lines are collinear, curve is not possible between them

% figNum = 10009;
% fprintf(1,'Figure %.0f: basic demo 1\n',figNum);
% figure(figNum); clf;
% 
% 
% path_1 = [4 9; 4 5; 2 2];
% path_2 = [4 9; 4 2];
% radius_of_curve = 2;
% number_of_points_on_curve = 50;
% 
% plot_color = [];
% plot_line_width = [];
% plot_text = '';
% 
% [ TransitionCurves, ...
%     closest_path_point1, closest_path_point2,...
%     angle_point1_radians,angle_point2_radians] = ...
%     fcn_Path_calculateTransitionCurves(...
%     path_1, ...
%     path_2,...
%     radius_of_curve, ...
%     number_of_points_on_curve,...
%     plot_color,...
%     plot_line_width,...
%     plot_text,...
%     figNum);
% 
% % Check the size of the vectors
% assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
% assert(length(closest_path_point1(:,1))==1);
% assert(length(closest_path_point2(:,1))==1);
% assert(length(angle_point1_radians(:,1))==1);
% assert(length(angle_point2_radians(:,1))==1);
% 
% % Make sure plot opened up
% assert(isequal(get(gcf,'Number'),figNum));


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
figNum = 80001;
fprintf(1,'Figure: %.0f: Demo of fast mode, empty figNum\n',figNum);
figure(figNum); close(figNum);

path_1 = [0 0; 4 4];
path_2 = [4 1; 1 4];
radius_of_curve = 1;
plot_color = [];
plot_line_width = [];
plot_text = '';
number_of_points_on_curve = 50;

[ TransitionCurves, ...
    closest_path_point1, closest_path_point2,...
    angle_point1_radians,angle_point2_radians] = ...
    fcn_Path_calculateTransitionCurves(...
    path_1, ...
    path_2,...
    radius_of_curve, ...
    number_of_points_on_curve,...
    plot_color,...
    plot_line_width,...
    plot_text,...
    []);

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(isscalar(closest_path_point1(:,1)));
assert(isscalar(closest_path_point2(:,1)));

% Check the angles
assert(isequal(round(angle_point1_radians,4),round(-pi/4,4)));
assert(abs(angle_point2_radians - pi/4)<eps*1000);

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Basic fast mode - NO FIGURE, FAST MODE
figNum = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, figNum=-1\n',figNum);
figure(figNum); close(figNum);

path_1 = [0 0; 4 4];
path_2 = [4 1; 1 4];
radius_of_curve = 1;
plot_color = [];
plot_line_width = [];
plot_text = '';
number_of_points_on_curve = 50;

[ TransitionCurves, ...
    closest_path_point1, closest_path_point2,...
    angle_point1_radians,angle_point2_radians] = ...
    fcn_Path_calculateTransitionCurves(...
    path_1, ...
    path_2,...
    radius_of_curve, ...
    number_of_points_on_curve,...
    plot_color,...
    plot_line_width,...
    plot_text,...
    -1);

% Check the size of the vectors
assert(length(TransitionCurves(:,1))==number_of_points_on_curve);
assert(isscalar(closest_path_point1(:,1)));
assert(isscalar(closest_path_point2(:,1)));

% Check the angles
assert(isequal(round(angle_point1_radians,4),round(-pi/4,4)));
assert(abs(angle_point2_radians - pi/4)<eps*1000);
% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
figNum = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',figNum);
figure(figNum);
close(figNum);

path_1 = [0 0; 4 4];
path_2 = [4 1; 1 4];
radius_of_curve = 1;
plot_color = [];
plot_line_width = [];
plot_text = '';
number_of_points_on_curve = 50;

Niterations = 10;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [ TransitionCurves, ...
        closest_path_point1, closest_path_point2,...
        angle_point1_radians,angle_point2_radians] = ...
        fcn_Path_calculateTransitionCurves(...
        path_1, ...
        path_2,...
        radius_of_curve, ...
        number_of_points_on_curve,...
        plot_color,...
        plot_line_width,...
        plot_text,...
        []);

end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [ TransitionCurves, ...
        closest_path_point1, closest_path_point2,...
        angle_point1_radians,angle_point2_radians] = ...
        fcn_Path_calculateTransitionCurves(...
        path_1, ...
        path_2,...
        radius_of_curve, ...
        number_of_points_on_curve,...
        plot_color,...
        plot_line_width,...
        plot_text,...
        -1);

end
fast_method = toc;

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));

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
assert(~any(figHandles==figNum));


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

%% FAIL CASES
if 1==0
    %% Case where paths completely overlap
    figNum = 4500;
    path_1 = [0 0; 4 4];
    path_2 = path_1;
    radius_of_curve = 1;
    plot_color = [];
    plot_line_width = [];
    plot_text = '';
    number_of_points_on_curve = 50;
    [ TransitionCurves, ...
        closest_path_point1, closest_path_point2,...
        angle_point1_radians,angle_point2_radians] = ...
        fcn_Path_calculateTransitionCurves(...
        path_1, ...
        path_2,...
        radius_of_curve, ...
        number_of_points_on_curve,...
        plot_color,...
        plot_line_width,...
        plot_text,...
        figNum);
end


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


