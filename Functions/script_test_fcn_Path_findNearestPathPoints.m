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

%% BASIC example: Small path with single query point
fig_num = 10001;
titleString = sprintf('BASIC example: start at zero');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

query_points = [0.8 0.4];
pathXY = [0 0; 1 0; 2 0; 2 1];

closest_path_point_indicies = ...
    fcn_Path_findNearestPathPoints(query_points, pathXY, fig_num);
fprintf(1,'Figure: %d,\n\t\t Closest point index is: %.0d %.2f\n',...
    fig_num, closest_path_point_indicies(1,1));
title('fcn_Path_findNearestPathPoints: tested with a single point query in 2D','Interpreter','none');

% Check variable types
assert(isnumeric(closest_path_point_indicies));

% Check variable sizes
assert(isequal(size(closest_path_point_indicies),[length(query_points(:,1)) 1]));

% Check variable values
assert(isequal(closest_path_point_indicies,2));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example: Small path with two query points
fig_num = 10002;
titleString = sprintf('BASIC example: start at zero');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

query_points = [0.8 0.4; 1.6 0.8];
pathXY = [0 0; 1 0; 2 0; 2 1];

closest_path_point_indicies = ...
    fcn_Path_findNearestPathPoints(query_points, pathXY,fig_num);
title('fcn_Path_findNearestPathPoints: tested with two point queries in 2D','Interpreter','none');
fprintf(1,'Figure: %d,\n\t\t Closest point index is: %.0d %.2f\n',...
    fig_num, closest_path_point_indicies(1,1));

% Check variable types
assert(isnumeric(closest_path_point_indicies));

% Check variable sizes
assert(isequal(size(closest_path_point_indicies),[length(query_points(:,1)) 1]));

% Check variable values
assert(isequal(closest_path_point_indicies,[2; 4]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example: Small path with many query points
fig_num = 10003;
titleString = sprintf('BASIC example: start at zero');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

query_points = rand(20,2);
query_points(:,1)=query_points(:,1)*4-1;
query_points(:,2)=query_points(:,2)*3-1.5;

pathXY = [0 0; 1 0; 2 0; 2 1];

closest_path_point_indicies = ...
    fcn_Path_findNearestPathPoints(query_points, pathXY,fig_num);
fprintf(1,'Figure: %d,\n\t\t Closest point index is: %.0d %.2f\n',...
    fig_num, closest_path_point_indicies(1,1));
title('fcn_Path_findNearestPathPoints: tested with many point queries in 2D','Interpreter','none');

% Check variable types
assert(isnumeric(closest_path_point_indicies));

% Check variable sizes
assert(isequal(size(closest_path_point_indicies),[length(query_points(:,1)) 1]));

% Check variable values
% assert(isequal(closest_path_point_indicies,[2; 4]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


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


%% BASIC example 3D: simple 3D snapping onto a vertex
fig_num = 20001;
titleString = sprintf('BASIC example: start at zero');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

query_points = [0.8 1.3 2.1];
pathXYZ = [0 0 0; 0.5 0.2 0.4; 0.9 0.9 0.8; 3 0 1];

closest_path_point_indicies = ...
    fcn_Path_findNearestPathPoints(query_points, pathXYZ,fig_num);
title('fcn_Path_findNearestPathPoints: tested with single point queries in 3D','Interpreter','none');

fprintf(1,'Figure: %d, Closest point index is: %.0d \n',...
    fig_num, closest_path_point_indicies(1,1));
view(3);

% Check variable types
assert(isnumeric(closest_path_point_indicies));

% Check variable sizes
assert(isequal(size(closest_path_point_indicies),[length(query_points(:,1)) 1]));

% Check variable values
assert(isequal(closest_path_point_indicies, 3 ));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example 3D: simple 3D snapping onto a vertex
fig_num = 20002;
titleString = sprintf('BASIC example: start at zero');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

query_points = [2 1.3 2.1];
pathXYZ = [0 0 0; 0.5 0.2 0.4; 0.9 0.9 0.8; 3 0 1];

closest_path_point_indicies = ...
    fcn_Path_findNearestPathPoints(query_points, pathXYZ,fig_num);
title('fcn_Path_findNearestPathPoints: tested with single point queries in 3D','Interpreter','none');
xlabel('X'); ylabel('Y'); zlabel('Z');

fprintf(1,'Figure: %d, Closest point index is: %.0d \n',...
    fig_num, closest_path_point_indicies(1,1));
view(3);

% Check variable types
assert(isnumeric(closest_path_point_indicies));

% Check variable sizes
assert(isequal(size(closest_path_point_indicies),[length(query_points(:,1)) 1]));

% Check variable values
assert(isequal(closest_path_point_indicies,3));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example 3D: simple 3D snapping with two query points
fig_num = 20003;
titleString = sprintf('BASIC example: start at zero');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

query_points = [0.8 1.3 2.1; 2 1.3 2.1];
pathXYZ = [0 0 0; 0.5 0.2 0.4; 0.9 0.9 0.8; 3 0 1];

closest_path_point_indicies = ...
    fcn_Path_findNearestPathPoints(query_points, pathXYZ,fig_num);
title('fcn_Path_findNearestPathPoints: tested with single point queries in 3D','Interpreter','none');
xlabel('X'); ylabel('Y'); zlabel('Z');

fprintf(1,'Figure: %d, Closest point index is: %.0d \n',...
    fig_num, closest_path_point_indicies(1,1));
view(3);

% Check variable types
assert(isnumeric(closest_path_point_indicies));

% Check variable sizes
assert(isequal(size(closest_path_point_indicies),[length(query_points(:,1)) 1]));

% Check variable values
assert(isequal(closest_path_point_indicies,[3;3]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example 3: Small path with many query points
fig_num = 20004;
titleString = sprintf('BASIC example: start at zero');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

query_points = rand(20,3);
query_points(:,1)=query_points(:,1)*5-2.5;
query_points(:,2)=query_points(:,2)*3-1.5;
query_points(:,3)=query_points(:,2)*3-1.5;

pathXYZ = [0 0 0; 0.5 0.2 0.4; 0.9 0.9 0.8; 3 0 1];

closest_path_point_indicies = ...
    fcn_Path_findNearestPathPoints(query_points, pathXYZ,fig_num);
title('fcn_Path_findNearestPathPoints: tested with single point queries in 3D','Interpreter','none');
xlabel('X'); ylabel('Y'); zlabel('Z');


fprintf(1,'Figure: %d, Closest point index is: %.0d \n',...
    fig_num, closest_path_point_indicies(1,1));
view(3);

% Check variable types
assert(isnumeric(closest_path_point_indicies));

% Check variable sizes
assert(isequal(size(closest_path_point_indicies),[length(query_points(:,1)) 1]));

% Check variable values
% assert(isequal(closest_path_point_indicies,3));

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

query_points = [0.8 0.4; 1.6 0.8];
pathXY = [0 0; 1 0; 2 0; 2 1];

closest_path_point_indicies = ...
    fcn_Path_findNearestPathPoints(query_points, pathXY,[]);

% Check variable types
assert(isnumeric(closest_path_point_indicies));

% Check variable sizes
assert(isequal(size(closest_path_point_indicies),[length(query_points(:,1)) 1]));

% Check variable values
assert(isequal(closest_path_point_indicies,[2; 4]));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

query_points = [0.8 0.4; 1.6 0.8];
pathXY = [0 0; 1 0; 2 0; 2 1];

closest_path_point_indicies = ...
    fcn_Path_findNearestPathPoints(query_points, pathXY,(-1));

% Check variable types
assert(isnumeric(closest_path_point_indicies));

% Check variable sizes
assert(isequal(size(closest_path_point_indicies),[length(query_points(:,1)) 1]));

% Check variable values
assert(isequal(closest_path_point_indicies,[2; 4]));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

query_points = [0.8 0.4; 1.6 0.8];
pathXY = [0 0; 1 0; 2 0; 2 1];


Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function

    closest_path_point_indicies = ...
        fcn_Path_findNearestPathPoints(query_points, pathXY,([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function

    closest_path_point_indicies = ...
        fcn_Path_findNearestPathPoints(query_points, pathXY,(-1));
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
