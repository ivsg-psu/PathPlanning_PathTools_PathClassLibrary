% script_test_fcn_Path_findClosestPointsToTraversal.m
% This is a script to exercise the function: fcn_Path_findClosestPointsToTraversal.m
% This function was written on 2020_11_15 by S. Brennan
%     Modified on 2020_11_15 to prep for Path class
% Questions or comments? sbrennan@psu.edu

% Here is the function format:
% function [closestXs,closestYs,closestZs,closestYaws] = ...
% fcn_Path_findClosestPointsToTraversal(...
% reference_traversal,data,flag_yaw,flag_3D)

% Revision history
% 2021_01_09
% -- added more comments


close all;

%% BASIC example - parallel lines, query is in middle area
fig_num = 10001;
titleString = sprintf('BASIC example - parallel lines, query is in middle area');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Create a dummy central path and convert it to a traversal
central_path = [0 0; 2 0];  
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [0 4; 2 4];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);
clear data
data.traversal{1} = nearby_traversal;
flag_yaw = 0;
flag_3D = 0;

% Calculate the closest point and distance on the nearby path
[closestXs,closestYs,closestZs,closestYaws] = ...
    fcn_Path_findClosestPointsToTraversal(central_traversal,data,flag_yaw,flag_3D,fig_num); %#ok<*ASGLU>

title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestZs));
assert(isnumeric(closestYaws));

% Check variable sizes
NreferenceTraversal = length(central_traversal.X(:,1));
assert(isequal(size(closestXs),[NreferenceTraversal 1]));
assert(isequal(size(closestYs),[NreferenceTraversal 1]));
assert(isequal(size(closestZs),[NreferenceTraversal 1]));
assert(isequal(size(closestYaws),[NreferenceTraversal 1]));

% Check variable values
assert(isequal(closestXs,[0; 2]));
assert(isequal(closestYs,[4; 4]));
assert(isequal(closestZs,[0; 0]));
assert(isequal(closestYaws,[0; 0]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example - angled line segment adjacent to endpoint query
fig_num = 10002;
titleString = sprintf('BASIC example - angled line segment adjacent to endpoint query');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Create a dummy central path and convert it to a traversal
nearby_path = [0 4; 2 7];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);
clear data
data.traversal{1} = nearby_traversal;
flag_yaw = 0;
flag_3D = 0;

% Calculate the closest point and distance on the nearby path
[closestXs,closestYs,closestZs,closestYaws] = ...
    fcn_Path_findClosestPointsToTraversal(central_traversal,data,flag_yaw,flag_3D,fig_num);

title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestZs));
assert(isnumeric(closestYaws));

% Check variable sizes
NreferenceTraversal = length(central_traversal.X(:,1));
assert(isequal(size(closestXs),[NreferenceTraversal 1]));
assert(isequal(size(closestYs),[NreferenceTraversal 1]));
assert(isequal(size(closestZs),[NreferenceTraversal 1]));
assert(isequal(size(closestYaws),[NreferenceTraversal 1]));

% Check variable values
assert(isequal(closestXs,[0; 0]));
assert(isequal(closestYs,[4; 4]));
assert(isequal(closestZs,[0; 0]));
assert(isequal(closestYaws,[0; 0]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example - angled line segment adjacent to endpoint query but near-miss
fig_num = 10003;
titleString = sprintf('BASIC example - angled line segment adjacent to endpoint query but near-miss');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Create a dummy central path and convert it to a traversal
central_path = [0 0; 10 0];
central_traversal = ...
    fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [0 4; 10 7];
nearby_traversal = ...
    fcn_Path_convertPathToTraversalStructure(nearby_path);
clear data
data.traversal{1} = nearby_traversal;
flag_yaw = 0;
flag_3D = 0;

% Calculate the closest point and distance on the nearby path
[closestXs,closestYs,closestZs,closestYaws] = ...
    fcn_Path_findClosestPointsToTraversal(central_traversal,data,flag_yaw,flag_3D,fig_num);

title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestZs));
assert(isnumeric(closestYaws));

% Check variable sizes
NreferenceTraversal = length(central_traversal.X(:,1));
assert(isequal(size(closestXs),[NreferenceTraversal 1]));
assert(isequal(size(closestYs),[NreferenceTraversal 1]));
assert(isequal(size(closestZs),[NreferenceTraversal 1]));
assert(isequal(size(closestYaws),[NreferenceTraversal 1]));

% Check variable values
% assert(isequal(closestXs,[0; 2]));
% assert(isequal(closestYs,[4; 4]));
% assert(isequal(closestZs,[0; 0]));
% assert(isequal(closestYaws,[0; 0]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));
%% BASIC example - angled line segment adjacent to startpoint query
fig_num = 10004;
titleString = sprintf('BASIC example - angled line segment adjacent to startpoint query');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

central_path = [0 0; 10 0];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [-1 4; 12 7];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);
clear data
data.traversal{1} = nearby_traversal;
flag_yaw = 0;
flag_3D = 0;

% Calculate the closest point and distance on the nearby path
[closestXs,closestYs,closestZs,closestYaws] = ...
    fcn_Path_findClosestPointsToTraversal(central_traversal,data,flag_yaw,flag_3D,fig_num);

title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestZs));
assert(isnumeric(closestYaws));

% Check variable sizes
NreferenceTraversal = length(central_traversal.X(:,1));
assert(isequal(size(closestXs),[NreferenceTraversal 1]));
assert(isequal(size(closestYs),[NreferenceTraversal 1]));
assert(isequal(size(closestZs),[NreferenceTraversal 1]));
assert(isequal(size(closestYaws),[NreferenceTraversal 1]));

% Check variable values
% assert(isequal(closestXs,[0; 2]));
% assert(isequal(closestYs,[4; 4]));
% assert(isequal(closestZs,[0; 0]));
% assert(isequal(closestYaws,[0; 0]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));
%% BASIC example - parallel line segment adjacent to startpoint query but near-miss
fig_num = 10005;
titleString = sprintf('BASIC example - parallel line segment adjacent to startpoint query but near-miss');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Create a dummy central path and convert it to a traversal
central_path = [0 0; 10 0];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [0 4; 10 4];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);
clear data
data.traversal{1} = nearby_traversal;
flag_yaw = 0;
flag_3D = 0;

% Calculate the closest point and distance on the nearby path
[closestXs,closestYs,closestZs,closestYaws] = ...
    fcn_Path_findClosestPointsToTraversal(central_traversal,data,flag_yaw,flag_3D,fig_num);

title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestZs));
assert(isnumeric(closestYaws));

% Check variable sizes
NreferenceTraversal = length(central_traversal.X(:,1));
assert(isequal(size(closestXs),[NreferenceTraversal 1]));
assert(isequal(size(closestYs),[NreferenceTraversal 1]));
assert(isequal(size(closestZs),[NreferenceTraversal 1]));
assert(isequal(size(closestYaws),[NreferenceTraversal 1]));

% Check variable values
assert(isequal(closestXs,[0; 10]));
assert(isequal(closestYs,[4; 4]));
assert(isequal(closestZs,[0; 0]));
assert(isequal(closestYaws,[0; 0]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% AVERAGING examples - default setting
fig_num = 20001;
titleString = sprintf('AVERAGING examples - default setting');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Set up data
close all
central_path = [0 0; 1 1; 2 0];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [-1 0.5; 0.5 2; 1.5 2; 3 0.5];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);
clear data
data.traversal{1} = nearby_traversal;
flag_yaw = 0;
flag_3D = 0;

% Calculate the closest point and distance on the nearby path
[closestXs,closestYs,closestZs,closestYaws] = ...
    fcn_Path_findClosestPointsToTraversal(central_traversal,data,flag_yaw,flag_3D,fig_num);

title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestZs));
assert(isnumeric(closestYaws));

% Check variable sizes
NreferenceTraversal = length(central_traversal.X(:,1));
assert(isequal(size(closestXs),[NreferenceTraversal 1]));
assert(isequal(size(closestYs),[NreferenceTraversal 1]));
assert(isequal(size(closestZs),[NreferenceTraversal 1]));
assert(isequal(size(closestYaws),[NreferenceTraversal 1]));

% Check variable values
% assert(isequal(closestXs,[0; 2]));
% assert(isequal(closestYs,[4; 4]));
% assert(isequal(closestZs,[0; 0]));
% assert(isequal(closestYaws,[0; 0]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));
%%

close all;

%% NEGATIVE examples - default setting
fig_num = 30001;
titleString = sprintf('NEGATIVE examples - default setting');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Prep the example
central_path = [-2 1; 1 4; 3 2];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [-1 0.5; 0.5 2; 1.5 2; 3 0.5];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);
clear data
data.traversal{1} = nearby_traversal;
flag_yaw = 0;
flag_3D = 0;

% Calculate the closest point and distance on the nearby path
[closestXs,closestYs,closestZs,closestYaws] = ...
    fcn_Path_findClosestPointsToTraversal(central_traversal,data,flag_yaw,flag_3D,fig_num);

title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestZs));
assert(isnumeric(closestYaws));

% Check variable sizes
NreferenceTraversal = length(central_traversal.X(:,1));
assert(isequal(size(closestXs),[NreferenceTraversal 1]));
assert(isequal(size(closestYs),[NreferenceTraversal 1]));
assert(isequal(size(closestZs),[NreferenceTraversal 1]));
assert(isequal(size(closestYaws),[NreferenceTraversal 1]));

% Check variable values
% assert(isequal(closestXs,[0; 2]));
% assert(isequal(closestYs,[4; 4]));
% assert(isequal(closestZs,[0; 0]));
% assert(isequal(closestYaws,[0; 0]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));
%%

close all;


%% MULTICROSS examples - default setting
fig_num = 40001;
titleString = sprintf('MULTICROSS examples - default setting');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Setup
central_path = [-2 1; 1 4; 3 2; 5 2; 6 3; 7 2];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [-1 0.5; 0.5 2; 1.5 2; 3 0.5; 4 3; 5 4; 6 3; 7 1];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);
clear data
data.traversal{1} = nearby_traversal;
flag_yaw = 0;
flag_3D = 0;

% Calculate the closest point and distance on the nearby path
[closestXs,closestYs,closestZs,closestYaws] = ...
    fcn_Path_findClosestPointsToTraversal(central_traversal,data,flag_yaw,flag_3D,fig_num);

title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestZs));
assert(isnumeric(closestYaws));

% Check variable sizes
NreferenceTraversal = length(central_traversal.X(:,1));
assert(isequal(size(closestXs),[NreferenceTraversal 1]));
assert(isequal(size(closestYs),[NreferenceTraversal 1]));
assert(isequal(size(closestZs),[NreferenceTraversal 1]));
assert(isequal(size(closestYaws),[NreferenceTraversal 1]));

% Check variable values
% assert(isequal(closestXs,[0; 2]));
% assert(isequal(closestYs,[4; 4]));
% assert(isequal(closestZs,[0; 0]));
% assert(isequal(closestYaws,[0; 0]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));
%%

close all;

%% REAL path example
fig_num = 40001;
titleString = sprintf('REAL path example');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

flag_yaw = 0;
flag_3D = 0;

% Fill in some dummy data
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
Ntraversals = 3;
clear data
for i_Path = 1:Ntraversals
    traversal = fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    data.traversal{i_Path} = traversal;
end

% Calculate the closest point and distance on the nearby path
central_traversal = data.traversal{1};
[closestXs,closestYs,closestZs,closestYaws] = ...
    fcn_Path_findClosestPointsToTraversal(central_traversal,data,flag_yaw,flag_3D,fig_num);

title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestZs));
assert(isnumeric(closestYaws));

% Check variable sizes
NreferenceTraversal = length(central_traversal.X(:,1));
assert(isequal(size(closestXs),[NreferenceTraversal Ntraversals]));
assert(isequal(size(closestYs),[NreferenceTraversal Ntraversals]));
assert(isequal(size(closestZs),[NreferenceTraversal Ntraversals]));
assert(isequal(size(closestYaws),[NreferenceTraversal Ntraversals]));

% Check variable values
% assert(isequal(closestXs,[0; 2]));
% assert(isequal(closestYs,[4; 4]));
% assert(isequal(closestZs,[0; 0]));
% assert(isequal(closestYaws,[0; 0]));

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


% Setup
central_path = [-2 1; 1 4; 3 2; 5 2; 6 3; 7 2];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [-1 0.5; 0.5 2; 1.5 2; 3 0.5; 4 3; 5 4; 6 3; 7 1];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);
clear data
data.traversal{1} = nearby_traversal;
flag_yaw = 0;
flag_3D = 0;

% Calculate the closest point and distance on the nearby path
[closestXs,closestYs,closestZs,closestYaws] = ...
    fcn_Path_findClosestPointsToTraversal(central_traversal,data,flag_yaw,flag_3D,([]));

% Check variable types
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestZs));
assert(isnumeric(closestYaws));

% Check variable sizes
NreferenceTraversal = length(central_traversal.X(:,1));
assert(isequal(size(closestXs),[NreferenceTraversal 1]));
assert(isequal(size(closestYs),[NreferenceTraversal 1]));
assert(isequal(size(closestZs),[NreferenceTraversal 1]));
assert(isequal(size(closestYaws),[NreferenceTraversal 1]));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);


% Setup
central_path = [-2 1; 1 4; 3 2; 5 2; 6 3; 7 2];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [-1 0.5; 0.5 2; 1.5 2; 3 0.5; 4 3; 5 4; 6 3; 7 1];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);
clear data
data.traversal{1} = nearby_traversal;
flag_yaw = 0;
flag_3D = 0;

% Calculate the closest point and distance on the nearby path
[closestXs,closestYs,closestZs,closestYaws] = ...
    fcn_Path_findClosestPointsToTraversal(central_traversal,data,flag_yaw,flag_3D,(-1));

% Check variable types
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestZs));
assert(isnumeric(closestYaws));

% Check variable sizes
NreferenceTraversal = length(central_traversal.X(:,1));
assert(isequal(size(closestXs),[NreferenceTraversal 1]));
assert(isequal(size(closestYs),[NreferenceTraversal 1]));
assert(isequal(size(closestZs),[NreferenceTraversal 1]));
assert(isequal(size(closestYaws),[NreferenceTraversal 1]));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);


% Setup
central_path = [-2 1; 1 4; 3 2; 5 2; 6 3; 7 2];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [-1 0.5; 0.5 2; 1.5 2; 3 0.5; 4 3; 5 4; 6 3; 7 1];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);
clear data
data.traversal{1} = nearby_traversal;
flag_yaw = 0;
flag_3D = 0;

Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [closestXs,closestYs,closestZs,closestYaws] = ...
        fcn_Path_findClosestPointsToTraversal(central_traversal,data,flag_yaw,flag_3D,([]));

end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [closestXs,closestYs,closestZs,closestYaws] = ...
        fcn_Path_findClosestPointsToTraversal(central_traversal,data,flag_yaw,flag_3D,(-1));

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
