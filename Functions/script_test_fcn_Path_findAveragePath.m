% script_test_fcn_Path_findAveragePath
% Tests the function: fcn_Path_findAveragePath

% Revision history
% 2020_11_10
% -- first write of the code
% 2021_01_07
% -- lots of bug fixes as we demo for the team (lol)
% 2020_01_09
% -- added more comments during clean-up
% 2021_12_27:
% -- corrected dependencies in comments
% 2022_01_03:
% -- corrected typos in comments
% -- fixed a bug where the Z value is not defined in loop
% 2022_01_06:
% -- refactored code, added weighted averaging to prevent iteration
% bouncing
% 2025_06_26 - Sean Brennan
% -- restructured test script to use function call for data loading
% -- added assertion testing
% 2025_07_01 - S. Brennan
% -- Removed traversal input type and replaced with cell array of paths
% -- Renamed function from fcn_Path_findAverageTraversalViaOrthoProjection


%      [path_average, closestXs, closestYs, closestDistances] = ...
%      fcn_Path_findAveragePath(...
%            data,
%            (reference_traversal),...
%            (num_iterations),
%            (weight_for_averaging),
%            (fig_num));

close all;


%% BASIC call: showing path averaging via ortho projection
fig_num = 10001;
titleString = sprintf('BASIC call: showing path averaging via ortho projection');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Load data
cellArrayOfPaths = fcn_INTERNAL_loadData;
stationInterval = 5;
max_num_iterations = 10;
exit_tolerance = []; averaging_weights = [];

% Find path average
% FORMAT:
%     [path_average, closestXs, closestYs, closestDistances] = ...
%         fcn_Path_findAveragePath(cellArrayOfPaths,...
%                 (stationInterval),...
%                 (max_num_iterations),...
%                 (exit_tolerance),(averaging_weights),......
%                 (averaging_weights),...
%                 (fig_num));
[path_average, closestXs, closestYs, closestDistances] = ...
    fcn_Path_findAveragePath(cellArrayOfPaths,...
    (stationInterval),...
    (max_num_iterations),...
    (exit_tolerance),(averaging_weights),......
    (fig_num));

title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(path_average));
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestDistances));

% Check variable sizes
NreferencePoints = 66;
assert(isequal(size(path_average),[NreferencePoints 2]));
assert(isequal(size(closestXs),[NreferencePoints length(cellArrayOfPaths)]));
assert(isequal(size(closestYs),[NreferencePoints length(cellArrayOfPaths)]));
assert(isequal(size(closestDistances),[NreferencePoints length(cellArrayOfPaths)]));

% Check  variable values
% Not possible - too many variables

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC call: adding more iterations
fig_num = 10002;
titleString = sprintf('BASIC call: adding more iterations');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Load data
cellArrayOfPaths = fcn_INTERNAL_loadData;
stationInterval = 1;
max_num_iterations = 50;
exit_tolerance = .001;
averaging_weights = [];

% Find path average
% FORMAT:
%     [path_average, closestXs, closestYs, closestDistances] = ...
%         fcn_Path_findAveragePath(cellArrayOfPaths,...
%                 (stationInterval),...
%                 (max_num_iterations),...
%                 (exit_tolerance),(averaging_weights),......
%                 (fig_num));
[path_average, closestXs, closestYs, closestDistances] = ...
    fcn_Path_findAveragePath(cellArrayOfPaths,...
    (stationInterval),...
    (max_num_iterations),...
    (exit_tolerance),(averaging_weights),......
    (fig_num));

title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(path_average));
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestDistances));

% Check variable sizes
NreferencePoints = 330;
assert(isequal(size(path_average),[NreferencePoints 2]));
assert(isequal(size(closestXs),[NreferencePoints length(cellArrayOfPaths)]));
assert(isequal(size(closestYs),[NreferencePoints length(cellArrayOfPaths)]));
assert(isequal(size(closestDistances),[NreferencePoints length(cellArrayOfPaths)]));

% Check  variable values
% Not possible - too many variables

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC call: show effect of weights
fig_num = 10003;
titleString = sprintf('BASIC call: show effect of weights');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Load data
clear cellArrayOfPaths
cellArrayOfPaths = cell(2,1);
cellArrayOfPaths{1} = [0 0; 10 0];
cellArrayOfPaths{2} = [1 1; 8 1];
stationInterval = 1;
max_num_iterations = 10;
exit_tolerance = []; 
averaging_weights = [];

% Find path average
% FORMAT:
%     [path_average, closestXs, closestYs, closestDistances] = ...
%         fcn_Path_findAveragePath(cellArrayOfPaths,...
%                 (stationInterval),...
%                 (max_num_iterations),...
%                 (exit_tolerance),(averaging_weights),......
%                 (averaging_weights),...
%                 (fig_num));
subplot(2,1,1)
[path_average, closestXs, closestYs, closestDistances] = ...
    fcn_Path_findAveragePath(cellArrayOfPaths,...
    (stationInterval),...
    (max_num_iterations),...
    (exit_tolerance),(averaging_weights),......
    (fig_num));

title('Default weights', 'Interpreter','none');

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(path_average));
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestDistances));

% Check variable sizes
NreferencePoints = 11;
assert(isequal(size(path_average),[NreferencePoints 2]));
assert(isequal(size(closestXs),[NreferencePoints length(cellArrayOfPaths)]));
assert(isequal(size(closestYs),[NreferencePoints length(cellArrayOfPaths)]));
assert(isequal(size(closestDistances),[NreferencePoints length(cellArrayOfPaths)]));

% Check  variable values
% Not possible - too many variables

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

subplot(2,1,2)
averaging_weights = [0.9; 0.1];
[path_average, closestXs, closestYs, closestDistances] = ...
    fcn_Path_findAveragePath(cellArrayOfPaths,...
    (stationInterval),...
    (max_num_iterations),...
    (exit_tolerance),(averaging_weights),......
    (fig_num));

title('90% weight on first', 'Interpreter','none');

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(path_average));
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestDistances));

% Check variable sizes
NreferencePoints = 11;
assert(isequal(size(path_average),[NreferencePoints 2]));
assert(isequal(size(closestXs),[NreferencePoints length(cellArrayOfPaths)]));
assert(isequal(size(closestYs),[NreferencePoints length(cellArrayOfPaths)]));
assert(isequal(size(closestDistances),[NreferencePoints length(cellArrayOfPaths)]));

% Check  variable values
% Not possible - too many variables

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

% Load data
cellArrayOfPaths = fcn_INTERNAL_loadData;
stationInterval = 5;
max_num_iterations = 10;
exit_tolerance = []; averaging_weights = [];

% Find path average
% FORMAT:
%     [path_average, closestXs, closestYs, closestDistances] = ...
%         fcn_Path_findAveragePath(cellArrayOfPaths,...
%                 (stationInterval),...
%                 (max_num_iterations),...
%                 (exit_tolerance),(averaging_weights),......
%                 (fig_num));
[path_average, closestXs, closestYs, closestDistances] = ...
    fcn_Path_findAveragePath(cellArrayOfPaths,...
    (stationInterval),...
    (max_num_iterations),...
    (exit_tolerance),(averaging_weights),......
    ([]));

% Check variable types
assert(isnumeric(path_average));
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestDistances));

% Check variable sizes
NreferencePoints = 66;
assert(isequal(size(path_average),[NreferencePoints 2]));
assert(isequal(size(closestXs),[NreferencePoints length(cellArrayOfPaths)]));
assert(isequal(size(closestYs),[NreferencePoints length(cellArrayOfPaths)]));
assert(isequal(size(closestDistances),[NreferencePoints length(cellArrayOfPaths)]));

% Check  variable values
% Not possible - too many variables

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

% Load data
cellArrayOfPaths = fcn_INTERNAL_loadData;
stationInterval = 5;
max_num_iterations = 10;
exit_tolerance = []; averaging_weights = [];

% Find path average
% FORMAT:
%     [path_average, closestXs, closestYs, closestDistances] = ...
%         fcn_Path_findAveragePath(cellArrayOfPaths,...
%                 (stationInterval),...
%                 (max_num_iterations),...
%                 (exit_tolerance),(averaging_weights),......
%                 (fig_num));
[path_average, closestXs, closestYs, closestDistances] = ...
    fcn_Path_findAveragePath(cellArrayOfPaths,...
    (stationInterval),...
    (max_num_iterations),...
    (exit_tolerance),(averaging_weights),......
    (-1));

% Check variable types
assert(isnumeric(path_average));
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestDistances));

% Check variable sizes
NreferencePoints = 66;
assert(isequal(size(path_average),[NreferencePoints 2]));
assert(isequal(size(closestXs),[NreferencePoints length(cellArrayOfPaths)]));
assert(isequal(size(closestYs),[NreferencePoints length(cellArrayOfPaths)]));
assert(isequal(size(closestDistances),[NreferencePoints length(cellArrayOfPaths)]));

% Check  variable values
% Not possible - too many variables

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

% Load data
cellArrayOfPaths = fcn_INTERNAL_loadData;
stationInterval = 5;
max_num_iterations = 10;
exit_tolerance = []; averaging_weights = [];

Niterations = 10;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [path_average, closestXs, closestYs, closestDistances] = ...
        fcn_Path_findAveragePath(cellArrayOfPaths,...
        (stationInterval),...
        (max_num_iterations),...
        (exit_tolerance),(averaging_weights),......
        ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [path_average, closestXs, closestYs, closestDistances] = ...
        fcn_Path_findAveragePath(cellArrayOfPaths,...
        (stationInterval),...
        (max_num_iterations),...
        (exit_tolerance),(averaging_weights),......
        (-1));
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

%% BUG case: found during HSOV testing, 18 laps in Reber parking lot (default)
fig_num = 90001;
titleString = sprintf('BUG case: found during HSOV testing, 18 laps in Reber parking lot (default)');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Load data - 18 laps in Reber parking lot
load('testData1_fcn_Path_findAveragePath.mat','paths');

cellArrayOfPaths = paths;
stationInterval = 1;
max_num_iterations = 20;
exit_tolerance = .01;
averaging_weights = [];

% Find path average
% FORMAT:
%     [path_average, closestXs, closestYs, closestDistances] = ...
%         fcn_Path_findAveragePath(cellArrayOfPaths,...
%                 (stationInterval),...
%                 (max_num_iterations),...
%                 (exit_tolerance),(averaging_weights),......
%                 (fig_num));
[path_average, closestXs, closestYs, closestDistances] = ...
    fcn_Path_findAveragePath(cellArrayOfPaths,...
    (stationInterval),...
    (max_num_iterations),...
    (exit_tolerance),(averaging_weights),......
    (fig_num));

title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(path_average));
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestDistances));

% Check variable sizes
NreferencePoints = 183;
assert(isequal(size(path_average),[NreferencePoints 2]));
assert(isequal(size(closestXs),[NreferencePoints length(cellArrayOfPaths)]));
assert(isequal(size(closestYs),[NreferencePoints length(cellArrayOfPaths)]));
assert(isequal(size(closestDistances),[NreferencePoints length(cellArrayOfPaths)]));

% Check  variable values
% Not possible - too many variables

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%%%%%%%%%%%%%%%%%%%%%%%%%
% BUG case: found during HSOV testing, 18 laps in Reber parking lot (iterative)
fig_num = 90002;
titleString = sprintf('BUG case: found during HSOV testing, 18 laps in Reber parking lot (iterative)');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Load data - 18 laps in Reber parking lot
load('testData1_fcn_Path_findAveragePath.mat','paths');

cellArrayOfPaths = paths;
stationInterval = 1;
max_num_iterations = 20;
exit_tolerance = .01;
averaging_weights = -1;

% Find path average
% FORMAT:
%     [path_average, closestXs, closestYs, closestDistances] = ...
%         fcn_Path_findAveragePath(cellArrayOfPaths,...
%                 (stationInterval),...
%                 (max_num_iterations),...
%                 (exit_tolerance),(averaging_weights),......
%                 (fig_num));
[path_average_iterative, closestXs, closestYs, closestDistances] = ...
    fcn_Path_findAveragePath(cellArrayOfPaths,...
    (stationInterval),...
    (max_num_iterations),...
    (exit_tolerance),(averaging_weights),......
    (fig_num));

title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(path_average_iterative));
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestDistances));

% Check variable sizes
NreferencePoints = 183;
assert(isequal(size(path_average_iterative),[NreferencePoints 2]));
assert(isequal(size(closestXs),[NreferencePoints length(cellArrayOfPaths)]));
assert(isequal(size(closestYs),[NreferencePoints length(cellArrayOfPaths)]));
assert(isequal(size(closestDistances),[NreferencePoints length(cellArrayOfPaths)]));

% Check  variable values
% Not possible - too many variables

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

% Compare results
if 1==1
    figure(38383);
    clf;
    plot(path_average(:,1),path_average(:,2),'b.-','LineWidth',1);
    hold on;
    plot(path_average_iterative(:,1),path_average_iterative(:,2),'r.-','LineWidth',1);
    axis equal
    
    % Calculate the difference in path points. These may be large as this
    % is determined by interpolation.
    difference = sum((path_average - path_average_iterative).^2,2).^0.5;
    maxDifference = max(difference);

    % Calculate the transverse distance error. This should be small, on
    % order of the tolerance (1 cm)
    St_points = fcn_Path_convertXY2St(path_average,path_average_iterative);
    maxTransverse = max(St_points(:,2));
    disp(maxTransverse);
end

%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 90003;
fprintf(1,'Figure: %.0f: Comparison of normal versus iterative timing\n',fig_num);
figure(fig_num);
close(fig_num);

% Load data - 18 laps in Reber parking lot
load('testData1_fcn_Path_findAveragePath.mat','paths');

cellArrayOfPaths = paths;
stationInterval = 1;
max_num_iterations = 20;
exit_tolerance = .01;

Niterations = 1;
averaging_weights = [];
% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [path_average, closestXs, closestYs, closestDistances] = ...
        fcn_Path_findAveragePath(cellArrayOfPaths,...
        (stationInterval),...
        (max_num_iterations),...
        (exit_tolerance),(averaging_weights),......
        (-1));
end
slow_method = toc;

averaging_weights = -1;
% Do calculation with pre-calculation, ITERATIVE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [path_average, closestXs, closestYs, closestDistances] = ...
        fcn_Path_findAveragePath(cellArrayOfPaths,...
        (stationInterval),...
        (max_num_iterations),...
        (exit_tolerance),(averaging_weights),......
        (-1));
end
fast_method = toc;

% Plot results as bar chart
figure(2343);
clf;
hold on;

X = categorical({'Normal mode','Iterative mode'});
X = reordercats(X,{'Normal mode','Iterative mode'}); % Forces bars to appear in this exact order, not alphabetized
Y = [slow_method fast_method ]/Niterations;
bar(X,Y)
ylabel('Execution time (Seconds)')


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


%% fcn_INTERNAL_loadData
function paths_array = fcn_INTERNAL_loadData

% Fill in sample paths (as a starter)
temp = fcn_Path_fillSamplePaths;
paths_array = cell(3,1);
for ith_path = 1:3
    paths_array{ith_path,1} = temp{ith_path};
end

% % Convert paths to traversal structures
% clear data
% for i_Path = 1:3
%     traversal = fcn_Path_convertPathToTraversalStructure(paths_array{i_Path}, -1);
%     data.traversal{i_Path} = traversal;
% end
end % Ends fcn_INTERNAL_loadData