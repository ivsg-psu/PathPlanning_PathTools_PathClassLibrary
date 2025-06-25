% script_test_fcn_Path_findOrthoScatterFromTraversalToTraversals.m
% Tests fcn_Path_findOrthoScatterFromTraversalToTraversals
       
% Revision history:
%      2021_01_02
%      -- first write of the code
%      2021_01_07:
%      -- renamed function to clarify paths versus traversals
%      2021_01_09
%      -- added more comments during clean-up
%      2022_01_03
%      -- found a bug in the constrainted search functionality
%      -- added assertion tests

close all;

%% BASIC example: 2 parallel paths with lower one as the central traversal
fig_num = 10001;
titleString = sprintf('BASIC example: 2 parallel paths with lower one as the central traversal');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

clear paths_array query_traversals all_traversals

paths_array{1} = [-2 0; 4 0; 8 0];
paths_array{2} = [0 2; 4 2; 6 2];


% Convert paths to traversal structures
for i_Path = 1:2
    traversal = ...
        fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    all_traversals.traversal{i_Path} = traversal;
end

% Pick the lower one to be the reference traversal
reference_traversal = all_traversals.traversal{1};

% Pick the upper one to be the query traversal
Nqueries = 0;
clear query_traversals
for i_Path = 2:length(paths_array)
    Nqueries = Nqueries+1;
    query_traversals.traversal{Nqueries} = all_traversals.traversal{i_Path};
end

% Set up the other function inputs
reference_station_points = [3; 7]; %(0:1:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 5;

[closestXs, closestYs, closestDistances] = fcn_Path_findOrthoScatterFromTraversalToTraversals(reference_station_points, reference_traversal, query_traversals, flag_rounding_type,search_radius,fig_num); %#ok<*ASGLU>

% Check variable types
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestDistances));

% Check variable sizes
NreferenceTraversal = length(reference_station_points(:,1));
Ntraversals = length(query_traversals.traversal);
assert(isequal(size(closestXs),[NreferenceTraversal Ntraversals]));
assert(isequal(size(closestYs),[NreferenceTraversal Ntraversals]));
assert(isequal(size(closestDistances),[NreferenceTraversal Ntraversals]));

% Check variable values
assert(isequal(round(closestXs,4),[ 1; 5]));
assert(isequal(round(closestYs,4),[ 2; 2]));
assert(isequal(round(closestDistances,4),[ 2; 2]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

% Peform histogram?
fcn_INTERNAL_plotHistogram(closestDistances)

%% BASIC test: 2 parallel paths with upper one as the central traversal
fig_num = 10002;
titleString = sprintf('BASIC test: 2 parallel paths with upper one as the central traversal');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

clear paths_array query_traversals all_traversals

paths_array{1} = [-2 0; 4 0; 8 0];
paths_array{2} = [0 -2; 4 -2; 6 -2];


% Convert paths to traversal structures
for i_Path = 1:2
    traversal = ...
        fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    all_traversals.traversal{i_Path} = traversal;
end

% Pick the upper one to be the reference traversal
reference_traversal = all_traversals.traversal{1};

% Pick the lower one to be the query traversal
Nqueries = 0;
clear query_traversals
for i_Path = 2:length(paths_array)
    Nqueries = Nqueries+1;
    query_traversals.traversal{Nqueries} = all_traversals.traversal{i_Path};
end

% Set up the other function inputs
reference_station_points = [3; 7]; %(0:1:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 5;

% Call the function
[closestXs, closestYs, closestDistances] = fcn_Path_findOrthoScatterFromTraversalToTraversals(reference_station_points, reference_traversal, query_traversals, flag_rounding_type,search_radius,fig_num); %#ok<*ASGLU>

% Check variable types
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestDistances));

% Check variable sizes
NreferenceTraversal = length(reference_station_points(:,1));
Ntraversals = length(query_traversals.traversal);
assert(isequal(size(closestXs),[NreferenceTraversal Ntraversals]));
assert(isequal(size(closestYs),[NreferenceTraversal Ntraversals]));
assert(isequal(size(closestDistances),[NreferenceTraversal Ntraversals]));

% Check variable values
assert(isequal(round(closestXs,4),[ 1; 5]));
assert(isequal(round(closestYs,4),[ -2; -2]));
assert(isequal(round(closestDistances,4),[ -2; -2]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

% Peform histogram?
fcn_INTERNAL_plotHistogram(closestDistances)

%% BASIC test: 3 parallel paths
fig_num = 10003;
titleString = sprintf('BASIC test: 3 parallel paths');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

clear paths_array all_traversals

paths_array{1} = [0 0; 4 0; 6 0];
paths_array{2} = [0 2; 4 2; 6 2];
paths_array{3} = [0 4; 4 4; 6 4];


% Convert paths to traversal structures
for i_Path = 1:3
    traversal = ...
        fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    all_traversals.traversal{i_Path} = traversal;
end

reference_traversal = all_traversals.traversal{2};
reference_station_points = [1; 5]; %(0:1:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 5;

% Call the function
[closestXs, closestYs, closestDistances] = fcn_Path_findOrthoScatterFromTraversalToTraversals(reference_station_points, reference_traversal, all_traversals, flag_rounding_type,search_radius,fig_num); %#ok<*ASGLU>

% Check variable types
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestDistances));

% Check variable sizes
NreferenceTraversal = length(reference_station_points(:,1));
Ntraversals = length(all_traversals.traversal);
assert(isequal(size(closestXs),[NreferenceTraversal Ntraversals]));
assert(isequal(size(closestYs),[NreferenceTraversal Ntraversals]));
assert(isequal(size(closestDistances),[NreferenceTraversal Ntraversals]));

% Check variable values
assert(isequal(round(closestXs,4),[ 1     1     1;   5     5     5]));
assert(isequal(round(closestYs,4),[ 0     2     4;   0     2     4]));
assert(isequal(round(closestDistances,4),[ -2     0     2;   -2     0     2]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

% Peform histogram?
fcn_INTERNAL_plotHistogram(closestDistances)


%% BASIC test: 4 parallel paths, one outside search radius
fig_num = 10004;
titleString = sprintf('BASIC test: 4 parallel paths, one outside search radius');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

clear paths_array all_traversals

paths_array{1} = [0 0; 4 0; 6 0];
paths_array{2} = [0 2; 4 2; 6 2];
paths_array{3} = [0 4; 4 4; 6 4];
paths_array{4} = [0 10; 4 10; 6 10];


% Convert paths to traversal structures
for i_Path = 1:3
    traversal = ...
        fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    all_traversals.traversal{i_Path} = traversal;
end

% Define the inputs
reference_traversal = all_traversals.traversal{2};
reference_station_points = [1; 5]; %(0:1:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 5;

% Call the function
[closestXs, closestYs, closestDistances] = fcn_Path_findOrthoScatterFromTraversalToTraversals(reference_station_points, reference_traversal, all_traversals, flag_rounding_type,search_radius,fig_num); %#ok<*ASGLU>

% Check variable types
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestDistances));

% Check variable sizes
NreferenceTraversal = length(reference_station_points(:,1));
Ntraversals = length(all_traversals.traversal);
assert(isequal(size(closestXs),[NreferenceTraversal Ntraversals]));
assert(isequal(size(closestYs),[NreferenceTraversal Ntraversals]));
assert(isequal(size(closestDistances),[NreferenceTraversal Ntraversals]));

% Check variable values
% assert(isequaln(round(closestXs,4),[ 1     1     1    NaN;   5     5     5    NaN]));
% assert(isequaln(round(closestYs,4),[ 0     2     4    NaN;   0     2     4    NaN]));
% assert(isequaln(round(closestDistances,4),[-2     0     2    NaN;  -2     0     2    NaN]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

% Peform histogram?
fcn_INTERNAL_plotHistogram(closestDistances)



%% REAL example: basic call
fig_num = 20001;
titleString = sprintf('REAL example: basic call');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

clear paths_array all_traversals

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
for i_Path = 1:3
    traversal = ...
        fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    all_traversals.traversal{i_Path} = traversal;
end

reference_traversal = all_traversals.traversal{2};
reference_station_points = (0:10:reference_traversal.Station(end))';
flag_rounding_type = 3; % Use average of projections at end points
search_radius = 7;

[closestXs, closestYs, closestDistances] = fcn_Path_findOrthoScatterFromTraversalToTraversals( reference_station_points, reference_traversal, all_traversals, flag_rounding_type,search_radius,fig_num); 


% Check variable types
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestDistances));

% Check variable sizes
NreferenceTraversal = length(reference_station_points(:,1));
Ntraversals = length(all_traversals.traversal);
assert(isequal(size(closestXs),[NreferenceTraversal Ntraversals]));
assert(isequal(size(closestYs),[NreferenceTraversal Ntraversals]));
assert(isequal(size(closestDistances),[NreferenceTraversal Ntraversals]));

% Check variable values
% assert(isequal(closestXs,[0; 2]));
% assert(isequal(closestYs,[4; 4]));
% assert(isequal(closestDistances,[0; 0]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

% Peform histogram?
fcn_INTERNAL_plotHistogram(closestDistances)

%% REAL example: basic call with finer resolution
fig_num = 20002;
titleString = sprintf('REAL example: basic call with finer resolution');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

clear paths_array all_traversals

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
for i_Path = 1:3
    traversal = ...
        fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    all_traversals.traversal{i_Path} = traversal;
end

reference_traversal = all_traversals.traversal{2};
reference_station_points = (0:1:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 7;

[closestXs, closestYs, closestDistances] = ...
       fcn_Path_findOrthoScatterFromTraversalToTraversals(...
       reference_station_points, reference_traversal, all_traversals,...
       flag_rounding_type,search_radius,fig_num);


% Check variable types
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestDistances));

% Check variable sizes
NreferenceTraversal = length(reference_station_points(:,1));
Ntraversals = length(all_traversals.traversal);
assert(isequal(size(closestXs),[NreferenceTraversal Ntraversals]));
assert(isequal(size(closestYs),[NreferenceTraversal Ntraversals]));
assert(isequal(size(closestDistances),[NreferenceTraversal Ntraversals]));

% Check variable values
% assert(isequal(closestXs,[0; 2]));
% assert(isequal(closestYs,[4; 4]));
% assert(isequal(closestDistances,[0; 0]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

% Peform histogram?
fcn_INTERNAL_plotHistogram(closestDistances)

%% REAL example: basic call but with random weird traversal nearby
% Used to teach the search radius criteria

fig_num = 20003;
titleString = sprintf('REAL example: basic call but with random weird traversal nearby');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

weird_path = [100 20; 90 50; 100 70; 100 90];

clear paths_array all_traversals

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
for i_Path = 1:3
    traversal = ...
        fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    all_traversals.traversal{i_Path} = traversal;
end

traversal = ...
    fcn_Path_convertPathToTraversalStructure(weird_path);

all_traversals.traversal{end+1} = traversal;


reference_traversal = all_traversals.traversal{2};
reference_station_points = (0:10:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 10;

[closestXs, closestYs, closestDistances] = ...
       fcn_Path_findOrthoScatterFromTraversalToTraversals(...
       reference_station_points, reference_traversal, all_traversals,...
       flag_rounding_type,search_radius,fig_num); %#ok<*ASGLU>

% Check variable types
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestDistances));

% Check variable sizes
NreferenceTraversal = length(reference_station_points(:,1));
Ntraversals = length(all_traversals.traversal);
assert(isequal(size(closestXs),[NreferenceTraversal Ntraversals]));
assert(isequal(size(closestYs),[NreferenceTraversal Ntraversals]));
assert(isequal(size(closestDistances),[NreferenceTraversal Ntraversals]));

% Check variable values
% assert(isequal(closestXs,[0; 2]));
% assert(isequal(closestYs,[4; 4]));
% assert(isequal(closestDistances,[0; 0]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

% Peform histogram?
fcn_INTERNAL_plotHistogram(closestDistances)


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

%% BASIC example - NO FIGURE
fig_num = 80001;
fprintf(1,'Figure: %.0f: Demo of fast mode, empty fig_num\n',fig_num);
figure(fig_num); close(fig_num);

clear paths_array all_traversals

paths_array{1} = [0 0; 4 0; 6 0];
paths_array{2} = [0 2; 4 2; 6 2];
paths_array{3} = [0 4; 4 4; 6 4];


% Convert paths to traversal structures
for i_Path = 1:3
    traversal = ...
        fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    all_traversals.traversal{i_Path} = traversal;
end

reference_traversal = all_traversals.traversal{2};
reference_station_points = [1; 5]; %(0:1:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 5;

% Call the function
[closestXs, closestYs, closestDistances] = ...
    fcn_Path_findOrthoScatterFromTraversalToTraversals(reference_station_points, reference_traversal, all_traversals, ...
    flag_rounding_type,search_radius,([])); 

% Check variable types
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestDistances));

% Check variable sizes
NreferenceTraversal = length(reference_station_points(:,1));
Ntraversals = length(all_traversals.traversal);
assert(isequal(size(closestXs),[NreferenceTraversal Ntraversals]));
assert(isequal(size(closestYs),[NreferenceTraversal Ntraversals]));
assert(isequal(size(closestDistances),[NreferenceTraversal Ntraversals]));

% Check variable values
assert(isequal(round(closestXs,4),[ 1     1     1;   5     5     5]));
assert(isequal(round(closestYs,4),[ 0     2     4;   0     2     4]));
assert(isequal(round(closestDistances,4),[ -2     0     2;   -2     0     2]));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% BASIC fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);


clear paths_array all_traversals

paths_array{1} = [0 0; 4 0; 6 0];
paths_array{2} = [0 2; 4 2; 6 2];
paths_array{3} = [0 4; 4 4; 6 4];


% Convert paths to traversal structures
for i_Path = 1:3
    traversal = ...
        fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    all_traversals.traversal{i_Path} = traversal;
end

reference_traversal = all_traversals.traversal{2};
reference_station_points = [1; 5]; %(0:1:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 5;

% Call the function
[closestXs, closestYs, closestDistances] = ...
    fcn_Path_findOrthoScatterFromTraversalToTraversals(reference_station_points, reference_traversal, all_traversals, ...
    flag_rounding_type,search_radius,(-1)); 

% Check variable types
assert(isnumeric(closestXs));
assert(isnumeric(closestYs));
assert(isnumeric(closestDistances));

% Check variable sizes
NreferenceTraversal = length(reference_station_points(:,1));
Ntraversals = length(all_traversals.traversal);
assert(isequal(size(closestXs),[NreferenceTraversal Ntraversals]));
assert(isequal(size(closestYs),[NreferenceTraversal Ntraversals]));
assert(isequal(size(closestDistances),[NreferenceTraversal Ntraversals]));

% Check variable values
assert(isequal(round(closestXs,4),[ 1     1     1;   5     5     5]));
assert(isequal(round(closestYs,4),[ 0     2     4;   0     2     4]));
assert(isequal(round(closestDistances,4),[ -2     0     2;   -2     0     2]));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

clear paths_array all_traversals

paths_array{1} = [0 0; 4 0; 6 0];
paths_array{2} = [0 2; 4 2; 6 2];
paths_array{3} = [0 4; 4 4; 6 4];


% Convert paths to traversal structures
for i_Path = 1:3
    traversal = ...
        fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    all_traversals.traversal{i_Path} = traversal;
end

reference_traversal = all_traversals.traversal{2};
reference_station_points = [1; 5]; %(0:1:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 5;

Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [closestXs, closestYs, closestDistances] = ...
        fcn_Path_findOrthoScatterFromTraversalToTraversals(reference_station_points, reference_traversal, all_traversals, ...
        flag_rounding_type,search_radius,([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [closestXs, closestYs, closestDistances] = ...
        fcn_Path_findOrthoScatterFromTraversalToTraversals(reference_station_points, reference_traversal, all_traversals, ...
        flag_rounding_type,search_radius,(-1));
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


function fcn_INTERNAL_plotHistogram(closestDistances)
Nrows = length(closestDistances(:,1));
Ncols = length(closestDistances(1,:));
allDistances = reshape(closestDistances,[Nrows*Ncols 1]);

figure(111);
histogram(allDistances,30);
title('Histogram of all orthogonal distance projections');
end