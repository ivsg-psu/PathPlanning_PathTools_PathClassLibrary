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

%% Basic test 1, 2 parallel paths with lower one as the central traversal
clear paths_array query_traversals all_traversals

paths_array{1} = [-2 0; 4 0; 8 0];
paths_array{2} = [0 2; 4 2; 6 2];


% Convert paths to traversal structures
for i_Path = 1:length(paths_array)
    traversal = ...
        fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    all_traversals.traversal{i_Path} = traversal;
end

% Plot the results? (Note: they are plotted below as well)
if 1==0
    fig_num = 12;
    fcn_Path_plotTraversalsYaw(all_traversals,fig_num);
    fig_num = 13;
    fcn_Path_plotTraversalsXY(all_traversals,fig_num);
end

% Pick the lower one to be the reference traversal
reference_traversal = all_traversals.traversal{1};

% Pick the upper one to be the query traversal
Nqueries = 0;
for i_Path = 2:length(paths_array)
    Nqueries = Nqueries+1;
    query_traversals.traversal{Nqueries} = all_traversals.traversal{i_Path};
end

% Set up the other function inputs
reference_station_points = [3; 7]; %(0:1:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 5;
fig_num = 1;

[closestXs, closestYs, closestDistances] = fcn_Path_findOrthoScatterFromTraversalToTraversals(reference_station_points, reference_traversal, query_traversals, flag_rounding_type,search_radius,fig_num); %#ok<*ASGLU>

% Make sure function worked!
assert(isequal(round(closestXs,4),[ 1; 5]));
assert(isequal(round(closestYs,4),[ 2; 2]));
assert(isequal(round(closestDistances,4),[ 2; 2]));

% Peform histogram
Nrows = length(closestDistances(:,1));
Ncols = length(closestDistances(1,:));
allDistances = reshape(closestDistances,[Nrows*Ncols 1]);

figure(111);
histogram(allDistances,30);
title('Histogram of all orthogonal distance projections');


%% Basic test 2, 2 parallel paths with upper one as the central traversal
clear paths_array query_traversals all_traversals

paths_array{1} = [-2 0; 4 0; 8 0];
paths_array{2} = [0 -2; 4 -2; 6 -2];


% Convert paths to traversal structures
for i_Path = 1:length(paths_array)
    traversal = ...
        fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    all_traversals.traversal{i_Path} = traversal;
end

% Pick the upper one to be the reference traversal
reference_traversal = all_traversals.traversal{1};

% Pick the lower one to be the query traversal
Nqueries = 0;
for i_Path = 2:length(paths_array)
    Nqueries = Nqueries+1;
    query_traversals.traversal{Nqueries} = all_traversals.traversal{i_Path};
end

% Set up the other function inputs
reference_station_points = [3; 7]; %(0:1:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 5;
fig_num = 1;

% Call the function
[closestXs, closestYs, closestDistances] = fcn_Path_findOrthoScatterFromTraversalToTraversals(reference_station_points, reference_traversal, query_traversals, flag_rounding_type,search_radius,fig_num); %#ok<*ASGLU>

% Make sure function worked!
assert(isequal(round(closestXs,4),[ 1; 5]));
assert(isequal(round(closestYs,4),[ -2; -2]));
assert(isequal(round(closestDistances,4),[ -2; -2]));

% Peform histogram
Nrows = length(closestDistances(:,1));
Ncols = length(closestDistances(1,:));
allDistances = reshape(closestDistances,[Nrows*Ncols 1]);

figure(111);
histogram(allDistances,30);
title('Histogram of all orthogonal distance projections');



%% Basic test 1, 3 parallel paths
clear paths_array all_traversals

paths_array{1} = [0 0; 4 0; 6 0];
paths_array{2} = [0 2; 4 2; 6 2];
paths_array{3} = [0 4; 4 4; 6 4];


% Convert paths to traversal structures
for i_Path = 1:length(paths_array)
    traversal = ...
        fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    all_traversals.traversal{i_Path} = traversal;
end

reference_traversal = all_traversals.traversal{2};
reference_station_points = [1; 5]; %(0:1:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 5;
fig_num = 1;

% Call the function
[closestXs, closestYs, closestDistances] = fcn_Path_findOrthoScatterFromTraversalToTraversals(reference_station_points, reference_traversal, all_traversals, flag_rounding_type,search_radius,fig_num); %#ok<*ASGLU>

% Make sure function worked!
assert(isequal(round(closestXs,4),[ 1     1     1;   5     5     5]));
assert(isequal(round(closestYs,4),[ 0     2     4;   0     2     4]));
assert(isequal(round(closestDistances,4),[ -2     0     2;   -2     0     2]));

% Peform histogram
Nrows = length(closestDistances(:,1));
Ncols = length(closestDistances(1,:));
allDistances = reshape(closestDistances,[Nrows*Ncols 1]);

figure(111);
histogram(allDistances,30);
title('Histogram of all orthogonal distance projections');


%% Basic test 2, 4 parallel paths, one outside search radius
clear paths_array all_traversals

paths_array{1} = [0 0; 4 0; 6 0];
paths_array{2} = [0 2; 4 2; 6 2];
paths_array{3} = [0 4; 4 4; 6 4];
paths_array{4} = [0 10; 4 10; 6 10];


% Convert paths to traversal structures
for i_Path = 1:length(paths_array)
    traversal = ...
        fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    all_traversals.traversal{i_Path} = traversal;
end

% Define the inputs
reference_traversal = all_traversals.traversal{2};
reference_station_points = [1; 5]; %(0:1:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 5;
fig_num = 1;

% Call the function
[closestXs, closestYs, closestDistances] = fcn_Path_findOrthoScatterFromTraversalToTraversals(reference_station_points, reference_traversal, all_traversals, flag_rounding_type,search_radius,fig_num); %#ok<*ASGLU>

% Make sure function worked!
assert(isequaln(round(closestXs,4),[ 1     1     1    NaN;   5     5     5    NaN]));
assert(isequaln(round(closestYs,4),[ 0     2     4    NaN;   0     2     4    NaN]));
assert(isequaln(round(closestDistances,4),[-2     0     2    NaN;  -2     0     2    NaN]));

% Peform histogram
Nrows = length(closestDistances(:,1));
Ncols = length(closestDistances(1,:));
allDistances = reshape(closestDistances,[Nrows*Ncols 1]);

figure(111);
histogram(allDistances,30);
title('Histogram of all orthogonal distance projections');



%% Prep for test cases
clear paths_array all_traversals

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
for i_Path = 1:length(paths_array)
    traversal = ...
        fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    all_traversals.traversal{i_Path} = traversal;
end

% Plot the results? (Note: they are plotted below as well)
if 1==0
    fig_num = 12;
    fcn_Path_plotTraversalsYaw(all_traversals,fig_num);
    fig_num = 13;
    fcn_Path_plotTraversalsXY(all_traversals,fig_num);
end

%% Test case 1: basic call
clear paths_array all_traversals

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
for i_Path = 1:length(paths_array)
    traversal = ...
        fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    all_traversals.traversal{i_Path} = traversal;
end

reference_traversal = all_traversals.traversal{2};
reference_station_points = (0:10:reference_traversal.Station(end))';
flag_rounding_type = 3; % Use average of projections at end points
search_radius = 7;
fig_num = 1;

[closestXs, closestYs, closestDistances] = fcn_Path_findOrthoScatterFromTraversalToTraversals( reference_station_points, reference_traversal, all_traversals, flag_rounding_type,search_radius,fig_num); 
   
figure(11);
histogram([closestDistances(:,1);closestDistances(:,3)],30);
title('Histogram of all orthogonal distance projections');


%% Test case 2: basic call with finer resolution
clear paths_array all_traversals

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
for i_Path = 1:length(paths_array)
    traversal = ...
        fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    all_traversals.traversal{i_Path} = traversal;
end

reference_traversal = all_traversals.traversal{2};
reference_station_points = (0:1:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 7;
fig_num = 2;

[closestXs, closestYs, closestDistances] = ...
       fcn_Path_findOrthoScatterFromTraversalToTraversals(...
       reference_station_points, reference_traversal, all_traversals,...
       flag_rounding_type,search_radius,fig_num);
   
figure(22);
histogram([closestDistances(:,1);closestDistances(:,3)],30);
title('Histogram of all orthogonal distance projections');


%% Test case 3: basic call but with random weird traversal nearby
% Used to teach the search radius criteria
weird_path = [100 20; 90 50; 100 70; 100 90];

clear paths_array all_traversals

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
for i_Path = 1:length(paths_array)
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
fig_num = 3;

[closestXs, closestYs, closestDistances] = ...
       fcn_Path_findOrthoScatterFromTraversalToTraversals(...
       reference_station_points, reference_traversal, all_traversals,...
       flag_rounding_type,search_radius,fig_num); %#ok<*ASGLU>
   
figure(33);
histogram([closestDistances(:,1);closestDistances(:,3)],30);
title('Histogram of all orthogonal distance projections');

