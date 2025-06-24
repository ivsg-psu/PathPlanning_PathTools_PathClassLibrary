% script_test_fcn_Path_convertTraversalXYtoSy.m
% Tests script_test_fcn_Path_convertTraversalXYtoSy

% Revision history:
%      2021_03_21
%      -- first write of the code (from
%      script_test_fcn_Path_findOrthoScatterFromTraversalToTraversals)


close all;


%% Test case 1: basic call
fig_num = 10001;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
clear all_traversals
for i_Path = 1:length(paths_array)
    traversal = ...
        fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    all_traversals.traversal{i_Path} = traversal;
end

reference_traversal = all_traversals.traversal{2};
reference_station_points = (0:30:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 40;

[closestXs, closestYs, closestDistances] = ...
    fcn_Path_convertTraversalXYtoSy(...
    reference_station_points, reference_traversal, all_traversals,...
    flag_rounding_type,search_radius,fig_num); %#ok<*ASGLU>

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Test case 2: basic call with different reference
fig_num = 10002;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
clear all_traversals
for i_Path = 1:length(paths_array)
    traversal = ...
        fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    all_traversals.traversal{i_Path} = traversal;
end

reference_traversal = all_traversals.traversal{1};
reference_station_points = (0:30:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 40;

[closestXs, closestYs, closestDistances] = ...
    fcn_Path_convertTraversalXYtoSy(...
    reference_station_points, reference_traversal, all_traversals,...
    flag_rounding_type,search_radius,fig_num); %#ok<*ASGLU>

reference_traversal = all_traversals.traversal{3};
reference_station_points = (0:10:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 40;

[closestXs, closestYs, closestDistances] = ...
    fcn_Path_convertTraversalXYtoSy(...
    reference_station_points, reference_traversal, all_traversals,...
    flag_rounding_type,search_radius,fig_num); %#ok<*ASGLU>


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Test case 3: basic call 1 with finer resolution
fig_num = 10003;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
clear all_traversals
for i_Path = 1:length(paths_array)
    traversal = ...
        fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    all_traversals.traversal{i_Path} = traversal;
end
reference_traversal = all_traversals.traversal{2};
reference_station_points = (0:10:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 40;

[closestXs, closestYs, closestDistances] = ...
    fcn_Path_convertTraversalXYtoSy(...
    reference_station_points, reference_traversal, all_traversals,...
    flag_rounding_type,search_radius,fig_num);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Advanced testing example of fcn_Path_convertTraversalXYtoSy
fig_num = 10004;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Set up data
reference_path = [0 0; 1 1; 2 0];
reference_traversal = fcn_Path_convertPathToTraversalStructure(reference_path);
% stations = [0; 1; 2^0.5-0.1; 2^0.5; 2^0.5+.1; 2; central_traversal.Station(end)];
reference_station_points_1st_half = linspace(0,sum(reference_path(2,:).^2,2).^0.5,11)'; 
reference_station_points_2nd_half = linspace(sum(reference_path(2,:).^2,2).^0.5,reference_traversal.Station(end),10)';
reference_station_points = [reference_station_points_1st_half(1:end-1,:);reference_station_points_2nd_half];

clear data 
% Load a test path that is challenging for this reference path
test_path = fcn_Path_fillSamplePaths;
test_path = test_path{4};
test_traversal.traversal{1} = fcn_Path_convertPathToTraversalStructure(test_path);

flag_rounding_type = 3;
search_radius = 40;

[closestXs, closestYs, closestDistances] = ...
    fcn_Path_convertTraversalXYtoSy(...
    reference_station_points, reference_traversal, test_traversal,...
    flag_rounding_type,search_radius,fig_num);

% Repeat with different projection type
flag_rounding_type = 4;
search_radius = 40;

[closestXs, closestYs, closestDistances] = ...
    fcn_Path_convertTraversalXYtoSy(...
    reference_station_points, reference_traversal, test_traversal,...
    flag_rounding_type,search_radius,fig_num);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Another advanced test, with crossings
fig_num = 10005;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Set up data
reference_path = [0 0; 1 2; 2 0];
reference_traversal = fcn_Path_convertPathToTraversalStructure(reference_path);
% stations = [0; 1; 2^0.5-0.1; 2^0.5; 2^0.5+.1; 2; central_traversal.Station(end)];
reference_station_points_1st_half = linspace(0,sum(reference_path(2,:).^2,2).^0.5,11)'; 
reference_station_points_2nd_half = linspace(sum(reference_path(2,:).^2,2).^0.5,reference_traversal.Station(end),10)';
reference_station_points = [reference_station_points_1st_half(1:end-1,:);reference_station_points_2nd_half];

clear data
% Load a test path that is challenging for this reference path
test_path = fcn_Path_fillSamplePaths;
test_path = test_path{4};
test_traversal.traversal{1} = fcn_Path_convertPathToTraversalStructure(test_path);

flag_rounding_type = 3;
search_radius = 40;

[closestXs, closestYs, closestDistances] = ...
    fcn_Path_convertTraversalXYtoSy(...
    reference_station_points, reference_traversal, test_traversal,...
    flag_rounding_type,search_radius,fig_num);

% Repeat with different projection type
flag_rounding_type = 4;
search_radius = 40;

[closestXs, closestYs, closestDistances] = ...
    fcn_Path_convertTraversalXYtoSy(...
    reference_station_points, reference_traversal, test_traversal,...
    flag_rounding_type,search_radius,fig_num);

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

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
clear all_traversals
for i_Path = 1:length(paths_array)
    traversal = ...
        fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    all_traversals.traversal{i_Path} = traversal;
end

reference_traversal = all_traversals.traversal{2};
reference_station_points = (0:30:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 40;

[closestXs, closestYs, closestDistances] = ...
    fcn_Path_convertTraversalXYtoSy(...
    reference_station_points, reference_traversal, all_traversals,...
    flag_rounding_type,search_radius,[]);

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);


% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
clear all_traversals
for i_Path = 1:length(paths_array)
    traversal = ...
        fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    all_traversals.traversal{i_Path} = traversal;
end

reference_traversal = all_traversals.traversal{2};
reference_station_points = (0:30:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 40;

[closestXs, closestYs, closestDistances] = ...
    fcn_Path_convertTraversalXYtoSy(...
    reference_station_points, reference_traversal, all_traversals,...
    flag_rounding_type,search_radius,-1);

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
clear all_traversals
for i_Path = 1:length(paths_array)
    traversal = ...
        fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    all_traversals.traversal{i_Path} = traversal;
end

reference_traversal = all_traversals.traversal{2};
reference_station_points = (0:30:reference_traversal.Station(end))';
flag_rounding_type = 3;
search_radius = 40;



Niterations = 20;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [closestXs, closestYs, closestDistances] = ...
        fcn_Path_convertTraversalXYtoSy(...
        reference_station_points, reference_traversal, all_traversals,...
        flag_rounding_type,search_radius,[]);
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [closestXs, closestYs, closestDistances] = ...
        fcn_Path_convertTraversalXYtoSy(...
        reference_station_points, reference_traversal, all_traversals,...
        flag_rounding_type,search_radius,-1);
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
