% script_test_fcn_fcn_Path_calcYawFromPathSegments
% Tests fcn_Path_calcYawFromPathSegments
       
% Revision history:
%      2021_01_06
%      -- first write of the code
%      2021_01_07
%      -- updated function calls to reflect paths vs traversals
%      -- minor comment clean-ups

close all

%% Test case 1: Basic call with one path
fig_num = 10001;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;
% Pick first path as reference_traversal structure
path_to_check = paths_array{1};

yaw_angles = fcn_Path_calcYawFromPathSegments(path_to_check,fig_num); 

% Check the variable type
assert(isnumeric(yaw_angles));

% Check the variable size
assert(isequal(size(yaw_angles),[length(path_to_check(:,1))-1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Test case 2: Multiple paths
fig_num = 10002;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;
figure(fig_num*10); clf;


% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;
clear all_traversals
for ith_path = 1:3
    % Pick first path as reference_traversal structure
    path_to_check = paths_array{ith_path};
    
    % Pick first path as reference_traversal structure
    traversal_to_check = fcn_Path_convertPathToTraversalStructure(path_to_check);
    all_traversals.traversal{1} = traversal_to_check;
    fcn_Path_plotTraversalsXY(all_traversals,fig_num*10);

    yaw_angles = fcn_Path_calcYawFromPathSegments(path_to_check,fig_num);

    % Check the variable type
    assert(isnumeric(yaw_angles));

    % Check the variable size
    assert(isequal(size(yaw_angles),[length(path_to_check(:,1))-1 1]));
end

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Test case 3: basic call with degenerate path
fig_num = 10003;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;


path_to_check = [1 1; 0 0];
yaw_angles = fcn_Path_calcYawFromPathSegments(path_to_check,fig_num);

% Check the variable type
assert(isnumeric(yaw_angles));

% Check the variable size
assert(isequal(size(yaw_angles),[length(path_to_check(:,1))-1 1]));

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
% Pick first path as reference_traversal structure
path_to_check = paths_array{1};

yaw_angles = fcn_Path_calcYawFromPathSegments(path_to_check,[]); 

% Check the variable type
assert(isnumeric(yaw_angles));

% Check the variable size
assert(isequal(size(yaw_angles),[length(path_to_check(:,1))-1 1]));


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;
% Pick first path as reference_traversal structure
path_to_check = paths_array{1};

yaw_angles = fcn_Path_calcYawFromPathSegments(path_to_check,-1); 

% Check the variable type
assert(isnumeric(yaw_angles));

% Check the variable size
assert(isequal(size(yaw_angles),[length(path_to_check(:,1))-1 1]));


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
% Pick first path as reference_traversal structure
path_to_check = paths_array{1};


Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    yaw_angles = fcn_Path_calcYawFromPathSegments(path_to_check,[]); 
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    yaw_angles = fcn_Path_calcYawFromPathSegments(path_to_check,-1); 
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

