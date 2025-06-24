% script_test_fcn_Path_fillOffsetTraversalsAboutTraversal
% Tests fcn_Path_fillOffsetTraversalsAboutTraversal
       
% Revision history:
%      2021_01_24
%      -- first write of the code


close all


%% Test case: basic call for one trajectory
% NOTE: the function itself does not plot since not given a figure number
fig_num = 10001;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

offsets = 2; 
offset_traversals = fcn_Path_fillOffsetTraversalsAboutTraversal(reference_traversal, offsets, [], fig_num);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Test case: basic call for one trajectory - specify figure
fig_num = 10002;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

offsets = 2; 
flag_rounding_type = [];
offset_traversals = fcn_Path_fillOffsetTraversalsAboutTraversal(reference_traversal, offsets, flag_rounding_type, fig_num); %#ok<*NASGU>

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Test case: basic call for two trajectories 
fig_num = 10003;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;


% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

offsets = [2; -2]; 
flag_rounding_type = [];
offset_traversals = fcn_Path_fillOffsetTraversalsAboutTraversal(reference_traversal, offsets,  flag_rounding_type, fig_num);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Test case: show how "pinching" can happen
fig_num = 10004;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;


% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

flag_rounding_type = [];
% Grab the "curve" of the path
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1}(13:20,:));
offsets = (-10:1:10)'; 
offset_traversals = fcn_Path_fillOffsetTraversalsAboutTraversal(reference_traversal, offsets, flag_rounding_type,  fig_num);
axis equal;

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Test case: Show how offsets can link lane markers
fig_num = 10005;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;


% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

flag_rounding_type = [];
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1}(30:end,:));
path2 = [30 10; 25 44];
second_traversal = fcn_Path_convertPathToTraversalStructure(path2);

offsets = 2; 
offset_traversal_1 = fcn_Path_fillOffsetTraversalsAboutTraversal(reference_traversal, offsets, flag_rounding_type,  fig_num);

offsets = [2; -2]; 
offset_traversal_2 = fcn_Path_fillOffsetTraversalsAboutTraversal(second_traversal, offsets, flag_rounding_type,  fig_num);

% Find intersections
[right_intersection_points,...
    ~,...
    ~] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    offset_traversal_1.traversal{1},...
    offset_traversal_2.traversal{1});
plot(right_intersection_points(:,1),right_intersection_points(:,2),'ro','Markersize',10);

[left_intersection_points,...
    ~,...
    ~] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    offset_traversal_1.traversal{1},...
    offset_traversal_2.traversal{2});
plot(left_intersection_points(:,1),left_intersection_points(:,2),'bo','Markersize',10);
title('Illustration of how to use offsets to link lane edge designations');

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Demonstration of effect of flag_rounding_type
fig_num = 10006;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

angles = (-45:45)'*pi/180;
path = 20*[cos(angles) sin(angles)];

reference_traversal = fcn_Path_convertPathToTraversalStructure(path);

offsets = [2; -2]; 

for flag_rounding_type = 1:4
    subplot(2,2,flag_rounding_type);
    hold on;
    grid on;
    axis equal;
    offset_traversals = fcn_Path_fillOffsetTraversalsAboutTraversal(reference_traversal, offsets, flag_rounding_type,  fig_num);
    title(sprintf('flag_rounding_type: %.0d',flag_rounding_type),'Interpreter','none');
end
sgtitle('Effect of flag_rounding_type','Interpreter','none');


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

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

offsets = [2; -2]; 
flag_rounding_type = [];
offset_traversals = fcn_Path_fillOffsetTraversalsAboutTraversal(reference_traversal, offsets,  flag_rounding_type, []);

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

offsets = [2; -2]; 
flag_rounding_type = [];
offset_traversals = fcn_Path_fillOffsetTraversalsAboutTraversal(reference_traversal, offsets,  flag_rounding_type, -1);

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

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

offsets = [2; -2]; 
flag_rounding_type = [];

Niterations = 20;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    offset_traversals = fcn_Path_fillOffsetTraversalsAboutTraversal(reference_traversal, offsets,  flag_rounding_type, []);

end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    offset_traversals = fcn_Path_fillOffsetTraversalsAboutTraversal(reference_traversal, offsets,  flag_rounding_type, -1);

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
