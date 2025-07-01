% script_test_fcn_Path_plotPathXYWithVarianceBands
% Tests fcn_Path_plotPathXYWithVarianceBands
       
% Revision history:
% 2021_01_05
% -- first write of the code
% 2021_01_07
% -- fixed naming convention on functions to reflect change from path to
% traversal notation
% 2021_01_08
% -- more fixes
% 2025_07_01 - S. Brennan
% -- Removed traversal types, redid script/function based on
% plotTraversalXYWithVarianceBounds

close all


%% Test case 1: basic call for one trajectory
fig_num = 10001;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Fill in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;

% Call the function
fcn_Path_plotPathXYWithVarianceBands(paths{1}, [], fig_num);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Test case 2: advanced call for one trajectory - specify figure
fig_num = 10002;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Fill in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;
std_deviation = [];

% Call the function
fcn_Path_plotPathXYWithVarianceBands(paths{1},...
    std_deviation,fig_num);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Test case 3: advanced call for one trajectory - specify std_deviation
fig_num = 10003;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Fill in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;
std_deviation = 1;

% Call the function
fcn_Path_plotPathXYWithVarianceBands(paths{1},...
    std_deviation,fig_num);

title(sprintf('Standard deviation: %.0d meters',std_deviation));

fig_num = 32;
std_deviation = 2;
fcn_Path_plotPathXYWithVarianceBands(paths{1},...
    std_deviation,fig_num);
title(sprintf('Standard deviation: %.0d meters',std_deviation));

fig_num = 35;
std_deviation = 5;
fcn_Path_plotPathXYWithVarianceBands(paths{1},...
    std_deviation,fig_num);
title(sprintf('Standard deviation: %.0d meters',std_deviation));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Test case 4: advanced call for multiple trajectories
fig_num = 10004;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Fill in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;

std_deviation = 2;
for i_Path = 1:3
    fcn_Path_plotPathXYWithVarianceBands(paths{i_Path},...
        std_deviation,fig_num);
end

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
paths = fcn_Path_fillSamplePaths;

% Call the function
fcn_Path_plotPathXYWithVarianceBands(paths{1}, [], ([]));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

% Fill in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;

% Call the function
fcn_Path_plotPathXYWithVarianceBands(paths{1}, [], (-1));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

% Fill in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;

Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    fcn_Path_plotPathXYWithVarianceBands(paths{1}, [], ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    fcn_Path_plotPathXYWithVarianceBands(paths{1}, [], ([]));
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
