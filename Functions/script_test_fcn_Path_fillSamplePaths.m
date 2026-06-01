% script_test_fcn_Path_fillSamplePaths.m
% tests fcn_Path_fillSamplePaths.m

% Revision history
% 2020_11_10
% - first write of the code
% 2021_01_07
% - cleanup to include trajectories functions

close all

%% Test case 1: Basic call loading default paths
figNum = 10001;
fprintf(1,'Figure %.0f: Test case 1: Basic call loading default paths\n',figNum);
figure(figNum); clf;


path_number = [];
paths_array = fcn_Path_fillSamplePaths((path_number),(figNum));

% Check the variable type
assert(iscell(paths_array));

% Check the variable size
assert(isequal(size(paths_array),[1 4]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% Demo case 2: Show how to load specific paths
figNum = 10002;
fprintf(1,'Figure %.0f: Demo case 2: Show how to load specific paths\n',figNum);
figure(figNum); clf;

path_number = 2;
paths_array = fcn_Path_fillSamplePaths((path_number),(figNum));

% Check the variable type
assert(iscell(paths_array));

% Check the variable size
assert(isequal(size(paths_array),[1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));


%% Demo case 3: a test path for testing the conversions from XY to Sy
figNum = 10003;
fprintf(1,'Figure %.0f: Demo case 3: a test path for testing the conversions from XY to Sy\n',figNum);
figure(figNum); clf;

path_number = 4;
paths_array = fcn_Path_fillSamplePaths((path_number),(figNum));

% Check the variable type
assert(iscell(paths_array));

% Check the variable size
assert(isequal(size(paths_array),[1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));


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

path_number = 4;
paths_array = fcn_Path_fillSamplePaths((path_number),([]));

% Check the variable type
assert(iscell(paths_array));

% Check the variable size
assert(isequal(size(paths_array),[1 1]));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Basic fast mode - NO FIGURE, FAST MODE
figNum = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, figNum=-1\n',figNum);
figure(figNum); close(figNum);

path_number = 4;
paths_array = fcn_Path_fillSamplePaths((path_number),(-1));

% Check the variable type
assert(iscell(paths_array));

% Check the variable size
assert(isequal(size(paths_array),[1 1]));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
figNum = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',figNum);
figure(figNum);
close(figNum);


path_number = 4;

Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    paths_array = fcn_Path_fillSamplePaths((path_number),([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    paths_array = fcn_Path_fillSamplePaths((path_number),(-1));
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

