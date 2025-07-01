% script_test_fcn_Path_cleanPathFromForwardBackwardJogs
% Tests the function: fcn_Path_cleanPathFromForwardBackwardJogs

% Revision history
% 2020_01_09
% -- first write of the code
% -- need to add more test cases
% 2025_06_30
% -- testing more with real-world cases

close all;

%% Demonstration cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% _____                                 _             _   _                _____
% |  __ \                               | |           | | (_)              / ____|
% | |  | | ___ _ __ ___   ___  _ __  ___| |_ _ __ __ _| |_ _  ___  _ __   | |     __ _ ___  ___  ___
% | |  | |/ _ \ '_ ` _ \ / _ \| '_ \/ __| __| '__/ _` | __| |/ _ \| '_ \  | |    / _` / __|/ _ \/ __|
% | |__| |  __/ | | | | | (_) | | | \__ \ |_| | | (_| | |_| | (_) | | | | | |___| (_| \__ \  __/\__ \
% |_____/ \___|_| |_| |_|\___/|_| |_|___/\__|_|  \__,_|\__|_|\___/|_| |_|  \_____\__,_|___/\___||___/
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Demonstration%20Cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All demonstration case figures start with the number 1

%% Basic call example
fig_num = 10001;
titleString = sprintf('Basic call example, eliminating first point');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

path_with_jogs = [0 0; 1 1; 2 2; 3.3 3.36; 2.5 2.9; 3.5 3.6; 5 5];
clean_path = fcn_Path_cleanPathFromForwardBackwardJogs...
    (path_with_jogs,(fig_num));
title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(clean_path));

% Check variable sizes
assert(length(clean_path(:,1))<=length(path_with_jogs(:,1)));

% Check variable values
assert(isequal(clean_path,[0 0; 1 1; 2 2; 3.3 3.36; 3.5 3.6; 5 5]));

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

path_with_jogs = [0 0; 1 1; 2 2.2; 3.3 3; 2.5 2.9; 3.5 3.6; 5 5];
clean_path = fcn_Path_cleanPathFromForwardBackwardJogs...
    (path_with_jogs,([]));

% Check variable types
assert(isnumeric(clean_path));

% Check variable sizes
assert(length(clean_path(:,1))<=length(path_with_jogs(:,1)));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

path_with_jogs = [0 0; 1 1; 2 2.2; 3.3 3; 2.5 2.9; 3.5 3.6; 5 5];
clean_path = fcn_Path_cleanPathFromForwardBackwardJogs...
    (path_with_jogs,(-1));

% Check variable types
assert(isnumeric(clean_path));

% Check variable sizes
assert(length(clean_path(:,1))<=length(path_with_jogs(:,1)));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

path_with_jogs = [0 0; 1 1; 2 2.2; 3.3 3; 2.5 2.9; 3.5 3.6; 5 5];

Niterations = 20;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    clean_path = fcn_Path_cleanPathFromForwardBackwardJogs...
        (path_with_jogs,([])); %#ok<NASGU>
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    clean_path = fcn_Path_cleanPathFromForwardBackwardJogs...
        (path_with_jogs,(-1)); %#ok<NASGU>
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

close all;

%% BUG - cuts off path when should not (FIXED)
fig_num = 90001;
titleString = sprintf('BUG - cuts off path when should not (FIXED)');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

load('testData1_fcn_Path_cleanPathFromForwardBackwardJogs.mat','rawPath');
path_with_jogs = rawPath;
clean_path = fcn_Path_cleanPathFromForwardBackwardJogs...
    (path_with_jogs,(fig_num));
title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(clean_path));

% Check variable sizes
assert(length(clean_path(:,1))<=length(path_with_jogs(:,1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BUG - cuts off path when should not (FIXED)
fig_num = 90002;
titleString = sprintf('BUG - cuts off path when should not (FIXED)');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

load('testData2_fcn_Path_cleanPathFromForwardBackwardJogs.mat','rawPath');
path_with_jogs = rawPath;
clean_path = fcn_Path_cleanPathFromForwardBackwardJogs...
    (path_with_jogs,(fig_num));
title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(clean_path));

% Check variable sizes
assert(length(clean_path(:,1))<=length(path_with_jogs(:,1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% BUG - cuts off path weirdly
fig_num = 90003;
titleString = sprintf('Real world example where path cuts wierdly');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

load('testData3_fcn_Path_cleanPathFromForwardBackwardJogs.mat','rawPath');
path_with_jogs = rawPath;
clean_path = fcn_Path_cleanPathFromForwardBackwardJogs...
    (path_with_jogs,(fig_num));
title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(clean_path));

% Check variable sizes
assert(length(clean_path(:,1))<=length(path_with_jogs(:,1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


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
