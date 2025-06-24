% script_test_fcn_Path_addElevationToPath.m
% This is a script to exercise the function: fcn_Path_addElevationToPath.m
% This function was written on 2021_03_06 by Satya Prasad, szm888@psu.edu

% Revision history:
% 2021_03_20 - by S. Brennan
% - Added another example that is a bit more clearly "above" the query!
% 2024_03_14 - by S. Brennan
% - Fixed bug with 3D vectors not working right
% NOTE: the results do not seem to be correct!
 
close all;

%% BASIC example 1
fignum = 111;
figure(fignum); clf;

point = [0.5 0.2; 1.4 1.3];
reference_elevated_path = [0 0 0.1; 1 0 0.2; 2 0 0.3; 2 1 0.4];

elevated_path = fcn_Path_addElevationToPath(point, reference_elevated_path, fignum);

% Assert that the result is 3D
assert(length(elevated_path(1,:))==3);

%% BASIC example 1.2 - works
fignum = 112; % Define the figure number
figure(fignum); clf;

point = [0.5 0.2; 1.4 1.3]; % Define the query point as an XY
reference_elevated_path = [0 0 0.1; 0.5 0.2 0.2; 0.9 0.9 0.3; 1.5 0.6 0.4; 3 0 0.5]; % Define an XYZ path

% Snap the point onto the path
elevated_path = fcn_Path_addElevationToPath(point, reference_elevated_path, fignum); %#ok<*NASGU>

% Assert that the result is 3D
assert(length(elevated_path(1,:))==3);


%% BASIC example 1.3 - works
fignum = 113;
figure(fignum); clf;

point = [0.5 0.2; 1.4 1.3]; % Define the query point as an XY
reference_elevated_path = [0 0 0.1; 0.25 0.2 0.2; 0.9 0.9 0.3; 1.1 1.1 0.4; 2.3 2.7 0.5]; % Define an XYZ path

% Snap the point onto the path
elevated_path = fcn_Path_addElevationToPath(point, reference_elevated_path, fignum);

% Assert that the result is 3D
assert(length(elevated_path(1,:))==3);

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

point = [0.5 0.2; 1.4 1.3];
reference_elevated_path = [0 0 0.1; 1 0 0.2; 2 0 0.3; 2 1 0.4];
elevated_path = fcn_Path_addElevationToPath(point, reference_elevated_path, []);

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

point = [0.5 0.2; 1.4 1.3];
reference_elevated_path = [0 0 0.1; 1 0 0.2; 2 0 0.3; 2 1 0.4];
elevated_path = fcn_Path_addElevationToPath(point, reference_elevated_path, -1);

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

point = [0.5 0.2; 1.4 1.3];
reference_elevated_path = [0 0 0.1; 1 0 0.2; 2 0 0.3; 2 1 0.4];

Niterations = 100;

% Do calculation without pre-calculation
fig_num = [];
tic;
for ith_test = 1:Niterations
    % Call the function
    elevated_path = fcn_Path_addElevationToPath(point, reference_elevated_path, []);
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
fig_num = -1;
tic;
for ith_test = 1:Niterations
    % Call the function
    elevated_path = fcn_Path_addElevationToPath(point, reference_elevated_path, []);
end
fast_method = toc;

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

% Plot results as bar chart
figure(373737);
clf;
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