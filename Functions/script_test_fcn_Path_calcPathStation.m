% script_test_fcn_Path_calcPathStation.m
% Tests fcn_Path_calcPathStation
       
% Revision history:
% 2025_06_29
% -- first write of the code, moving code set away from traversal(s)

close all

% FORMAT: 
%       [stations, differences] = fcn_Path_calcPathStation(path,(fig_num))

%% BASIC example: demonstration of station calculations for a 2D path
fig_num = 10001;
titleString = sprintf('BASIC example: demonstration of station calculations for a 2D path');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Fill in sample paths 
clear path
path = [0 0; 5 0; 5 3];
 
[stations, differences] = fcn_Path_calcPathStation(path,(fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(stations));
assert(isnumeric(differences));

% Check variable sizes
Nstations = length(path(:,1));
dimensionOfPath = length(path(1,:));
assert(isequal(size(stations),[Nstations 1]));
assert(isequal(size(differences),[Nstations dimensionOfPath]));

% Check variable values
assert(isequal(stations,[0; 5; 8]));
assert(isequal(differences,[0 0; 5 0; 0 3]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example: demonstration of station calculations for a path
fig_num = 10002;
titleString = sprintf('BASIC example: demonstration of station calculations for a 3D path');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;


% Fill in sample paths 
clear path
path = [0 0 0; 5 0 0; 5 3 0; 5 3 2];
 
[stations, differences] = fcn_Path_calcPathStation(path,(fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(stations));
assert(isnumeric(differences));

% Check variable sizes
Nstations = length(path(:,1));
dimensionOfPath = length(path(1,:));
assert(isequal(size(stations),[Nstations 1]));
assert(isequal(size(differences),[Nstations dimensionOfPath]));

% Check variable values
assert(isequal(stations,[0; 5; 8; 10]));
assert(isequal(differences,[0 0 0; 5 0 0; 0 3 0; 0 0 2]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example: demonstration of cell array of 2D paths
fig_num = 10003;
titleString = sprintf('BASIC example: demonstration of cell array of 2D paths');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;


% Fill in sample paths 
clear path
path{1,1} = [0 0; 5 0; 5 3];
path{2,1} = [-1 -1; 4 -1; 4 4; 4 5];
 
[stations, differences] = fcn_Path_calcPathStation(path,(fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(iscell(stations));
assert(iscell(differences));
assert(isnumeric(stations{1}));
assert(isnumeric(differences{1}));

% Check variable sizes
for ith_path = 1:length(path)
    Nstations = length(path{ith_path}(:,1));
    dimensionOfPath = length(path{ith_path}(1,:));
    assert(isequal(size(stations{ith_path}),[Nstations 1]));
    assert(isequal(size(differences{ith_path}),[Nstations dimensionOfPath]));
end

% Check variable values
assert(isequal(stations{1},[0; 5; 8]));
assert(isequal(differences{1},[0 0; 5 0; 0 3]));
assert(isequal(stations{2},[0; 5; 10; 11]));
assert(isequal(differences{2},[0 0; 5 0; 0 5; 0 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example: demonstration of cell array of 3D paths
fig_num = 10004;
titleString = sprintf('BASIC example: demonstration of cell array of 3D paths');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;


% Fill in sample paths 
clear path
path{1,1} = [0 0 0; 5 0 0 ; 5 3 0; 5 3 2];
path{2,1} = [-1 -1 0; 4 -1 0; 4 -1 2; 4 4 2; 4 5 2];
 
[stations, differences] = fcn_Path_calcPathStation(path,(fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(iscell(stations));
assert(iscell(differences));
assert(isnumeric(stations{1}));
assert(isnumeric(differences{1}));

% Check variable sizes
for ith_path = 1:length(path)
    Nstations = length(path{ith_path}(:,1));
    dimensionOfPath = length(path{ith_path}(1,:));
    assert(isequal(size(stations{ith_path}),[Nstations 1]));
    assert(isequal(size(differences{ith_path}),[Nstations dimensionOfPath]));
end

% Check variable values
assert(isequal(stations{1},[0; 5; 8; 10]));
assert(isequal(differences{1},[0 0 0; 5 0 0; 0 3 0; 0 0 2]));
assert(isequal(stations{2},[0; 5; 7; 12; 13]));
assert(isequal(differences{2},[0 0 0; 5 0 0; 0 0 2; 0 5 0; 0 1 0]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% PATH example: demonstration of sample path
fig_num = 20001;
titleString = sprintf('PATH example: demonstration of sample path');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;


% Fill in sample paths 
clear path
paths_array = fcn_Path_fillSamplePaths;
path = paths_array{1}; % Pick first path as reference_traversal structure

[stations, differences] = fcn_Path_calcPathStation(path,(fig_num));

% Check variable types
assert(isnumeric(stations));
assert(isnumeric(differences));

% Check variable sizes
Nstations = length(path(:,1));
dimensionOfPath = length(path(1,:));
assert(isequal(size(stations),[Nstations 1]));
assert(isequal(size(differences),[Nstations dimensionOfPath]));

% Check variable values
% assert(isequal(stations,[0; 5; 8]));
% assert(isequal(differences,[0 0; 5 0; 0 3]));

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

% Fill in sample paths 
clear path
path = [0 0; 5 0; 5 3];
 
[stations, differences] = fcn_Path_calcPathStation(path,([]));

% Check variable types
assert(isnumeric(stations));
assert(isnumeric(differences));

% Check variable sizes
Nstations = length(path(:,1));
dimensionOfPath = length(path(1,:));
assert(isequal(size(stations),[Nstations 1]));
assert(isequal(size(differences),[Nstations dimensionOfPath]));

% Check variable values
assert(isequal(stations,[0; 5; 8]));
assert(isequal(differences,[0 0; 5 0; 0 3]));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

% Fill in sample paths 
clear path
path = [0 0; 5 0; 5 3];
 
[stations, differences] = fcn_Path_calcPathStation(path,(-1));

% Check variable types
assert(isnumeric(stations));
assert(isnumeric(differences));

% Check variable sizes
Nstations = length(path(:,1));
dimensionOfPath = length(path(1,:));
assert(isequal(size(stations),[Nstations 1]));
assert(isequal(size(differences),[Nstations dimensionOfPath]));

% Check variable values
assert(isequal(stations,[0; 5; 8]));
assert(isequal(differences,[0 0; 5 0; 0 3]));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

% Fill in sample paths 
clear path
path = [0 0; 5 0; 5 3];
 
Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
[stations, differences] = fcn_Path_calcPathStation(path,([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
[stations, differences] = fcn_Path_calcPathStation(path,(-1));
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
