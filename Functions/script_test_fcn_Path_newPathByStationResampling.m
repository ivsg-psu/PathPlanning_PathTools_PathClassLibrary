% script_test_fcn_Path_newPathByStationResampling
% Tests the function: fcn_Path_newPathByStationResampling

% Revision history
%     2025_07_01: 
%     -- wrote the code originally in support of fcn_Path_findAveragePath

close all;


%% BASIC example: - start at zero
fig_num = 10001;
titleString = sprintf('BASIC example: start at zero');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Fill in sample paths (as a starter)
input_path = [0 0; 10 0; 20 0];


% Redecimate the path at 1-meter increments
interval = 5;
new_stations    = (0:interval:5)';

% Call the function
new_path = fcn_Path_newPathByStationResampling(input_path, new_stations, fig_num);

% Update title
title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(new_path));

% Check variable sizes
NreferencePoints = length(new_stations(:,1));
assert(isequal(size(new_path),[NreferencePoints 2]));

% Check variable values
assert(isequal(round(new_path(:,1),4),[0; 5]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% BASIC example: start at one
fig_num = 10002;
titleString = sprintf('BASIC example: start at one');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); 

% Fill in sample paths (as a starter)
input_path = [0 0; 10 0; 20 0];


% Redecimate the path at 1-meter increments
interval = 5;
new_stations    = (1:interval:6)';

% Call the function
new_path = fcn_Path_newPathByStationResampling(input_path, new_stations, fig_num);

% Update title
title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(new_path));

% Check variable sizes
NreferencePoints = length(new_stations(:,1));
assert(isequal(size(new_path),[NreferencePoints 2]));

% Check variable values
assert(isequal(round(new_path(:,1),4),[1; 6]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example: first point before start
fig_num = 10003;
titleString = sprintf('BASIC example: first point before start');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); 

% Fill in sample paths (as a starter)
input_path = [0 0; 10 0; 20 0];


% Redecimate the path
interval = 5;
new_stations    = (-5:interval:0)';

% Call the function
new_path = fcn_Path_newPathByStationResampling(input_path, new_stations, fig_num);

% Update title
title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(new_path));

% Check variable sizes
NreferencePoints = length(new_stations(:,1));
assert(isequal(size(new_path),[NreferencePoints 2]));
% Check variable values
assert(isequal(round(new_path(:,1),4),new_stations));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% BASIC example: perfect overlap
fig_num = 10004;
titleString = sprintf('BASIC example: perfect overlap');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num);

% Fill in sample paths (as a starter)
input_path = [0 0; 10 0; 20 0];


% Redecimate the path at 10-meter increments
interval = 10;
new_stations    = (0:interval:20)';

% Call the function
new_path = fcn_Path_newPathByStationResampling(input_path, new_stations, fig_num);

% Update title
title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(new_path));

% Check variable sizes
NreferencePoints = length(new_stations(:,1));
assert(isequal(size(new_path),[NreferencePoints 2]));

% Check variable values
assert(isequal(round(new_path(:,1),4),new_stations));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example: none on path
fig_num = 10005;
titleString = sprintf('BASIC example: none on path');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num);

% Fill in sample paths (as a starter)
input_path = [0 0; 10 0; 20 0];


% Redecimate the path at 1-meter increments
interval = 5;
new_stations    = (-10:interval:-5)';

% Call the function
new_path = fcn_Path_newPathByStationResampling(input_path, new_stations, fig_num);

% Update title
title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(new_path));

% Check variable sizes
NreferencePoints = length(new_stations(:,1));
assert(isequal(size(new_path),[NreferencePoints 2]));

% Check variable values
assert(isequal(round(new_path(:,1),4),new_stations));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% ADVANCED test: using a real path
fig_num = 20001;
titleString = sprintf('ADVANCED test: using a real path');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num);

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;
input_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1}, -1);

% Redecimate the path at 1-meter increments
interval = 10;
new_stations    = (0:interval:input_traversal.Station(end))';

% Call the function
new_path = fcn_Path_newPathByStationResampling(input_path, new_stations, fig_num);

% Update title
title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(new_path));

% Check variable sizes
NreferencePoints = length(new_stations(:,1));
assert(isequal(size(new_path),[NreferencePoints 2]));

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
input_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1}, -1);

% Redecimate the path at 10-meter increments
interval = 10;
new_stations    = (0:interval:input_traversal.Station(end))';

% Call the function
new_path = fcn_Path_newPathByStationResampling(input_path, new_stations, ([]));

% Check variable types
assert(isnumeric(new_path));

% Check variable sizes
NreferencePoints = length(new_stations(:,1));
assert(isequal(size(new_path),[NreferencePoints 2]));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);


% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;
input_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1}, -1);

% Redecimate the path at 10-meter increments
interval = 10;
new_stations    = (0:interval:input_traversal.Station(end))';

% Call the function
new_path = fcn_Path_newPathByStationResampling(input_path, new_stations, (-1));

% Check variable types
assert(isnumeric(new_path));

% Check variable sizes
NreferencePoints = length(new_stations(:,1));
assert(isequal(size(new_path),[NreferencePoints 2]));

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
input_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1}, -1);

% Redecimate the path at 10-meter increments
interval = 10;
new_stations    = (0:interval:input_traversal.Station(end))';

Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    new_path = fcn_Path_newPathByStationResampling(input_path, new_stations, ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    new_path = fcn_Path_newPathByStationResampling(input_path, new_stations, (-1));
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
