% script_test_fcn_Path_convertPathToTraversalStructure.m
% tests fcn_Path_convertPathToTraversalStructure.m

% Revision history:
% 2020_11_12
% -- first wrote the code
% 2021_01_06
% -- a bit more comments, added more plotting, renamed
% 2021_01_07
% -- changed name to reflect that the input is a Path, and that we are
% using an array of paths
% -- more descriptive comments
% 2021_01_09
% -- added another example to illustrate station errors
% 2021_03_21
% -- added examples to illustrate 2D and 3D basic plotting

close all

%% Example 1.1 - show how it works for 2D
% Basic call 
fig_num = 10001;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

path_simple_2D = [1 1; 1 2; 3 4; 4 5; 7 7];
traversal = fcn_Path_convertPathToTraversalStructure(path_simple_2D, (fig_num));

% Check variable types
assert(isstruct(traversal));
assert(isfield(traversal,'X'))
assert(isfield(traversal,'Y'))
assert(isfield(traversal,'Z'))
assert(isfield(traversal,'Diff'))
assert(isfield(traversal,'Station'))
assert(isfield(traversal,'Yaw'))

% Check variable sizes
Npoints = length(path_simple_2D(:,1));
assert(length(traversal.X(:,1))==Npoints)
assert(length(traversal.Y(:,1))==Npoints)
assert(length(traversal.Z(:,1))==Npoints)
assert(length(traversal.Diff(:,1))==Npoints)
assert(length(traversal.Station(:,1))==Npoints)
assert(length(traversal.Yaw(:,1))==Npoints-1)

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Example 1.2 - show how it works for 3D
% Basic call 
fig_num = 10002;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

path_simple_3D = [1 1 1; 1 2 2; 3 4 3; 4 5 3; 7 7 4];
traversal = fcn_Path_convertPathToTraversalStructure(path_simple_3D,fig_num);

% Check variable types
assert(isstruct(traversal));
assert(isfield(traversal,'X'))
assert(isfield(traversal,'Y'))
assert(isfield(traversal,'Z'))
assert(isfield(traversal,'Diff'))
assert(isfield(traversal,'Station'))
assert(isfield(traversal,'Yaw'))

% Check variable sizes
Npoints = length(path_simple_3D(:,1));
assert(length(traversal.X(:,1))==Npoints)
assert(length(traversal.Y(:,1))==Npoints)
assert(length(traversal.Z(:,1))==Npoints)
assert(length(traversal.Diff(:,1))==Npoints)
assert(length(traversal.Station(:,1))==Npoints)
assert(length(traversal.Yaw(:,1))==Npoints-1)

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Example 1.3 - show how it works
fig_num = 10003;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;
figure(fig_num*10); clf;

% Fill in some dummy data
paths_array = fcn_Path_fillSamplePaths;
paths_array = {paths_array{1}, paths_array{2}, paths_array{3}};

% Basic call with a figure option to plot output repeatedly onto figure
clear data
for i_traveral = 1:length(paths_array)
    traversal = fcn_Path_convertPathToTraversalStructure(paths_array{i_traveral},fig_num);
    data.traversal{i_traveral} = traversal;
end

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

% Call the plot command to show ALL results in XY, which gives same result
fcn_Path_plotTraversalsXY(data,fig_num*10);
xlabel('X [m]');
ylabel('Y [m]');


%% Example 2: Show station discrepancies
fig_num = 10004;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;
figure(fig_num*10); clf;

paths_array = fcn_Path_fillSamplePaths;
paths_array = {paths_array{1}, paths_array{2}, paths_array{3}};

% Basic call with a figure option to plot output repeatedly onto figure
clear data
for i_traveral = 1:length(paths_array)
    traversal = fcn_Path_convertPathToTraversalStructure(paths_array{i_traveral},fig_num);
    data.traversal{i_traveral} = traversal;
end

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

figure(fig_num*10);
clf;
hold on;
grid on;
grid minor;

Station_step = 40;

fcn_Path_plotTraversalsXY(data,fig_num*10);
xlabel('X [m]');
ylabel('Y [m]');

for i_traveral = 1:length(data.traversal)
    traversal_stations = data.traversal{i_traveral}.Station;
    for i_station = Station_step:Station_step:traversal_stations(end)
        index = find(traversal_stations >= i_station,1);
        plot(data.traversal{i_traveral}.X(index),...
            data.traversal{i_traveral}.Y(index),...
            'k.','Markersize',15);
        text(data.traversal{i_traveral}.X(index),...
            data.traversal{i_traveral}.Y(index),...
            sprintf('Station: %.2f',...
            data.traversal{i_traveral}.Station(index)));
    end
end


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

% Fill in some dummy data
paths_array = fcn_Path_fillSamplePaths;
traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1},[]);

% Check variable types
assert(isstruct(traversal));
assert(isfield(traversal,'X'))
assert(isfield(traversal,'Y'))
assert(isfield(traversal,'Z'))
assert(isfield(traversal,'Diff'))
assert(isfield(traversal,'Station'))
assert(isfield(traversal,'Yaw'))

% Check variable sizes
Npoints = length(paths_array{1}(:,1));
assert(length(traversal.X(:,1))==Npoints)
assert(length(traversal.Y(:,1))==Npoints)
assert(length(traversal.Z(:,1))==Npoints)
assert(length(traversal.Diff(:,1))==Npoints)
assert(length(traversal.Station(:,1))==Npoints)
assert(length(traversal.Yaw(:,1))==Npoints-1)

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

% Fill in some dummy data
paths_array = fcn_Path_fillSamplePaths;
traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1},-1);

% Check variable types
assert(isstruct(traversal));
assert(isfield(traversal,'X'))
assert(isfield(traversal,'Y'))
assert(isfield(traversal,'Z'))
assert(isfield(traversal,'Diff'))
assert(isfield(traversal,'Station'))
assert(isfield(traversal,'Yaw'))

% Check variable sizes
Npoints = length(paths_array{1}(:,1));
assert(length(traversal.X(:,1))==Npoints)
assert(length(traversal.Y(:,1))==Npoints)
assert(length(traversal.Z(:,1))==Npoints)
assert(length(traversal.Diff(:,1))==Npoints)
assert(length(traversal.Station(:,1))==Npoints)
assert(length(traversal.Yaw(:,1))==Npoints-1)

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

paths_array = fcn_Path_fillSamplePaths;

Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1},[]);
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1},-1);
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

