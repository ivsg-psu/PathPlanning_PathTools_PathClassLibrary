% script_test_fcn_Path_findOrthogonalPathVectorsAtStations
% This is a script to exercise the function: fcn_Path_findOrthogonalPathVectorsAtStations.m.m
% This function was written on 2020_12_31 by S. Brennan
%     Modified on 2020_12_31 using script_test_fcn_Path_FindOrthogonalHitFromPathToPath
% Questions or comments? sbrennan@psu.edu

% Revision history
% 2020_12_31
% -- first write of the code
% 2021_01_09
% -- added more comments during clean-up
% 2022_01_03
% -- added assertion tests

close all;

%% BASIC example - simple horizontal line
fig_num = 10001;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

stations = 1; % Define the station
flag_rounding_type = 1; % Define the rounding type

% Create a dummy central path
central_path = [0 0; 4 0];  

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalPathVectorsAtStations(...
    stations,central_path,flag_rounding_type,fig_num); %#ok<*ASGLU>

% Make sure function worked
assert(isequal(round(unit_normal_vector_start,4),[ 1 0]));
assert(isequal(round(unit_normal_vector_end,4),[ 1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example - simple horizontal line with flag type 4
fig_num = 10002;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

stations = linspace(0,4,10)'; % Define the stations
flag_rounding_type = 4; % Define the rounding type

% Create a dummy central path
central_path = [0 0; 4 0];  


% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalPathVectorsAtStations(...
    stations,central_path,flag_rounding_type,fig_num); 

% Make sure function worked
% assert(isequal(round(unit_normal_vector_start,4),[ 1 0]));
% assert(isequal(round(unit_normal_vector_end,4),[ 1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% BASIC example - angled line segment - flag 1
fig_num = 10003;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

stations = 2;
flag_rounding_type = 1; % Define the rounding type
central_path = [0 0; 2 0; 2 -2];



% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalPathVectorsAtStations(...
    stations,central_path,flag_rounding_type,fig_num);

% Make sure function worked
assert(isequal(round(unit_normal_vector_start,4),[2 0]));
assert(isequal(round(unit_normal_vector_end,4),[ 2 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% BASIC example 2 - angled line segment - flag 2
fig_num = 10004;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

stations = 2;
flag_rounding_type = 2; % Define the rounding type
central_path = [0 0; 2 0; 2 -2];


% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalPathVectorsAtStations(...
    stations,central_path,flag_rounding_type,fig_num);

% Make sure function worked
assert(isequal(round(unit_normal_vector_start,4),[2 0]));
assert(isequal(round(unit_normal_vector_end,4),[ 3 0]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% BASIC example 2 - angled line segment - flag 3
fig_num = 10005;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

stations = 2;
flag_rounding_type = 3; % Define the rounding type
central_path = [0 0; 2 0; 2 -2];


% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalPathVectorsAtStations(...
    stations,central_path,flag_rounding_type,fig_num);

% Make sure function worked
assert(isequal(round(unit_normal_vector_start,4),[2 0]));
assert(isequal(round(unit_normal_vector_end,4),[2.7071 .7071]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% BASIC example 2 - angled line segment - flag 4
fig_num = 10005;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

stations = 2;
flag_rounding_type = 4; % Define the rounding type
central_path = [0 0; 2 0; 2 -2];


% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalPathVectorsAtStations(...
    stations,central_path,flag_rounding_type,fig_num);

% Make sure function worked
assert(isequal(round(unit_normal_vector_start,4),[2 0]));
assert(isequal(round(unit_normal_vector_end,4),[2.7071 .7071]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));



%% BASIC example - angled line segment adjacent to endpoint 
fig_num = 10006;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

stations = 2*2^0.5;
flag_rounding_type = 1; % Define the rounding type
central_path = [0 0; 2 2];


% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalPathVectorsAtStations(...
    stations,central_path,flag_rounding_type,fig_num);

% Make sure function worked
assert(isequal(round(unit_normal_vector_start,4),[     2     2]));
assert(isequal(round(unit_normal_vector_end,4),  [1.2929    2.7071]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% BASIC example - angled line segment adjacent to startpoint
fig_num = 10007;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

stations = 0;
flag_rounding_type = 1; % Define the rounding type
central_path = [0 0; 2 2];


% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalPathVectorsAtStations(...
    stations,central_path,flag_rounding_type,fig_num);

% Make sure function worked
assert(isequal(round(unit_normal_vector_start,4),[      0    0]));
assert(isequal(round(unit_normal_vector_end,4),  [-0.7071    0.7071]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% AVERAGING examples
fig_num = 10008;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Set up data
% close all
central_path = [0 0; 1 1; 2 0];

central_path_stations = fcn_Path_calcPathStation(central_path,-1);
stations = [linspace(0,central_path_stations(end),20)'; 2^0.5];

% AVERAGING example 1 - default setting
flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

subplot(1,4,1);
% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalPathVectorsAtStations(...
    stations,central_path,flag_rounding_type,fig_num);
title('Using prior segment (default, flag=1)');
axis equal
xlim([-1 3]);
ylim([-1 3]);

% AVERAGING example 2 - use following segment
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use averagae projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

%fig_num = 12;  % Define the figure
subplot(1,4,2);

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalPathVectorsAtStations(...
    stations,central_path,flag_rounding_type,fig_num);
title('Using following segment (flag=2)');
axis equal
xlim([-1 3]);
ylim([-1 3]);
legend off

% AVERAGING example 3 - use average of both segments
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

% fig_num = 13;  % Define the figure
subplot(1,4,3);

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalPathVectorsAtStations(...
    stations,central_path,flag_rounding_type,fig_num);
title('Averaging only at vertex (flag=3)');
axis equal
xlim([-1 3]);
ylim([-1 3]);
legend off

% AVERAGING example 4 - use average always
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

subplot(1,4,4);
% fig_num = 14;  % Define the figure

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalPathVectorsAtStations(...
    stations,central_path,flag_rounding_type,fig_num);
title('Continous averaging (flag=4)');
axis equal
xlim([-1 3]);
ylim([-1 3]);
legend off;

sgtitle('Comparison of projection types');


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));



%% NEGATIVE examples
fig_num = 10009;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Prep the example and workspace
% close all;
central_path = [-2 1; 1 4; 3 2];
central_path_stations = fcn_Path_calcPathStation(central_path,-1);
stations = [0; 1.5; 3; 3.5; 18^0.5-0.1; 18^0.5; 18^0.5+.1; 5; 5.5; 6.5; central_path_stations(end)];

% NEGATIVE example 1 - default setting
flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

% Calculate the unit normal vectors at given stations and put results into
% the figure.
subplot(2,2,1)
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalPathVectorsAtStations(...
    stations,central_path,flag_rounding_type,fig_num);
title('Vertex projection via prior segment (default, flag=1)');

% NEGATIVE example 2 - using following
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

% Calculate the unit normal vectors at given stations and put results into
% the figure.
subplot(2,2,2);
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalPathVectorsAtStations(...
    stations,central_path,flag_rounding_type,fig_num);
title('Vertex projection via following segment (flag=2)');

% NEGATIVE example 3 - using average at apex only
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

% Calculate the unit normal vectors at given stations and put results into
% the figure.
subplot(2,2,3);
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalPathVectorsAtStations(...
    stations,central_path,flag_rounding_type,fig_num);
title('Vertex projection via averaging prior and following segment at vertex (flag=3)');

% NEGATIVE example 4 - using average always
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

% Calculate the unit normal vectors at given stations and put results into
% the figure.
subplot(2,2,4)
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalPathVectorsAtStations(...
    stations,central_path,flag_rounding_type,fig_num);
title('Vertex projection via averaging everywhere (flag=4)');

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));



%% MULTICROSS examples
fig_num = 10010;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Setup
central_path = [-2 1; 1 4; 3 2; 5 2; 6 3; 7 2];

step_size = 0.2;
central_path_stations = fcn_Path_calcPathStation(central_path,-1);
stations = sort([(0:step_size:central_path_stations(end))'; central_path_stations]);
stations = unique(stations);

% MULTICROSS example 1 - default setting
flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

% Calculate the unit normal vectors at given stations and put results into
% the figure.
subplot(2,2,1);
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalPathVectorsAtStations(...
    stations,central_path,flag_rounding_type,fig_num);
title('Multicross example using projection via prior segment (default)');

% MULTICROSS example 2 - using following
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

% Calculate the unit normal vectors at given stations and put results into
% the figure.
subplot(2,2,2);
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalPathVectorsAtStations(...
    stations,central_path,flag_rounding_type,fig_num);
title('Multicross example using projection via following segment');

% MULTICROSS example 3 - using average at apex only
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

% Calculate the unit normal vectors at given stations and put results into
% the figure.
subplot(2,2,3);
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalPathVectorsAtStations(...
    stations,central_path,flag_rounding_type,fig_num);
title('Multicross example using projection via averaging of prior and following segment only at apex');


% MULTICROSS example 4 - using average always
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

% Calculate the unit normal vectors at given stations and put results into
% the figure.
subplot(2,2,4);
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalPathVectorsAtStations(...
    stations,central_path,flag_rounding_type,fig_num);
title('Multicross example using projection via averaging of prior and following segment always');

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%%

close all;

%% Real path examples
fig_num = 20001;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Fill in some dummy data
paths = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
clear data
for i_Path = 1:3
    traversal = fcn_Path_convertPathToTraversalStructure(paths{i_Path});
    data.traversal{i_Path} = traversal;
end


flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints

central_path = paths{1};
step_size = 10;
central_path_stations = fcn_Path_calcPathStation(central_path,-1);
stations = (0:step_size:central_path_stations(end))';
stations = unique(stations);

for i_Path = 1:3
    nearby_traversal  = data.traversal{i_Path};


    temp_fig_num = fig_num - 1 +i_Path;  % Define the figure
    
    % Calculate the unit normal vectors at given stations and put results into
    % the figure.
    [unit_normal_vector_start, unit_normal_vector_end] = ...
        fcn_Path_findOrthogonalPathVectorsAtStations(...
        stations,central_path,flag_rounding_type,temp_fig_num);
end

% % Make sure plot opened up
% assert(isequal(get(gcf,'Number'),fig_num));

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

stations = 1; % Define the station
flag_rounding_type = 1; % Define the rounding type

% Create a dummy central path
central_path = [0 0; 4 0];  


% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalPathVectorsAtStations(...
    stations,central_path,flag_rounding_type,[]); %#ok<*ASGLU>

% Make sure function worked
assert(isequal(round(unit_normal_vector_start,4),[ 1 0]));
assert(isequal(round(unit_normal_vector_end,4),[ 1 1]));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

stations = 1; % Define the station
flag_rounding_type = 1; % Define the rounding type

% Create a dummy central path
central_path = [0 0; 4 0];  


% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalPathVectorsAtStations(...
    stations,central_path,flag_rounding_type,-1);

% Make sure function worked
assert(isequal(round(unit_normal_vector_start,4),[ 1 0]));
assert(isequal(round(unit_normal_vector_end,4),[ 1 1]));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

stations = 1; % Define the station
flag_rounding_type = 1; % Define the rounding type

% Create a dummy central path
central_path = [0 0; 4 0];  


Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [unit_normal_vector_start, unit_normal_vector_end] = ...
        fcn_Path_findOrthogonalPathVectorsAtStations(...
        stations,central_path,flag_rounding_type,[]);
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [unit_normal_vector_start, unit_normal_vector_end] = ...
        fcn_Path_findOrthogonalPathVectorsAtStations(...
        stations,central_path,flag_rounding_type,-1);
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