% script_test_fcn_Path_findOrthogonalHitFromPathToPath.m
% This is a script to exercise the function: fcn_Path_findOrthogonalHitFromPathToPath.m
% This function was written on 2020_11_14 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history
% 2020_11_10
% -- first write of the code
% 2020_11_14
% -- prep for Path class
% 2020_12_25
% -- include situation where central path and nearby path are the same
% 2021_01_07
% -- lots of bug fixes as we demo for the team (lol)
% 2021_01_09
% -- added more comments during clean-up
% 2022_01_03
% -- added more checks for positive and neg case distances
% -- output negative distances if in negative direction

close all;

%% BASIC example: Parallel lines, query is in middle area
fig_num = 10001;
fprintf(1,'Figure %.0f: BASIC example: Parallel lines, query is in middle area\n',fig_num);
figure(fig_num); clf;

stations = 1; % Define the station

% Create a dummy central path and convert it to a traversal
central_path = [0 0; 2 0];  


% Define a "nearby" path and convert it to a traversal
nearby_path = [0 4; 2 4];


% Set default values
flag_rounding_type = 3;
search_radius = 5;


% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,...
    central_path,nearby_path,...
    flag_rounding_type,search_radius,fig_num);

% title(titleString, 'Interpreter','none');

print_results(stations,closest_path_point,distances);

% Check variable types
assert(isnumeric(closest_path_point));
assert(isnumeric(distances));

% Check variable sizes
assert(isequal(size(closest_path_point),[1 2]));
assert(isequal(size(distances),[1 1]));

% Check variable values
assert(isequal(round(closest_path_point,4),[1     4]));
assert(isequal(round(distances,4),[4])); %#ok<*NBRAK>

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example: parallel lines, negative
fig_num = 10002;
fprintf(1,'Figure %.0f: BASIC example: parallel lines, negative\n',fig_num);
figure(fig_num); clf;

stations = 1; % Define the station

% Create a dummy central path and convert it to a traversal
central_path = [0 0; 2 0];  


% Define a "nearby" path and convert it to a traversal
nearby_path = [0 -4; 2 -4];


% Set default values
flag_rounding_type = 3;
search_radius = 5;

% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,...
    central_path,nearby_path,...
    flag_rounding_type,search_radius,fig_num);

% title(titleString, 'Interpreter','none');

print_results(stations,closest_path_point,distances);

% Check variable types
assert(isnumeric(closest_path_point));
assert(isnumeric(distances));

% Check variable sizes
assert(isequal(size(closest_path_point),[1 2]));
assert(isequal(size(distances),[1 1]));

% Check variable values
assert(isequal(round(closest_path_point,4),[1     -4]));
assert(isequal(round(distances,4),[-4]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% BASIC example: angled line segment adjacent to endpoint query
fig_num = 10003;
fprintf(1,'Figure %.0f: BASIC example: angled line segment adjacent to endpoint query\n',fig_num);
figure(fig_num); clf;

stations = 1;
central_path = [0 0; 2 0];


nearby_path = [0 4; 2 7];


% Set default values
flag_rounding_type = 3;
search_radius = 10;

% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,...
    central_path,nearby_path,...
    flag_rounding_type,search_radius,fig_num);

% title(titleString, 'Interpreter','none');

print_results(stations,closest_path_point,distances);

% Check variable types
assert(isnumeric(closest_path_point));
assert(isnumeric(distances));

% Check variable sizes
assert(isequal(size(closest_path_point),[1 2]));
assert(isequal(size(distances),[1 1]));

% Check variable values
assert(isequal(round(closest_path_point,4),[1.0000    5.5000]));
assert(isequal(round(distances,4),[5.5000]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% BASIC example: angled line segment adjacent to endpoint query 
fig_num = 10004;
fprintf(1,'Figure %.0f: BASIC example: angled line segment adjacent to endpoint query \n',fig_num);
figure(fig_num); clf;

stations = 10;
central_path = [0 0; 10 0];
nearby_path = [0 4; 10 7];

% Set default values
flag_rounding_type = 3;
search_radius = 20;

% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,...
    central_path,nearby_path,...
    flag_rounding_type,search_radius,fig_num);

% title(titleString, 'Interpreter','none');

print_results(stations,closest_path_point,distances);

% Check variable types
assert(isnumeric(closest_path_point));
assert(isnumeric(distances));

% Check variable sizes
assert(isequal(size(closest_path_point),[1 2]));
assert(isequal(size(distances),[1 1]));

% Check variable values
assert(isequal(round(closest_path_point,4),[10     7]));
assert(isequal(round(distances,4),[7]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% BASIC example: Angled line segment adjacent to startpoint query
fig_num = 10005;
fprintf(1,'Figure %.0f: BASIC example: Angled line segment adjacent to startpoint query\n',fig_num);
figure(fig_num); clf;

stations = 0;
central_path = [0 0; 10 0];
nearby_path = [-1 4; 12 7];


% Set default values
flag_rounding_type = 3;
search_radius = 20;

% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,...
    central_path,nearby_path,...
    flag_rounding_type,search_radius,fig_num);

% title(titleString, 'Interpreter','none');

print_results(stations,closest_path_point,distances);

% Check variable types
assert(isnumeric(closest_path_point));
assert(isnumeric(distances));

% Check variable sizes
assert(isequal(size(closest_path_point),[1 2]));
assert(isequal(size(distances),[1 1]));

% Check variable values
assert(isequal(round(closest_path_point,4),[0    4.2308]));
assert(isequal(round(distances,4),[4.2308]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% BASIC example: parallel line segment adjacent to startpoint query
% Query point is right at start, so need to check it will not "miss"
fig_num = 10006;
fprintf(1,'Figure %.0f: BASIC example: parallel line segment adjacent to startpoint query\n',fig_num);
figure(fig_num); clf;

stations = 0;
central_path = [0 0; 10 0];


nearby_path = [0 4; 10 4];


% Set default values
flag_rounding_type = 3;
search_radius = 20;

% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,...
    central_path,nearby_path,...
    flag_rounding_type,search_radius,fig_num);

% title(titleString, 'Interpreter','none');

print_results(stations,closest_path_point,distances);

% Check variable types
assert(isnumeric(closest_path_point));
assert(isnumeric(distances));

% Check variable sizes
assert(isequal(size(closest_path_point),[1 2]));
assert(isequal(size(distances),[1 1]));

% Check variable values
assert(isequal(round(closest_path_point,4),[0     4]));
assert(isequal(round(distances,4),[4]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% BASIC example: central path and nearby path are the same
% We should get that the very first point is the station point
fig_num = 10007;
fprintf(1,'Figure %.0f: BASIC example: central path and nearby path are the same\n',fig_num);
figure(fig_num); clf;

stations = 1;
flag_rounding_type = 3;
search_radius = 10;

central_path = [0 0; 2 2];


nearby_path = central_path;


% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,...
    central_path,nearby_path,...
    flag_rounding_type,search_radius,fig_num);

% title(titleString, 'Interpreter','none');

print_results(stations,closest_path_point,distances);

% Check variable types
assert(isnumeric(closest_path_point));
assert(isnumeric(distances));

% Check variable sizes
assert(isequal(size(closest_path_point),[1 2]));
assert(isequal(size(distances),[1 1]));

% Check variable values
assert(isequal(round(closest_path_point,4),[0.7071    0.7071]));
assert(isequal(round(distances,4),[0]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% BASIC example: central path and query point do not hit anything
% We should get NaN values

fig_num = 10008;
fprintf(1,'Figure %.0f: BASIC example: central path and query point do not hit anything\n',fig_num);
figure(fig_num); clf;

stations = 1;
flag_rounding_type = 3;
search_radius = 10;

central_path = [0 0; 2 0];


nearby_path = [5 5; 8 5];


% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,...
    central_path,nearby_path,...
    flag_rounding_type,search_radius,fig_num);

% title(titleString, 'Interpreter','none');

print_results(stations,closest_path_point,distances);

% Check variable types
assert(all(isnan(closest_path_point)));
assert(all(isnan(distances)));

% Check variable sizes
assert(isequal(size(closest_path_point),[1 2]));
assert(isequal(size(distances),[1 1]));

% Check variable values
% assert(isequal(round(closest_path_point,4),[0.7071    0.7071]));
% assert(isequal(round(distances,4),[0]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));



%% 
close all

%% Showing flag_rounding_type
fig_num = 20001;
fprintf(1,'Figure %.0f: Showing flag_rounding_type\n',fig_num);
figure(fig_num); clf;

% Set up data
central_path = [0 0; 1 1; 2 0];

nearby_path = [-1 0.5; 0.5 2; 1.5 2; 3 0.5];
[central_path_Station, ~] = fcn_Path_calcPathStation(central_path,-1);
stations = [0; 1; 2^0.5-0.1; 2^0.5; 2^0.5+.1; 2; central_path_Station(end)];
search_radius = 1.5; % Distance to search for nearby segments

% AVERAGING example 1 - default setting
subplot(2,2,1);
flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,central_path,nearby_path,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via prior segment (default, flag=1)');

% AVERAGING example 2 - use following segment
subplot(2,2,2);
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,central_path,nearby_path,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via following segment (flag=2)');

% AVERAGING example 3 - use average of both segments
subplot(2,2,3)
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,central_path,nearby_path,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via averaging prior and following segment at vertex (flag=3)');

% AVERAGING example 4 - use average always
subplot(2,2,4)
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,central_path,nearby_path,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via interpolation (flag=4)');

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% NEGATIVE examples
fig_num = 20002;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Prep the example and workspace
central_path = [-2 1; 1 4; 3 2];

nearby_path = [-1 0.5; 0.5 2; 1.5 2; 3 0.5];
[central_path_Station, ~] = fcn_Path_calcPathStation(central_path,-1);
stations = [0; 1.5; 3; 3.5; 18^0.5-0.1; 18^0.5; 18^0.5+.1; 5; 5.5; 6.5; central_path_Station(end)];
search_radius = 1.5;

% NEGATIVE example 1 - default setting
subplot(2,2,1);
flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,central_path,nearby_path,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via prior segment (default, flag=1)');

% NEGATIVE example 2 - using following
subplot(2,2,2);
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,central_path,nearby_path,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via following segment (flag=2)');

% NEGATIVE example 3 - using average at apex only
subplot(2,2,3);
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,central_path,nearby_path,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via averaging prior and following segment at vertex (flag=3)');

% NEGATIVE example 4 - using average always
subplot(2,2,4);
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,central_path,nearby_path,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via averaging everywhere (flag=4)');

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% AVERAGING examples with search radius limitation
fig_num = 20003;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Set up data
central_path = [0 0; 1 1; 2 0];

nearby_path = [-1 0.5; 0 2.5; 0.5 2; 1.5 2; 2 1; 3 3; 3 0.5];
[central_path_Station, ~] = fcn_Path_calcPathStation(central_path,-1);
stations = [0; 1; 2^0.5-0.1; 2^0.5; 2^0.5+.1; 2; central_path_Station(end)];
search_radius = 1.5; % Distance to search for nearby segments

% AVERAGING example 1 - default setting
subplot(2,2,1);
flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,central_path,nearby_path,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via prior segment (default, flag=1), search radius limited to 1.5');

% AVERAGING example 2 - use following segment
subplot(2,2,2);
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,central_path,nearby_path,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via following segment (flag=2), search radius limited to 1.5');

% AVERAGING example 3 - use average of both segments
subplot(2,2,3);
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,central_path,nearby_path,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via averaging prior and following segment at vertex (flag=3), search radius limited to 1.5');

% AVERAGING example 4 - use average always
subplot(2,2,4);
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,central_path,nearby_path,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via averaging everywhere (flag=4), search radius limited to 1.5');

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% NEGATIVE examples with search radius limitation
fig_num = 20004;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Prep the example and workspace
central_path = [-2 1; 1 4; 3 2];

nearby_path = [-1 0.5; 0.5 2; 1.5 2; 3 0.5];
[central_path_Station, ~] = fcn_Path_calcPathStation(central_path,-1);
stations = [0; 1.5; 3; 3.5; 18^0.5-0.1; 18^0.5; 18^0.5+.1; 5; 5.5; 6.5; central_path_Station(end)];
search_radius = 1.5;

% NEGATIVE example 1 - default setting
subplot(2,2,1);
flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,central_path,nearby_path,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via prior segment (default, flag=1)');

% NEGATIVE example 2 - using following
subplot(2,2,2);
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,central_path,nearby_path,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via following segment (flag=2)');

% NEGATIVE example 3 - using average at apex only
subplot(2,2,3);
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,central_path,nearby_path,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via averaging prior and following segment at vertex (flag=3)');

% NEGATIVE example 4 - using average always
subplot(2,2,4);
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,central_path,nearby_path,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via averaging everywhere (flag=4)');

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%%
close all;

%% MULTICROSS examples
fig_num = 30001;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

search_radius = 1.5;

% Setup
central_path = [-2 1; 1 4; 3 2; 5 2; 6 3; 7 2];

nearby_path = [-1 0.5; 0.5 2; 1.5 2; 3 0.5; 4 3; 5 4; 6 3; 7 1];
[central_path_Station, ~] = fcn_Path_calcPathStation(central_path,-1);
step_size = 0.2;
stations = sort([[0:step_size:central_path_Station(end)]'; central_path_Station]);
stations = unique(stations);

% MULTICROSS example 1 - default setting
subplot(2,2,1);
flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,central_path,nearby_path,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Multicross example using projection via prior segment (default)');

% MULTICROSS example 2 - using following
subplot(2,2,2);
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,central_path,nearby_path,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Multicross example using projection via following segment');

% MULTICROSS example 3 - using average at apex only
subplot(2,2,3);
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,central_path,nearby_path,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Multicross example using projection via averaging of prior and following segment only at apex');


% MULTICROSS example 4 - using average always
subplot(2,2,4);
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,central_path,nearby_path,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Multicross example using projection via averaging of prior and following segment always');

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%%

close all;

%% Real path examples
fig_num = 40001;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Fill in some dummy data
paths_array = fcn_Path_fillSamplePaths;

% Call the plot command to show results in XY
flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints

central_path = paths_array{1};
step_size = 10;
[central_path_Station, ~] = fcn_Path_calcPathStation(central_path,-1);
stations = [0:step_size:central_path_Station(end)]';
stations = unique(stations);

for i_Path = 1:length(paths_array)
    nearby_path  = paths_array{i_Path};

    [closest_path_point,distances] = ...
        fcn_Path_findOrthogonalHitFromPathToPath(stations,central_path,nearby_path,flag_rounding_type,30,fig_num); %#ok<ASGLU>
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


stations = 1;
central_path = [0 0; 2 0];


nearby_path = [0 4; 2 7];


% Set default values
flag_rounding_type = 3;
search_radius = 10;

% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,...
    central_path,nearby_path,...
    flag_rounding_type,search_radius,[]);

% Check variable types
assert(isnumeric(closest_path_point));
assert(isnumeric(distances));

% Check variable sizes
assert(isequal(size(closest_path_point),[1 2]));
assert(isequal(size(distances),[1 1]));

% Check variable values
assert(isequal(round(closest_path_point,4),[1.0000    5.5000]));
assert(isequal(round(distances,4),[5.5000]));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);


stations = 1;
central_path = [0 0; 2 0];
nearby_path = [0 4; 2 7];


% Set default values
flag_rounding_type = 3;
search_radius = 10;

% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromPathToPath(stations,...
    central_path,nearby_path,...
    flag_rounding_type,search_radius,-1);

% Check variable types
assert(isnumeric(closest_path_point));
assert(isnumeric(distances));

% Check variable sizes
assert(isequal(size(closest_path_point),[1 2]));
assert(isequal(size(distances),[1 1]));

% Check variable values
assert(isequal(round(closest_path_point,4),[1.0000    5.5000]));
assert(isequal(round(distances,4),[5.5000]));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);


stations = 1;
central_path = [0 0; 2 0];


nearby_path = [0 4; 2 7];


% Set default values
flag_rounding_type = 3;
search_radius = 10;

Niterations = 20;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [closest_path_point,distances] = ...
        fcn_Path_findOrthogonalHitFromPathToPath(stations,...
        central_path,nearby_path,...
        flag_rounding_type,search_radius,[]);
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [closest_path_point,distances] = ...
        fcn_Path_findOrthogonalHitFromPathToPath(stations,...
        central_path,nearby_path,...
        flag_rounding_type,search_radius,-1);
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

%%
function print_results(stations,closest_path_point,distances)
if 1==0
    fprintf(1,'\n\nStation \t Location X \t Location Y \t Distance \n');
    for i_station =1:length(stations)
        fprintf(1,'%.2f \t\t %.2f \t\t\t %.2f \t\t\t %.2f\n',stations(i_station),closest_path_point(i_station,1),closest_path_point(i_station,2),distances(i_station));
    end
end
end