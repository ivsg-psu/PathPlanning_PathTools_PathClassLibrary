% script_test_fcn_Path_findOrthogonalHitFromTraversalToTraversal.m
% This is a script to exercise the function: fcn_Path_findOrthogonalHitFromTraversalToTraversal.m
% This function was written on 2020_11_14 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history
%     2020_11_10
%     -- first write of the code
%     2020_11_14 
%     -- prep for Path class
%     2020_12_25 
%     -- include situation where central path and nearby path are the same
%     2021_01_07
%     -- lots of bug fixes as we demo for the team (lol)
%     2021_01_09
%     -- added more comments during clean-up
%     2022_01_03
%     -- added more checks for positive and neg case distances
%     -- output negative distances if in negative direction

close all;

%% BASIC example 1 - parallel lines, query is in middle area
stations = 1; % Define the station

% Create a dummy central path and convert it to a traversal
central_path = [0 0; 2 0];  
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);

% Define a "nearby" path and convert it to a traversal
nearby_path = [0 4; 2 4];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);

% Set default values
flag_rounding_type = 3;
search_radius = 5;
fig_num = 1;

% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,...
    central_traversal,nearby_traversal,...
    flag_rounding_type,search_radius,fig_num);

% Make sure function worked
assert(isequal(round(closest_path_point,4),[1     4]));
assert(isequal(round(distances,4),[4])); %#ok<*NBRAK>

print_results(stations,closest_path_point,distances);

%% BASIC example 1.5 - parallel lines, negative
stations = 1; % Define the station

% Create a dummy central path and convert it to a traversal
central_path = [0 0; 2 0];  
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);

% Define a "nearby" path and convert it to a traversal
nearby_path = [0 -4; 2 -4];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);

% Set default values
flag_rounding_type = 3;
search_radius = 5;
fig_num = 1;

% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,...
    central_traversal,nearby_traversal,...
    flag_rounding_type,search_radius,fig_num);

% Make sure function worked
assert(isequal(round(closest_path_point,4),[1     -4]));
assert(isequal(round(distances,4),[-4]));

print_results(stations,closest_path_point,distances);

%% BASIC example 2 - angled line segment adjacent to endpoint query
stations = 1;
central_path = [0 0; 2 0];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);

nearby_path = [0 4; 2 7];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);

% Set default values
flag_rounding_type = 3;
search_radius = 10;
fig_num = 2;

% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,...
    central_traversal,nearby_traversal,...
    flag_rounding_type,search_radius,fig_num);

% Make sure function worked
assert(isequal(round(closest_path_point,4),[1.0000    5.5000]));
assert(isequal(round(distances,4),[5.5000]));

print_results(stations,closest_path_point,distances);
%% BASIC example 3 - angled line segment adjacent to endpoint query 
stations = 10;
central_path = [0 0; 10 0];
central_traversal = ...
    fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [0 4; 10 7];
nearby_traversal = ...
    fcn_Path_convertPathToTraversalStructure(nearby_path);

% Set default values
flag_rounding_type = 3;
search_radius = 20;
fig_num = 3;

% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,...
    central_traversal,nearby_traversal,...
    flag_rounding_type,search_radius,fig_num);

% Make sure function worked
assert(isequal(round(closest_path_point,4),[10     7]));
assert(isequal(round(distances,4),[7]));

print_results(stations,closest_path_point,distances);
%% BASIC example 4 - angled line segment adjacent to startpoint query
stations = 0;
central_path = [0 0; 10 0];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);

nearby_path = [-1 4; 12 7];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);

% Set default values
flag_rounding_type = 3;
search_radius = 20;
fig_num = 4;

% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,...
    central_traversal,nearby_traversal,...
    flag_rounding_type,search_radius,fig_num);

% Make sure function worked
assert(isequal(round(closest_path_point,4),[0    4.2308]));
assert(isequal(round(distances,4),[4.2308]));

print_results(stations,closest_path_point,distances);

%% BASIC example 5 - parallel line segment adjacent to startpoint query
% Query point is right at start, so need to check it will not "miss"
stations = 0;
central_path = [0 0; 10 0];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);

nearby_path = [0 4; 10 4];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);

% Set default values
flag_rounding_type = 3;
search_radius = 20;
fig_num = 5;

% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,...
    central_traversal,nearby_traversal,...
    flag_rounding_type,search_radius,fig_num);

% Make sure function worked
assert(isequal(round(closest_path_point,4),[0     4]));
assert(isequal(round(distances,4),[4]));

print_results(stations,closest_path_point,distances);

%% BASIC example 6 - central path and nearby path are the same
% We should get that the very first point is the station point
stations = 1;
flag_rounding_type = 3;
search_radius = 10;
fig_num = 6;

central_path = [0 0; 2 2];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);

nearby_path = central_path;
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);

% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,...
    central_traversal,nearby_traversal,...
    flag_rounding_type,search_radius,fig_num);

% Make sure function worked
assert(isequal(round(closest_path_point,4),[0.7071    0.7071]));
assert(isequal(round(distances,4),[0]));

print_results(stations,closest_path_point,distances);


%% AVERAGING examples

% Set up data
central_path = [0 0; 1 1; 2 0];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [-1 0.5; 0.5 2; 1.5 2; 3 0.5];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);
stations = [0; 1; 2^0.5-0.1; 2^0.5; 2^0.5+.1; 2; central_traversal.Station(end)];
search_radius = 1.5; % Distance to search for nearby segments

% AVERAGING example 1 - default setting
fig_num = 11;
flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via prior segment (default, flag=1)');

% AVERAGING example 2 - use following segment
fig_num = 12;
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via following segment (flag=2)');

% AVERAGING example 3 - use average of both segments
fig_num = 13;
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via averaging prior and following segment at vertex (flag=3)');

% AVERAGING example 4 - use average always
fig_num = 14;
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via averaging everywhere (flag=4)');



%% NEGATIVE examples
fig_num = 20;

% Prep the example and workspace
central_path = [-2 1; 1 4; 3 2];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [-1 0.5; 0.5 2; 1.5 2; 3 0.5];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);
stations = [0; 1.5; 3; 3.5; 18^0.5-0.1; 18^0.5; 18^0.5+.1; 5; 5.5; 6.5; central_traversal.Station(end)];
search_radius = 1.5;

% NEGATIVE example 1 - default setting
flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via prior segment (default, flag=1)');

% NEGATIVE example 2 - using following
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via following segment (flag=2)');

% NEGATIVE example 3 - using average at apex only
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via averaging prior and following segment at vertex (flag=3)');

% NEGATIVE example 4 - using average always
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via averaging everywhere (flag=4)');

%% AVERAGING examples with search radius limitation

% Set up data

central_path = [0 0; 1 1; 2 0];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [-1 0.5; 0 2.5; 0.5 2; 1.5 2; 2 1; 3 3; 3 0.5];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);
stations = [0; 1; 2^0.5-0.1; 2^0.5; 2^0.5+.1; 2; central_traversal.Station(end)];
search_radius = 1.5; % Distance to search for nearby segments

% AVERAGING example 1 - default setting
fig_num = 1111;
flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via prior segment (default, flag=1), search radius limited to 1.5');

% AVERAGING example 2 - use following segment
fig_num = 2222;
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via following segment (flag=2), search radius limited to 1.5');

% AVERAGING example 3 - use average of both segments
fig_num = 3333;
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via averaging prior and following segment at vertex (flag=3), search radius limited to 1.5');

% AVERAGING example 4 - use average always
fig_num = 4444;
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via averaging everywhere (flag=4), search radius limited to 1.5');


%% NEGATIVE examples with search radius limitation


% Prep the example and workspace
central_path = [-2 1; 1 4; 3 2];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [-1 0.5; 0.5 2; 1.5 2; 3 0.5];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);
stations = [0; 1.5; 3; 3.5; 18^0.5-0.1; 18^0.5; 18^0.5+.1; 5; 5.5; 6.5; central_traversal.Station(end)];
search_radius = 1.5;

% NEGATIVE example 1 - default setting
fig_num = 11112;
flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via prior segment (default, flag=1)');

% NEGATIVE example 2 - using following
fig_num = 22223;
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via following segment (flag=2)');

% NEGATIVE example 3 - using average at apex only
fig_num = 33334;
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via averaging prior and following segment at vertex (flag=3)');

% NEGATIVE example 4 - using average always
fig_num = 44445;
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Vertex projection via averaging everywhere (flag=4)');

%% MULTICROSS examples
close all;
search_radius = 1.5;

% Setup
central_path = [-2 1; 1 4; 3 2; 5 2; 6 3; 7 2];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
nearby_path = [-1 0.5; 0.5 2; 1.5 2; 3 0.5; 4 3; 5 4; 6 3; 7 1];
nearby_traversal =  fcn_Path_convertPathToTraversalStructure(nearby_path);
step_size = 0.2;
stations = sort([[0:step_size:central_traversal.Station(end)]'; central_traversal.Station]);
stations = unique(stations);

% MULTICROSS example 1 - default setting
fig_num = 111;
flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Multicross example using projection via prior segment (default)');

% MULTICROSS example 2 - using following
fig_num = 222;
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Multicross example using projection via following segment');

% MULTICROSS example 3 - using average at apex only
fig_num = 333;
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Multicross example using projection via averaging of prior and following segment only at apex');


% MULTICROSS example 4 - using average always
fig_num = 444;
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type, search_radius, fig_num);
print_results(stations,closest_path_point,distances);
title('Multicross example using projection via averaging of prior and following segment always');


%% Real path examples
close all;

search_radius = 1.5;

% Fill in some dummy data
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
for i_Path = 1:length(paths_array)
    traversal = fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    data.traversal{i_Path} = traversal;
end

% Call the plot command to show results in XY
fig_num = 12;
fcn_Path_plotTraversalsXY(data,fig_num);

fig_num = 13;
flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints

central_traversal = data.traversal{1};
step_size = 10;
% stations = sort([[0:step_size:central_traversal.Station(end)]'; central_traversal.Station]);
stations = [0:step_size:central_traversal.Station(end)]';
stations = unique(stations);

for i_Path = 1:length(paths_array)
    nearby_traversal  = data.traversal{i_Path};

    [closest_path_point,distances] = ...
        fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,central_traversal,nearby_traversal,flag_rounding_type,30,fig_num);
end




function print_results(stations,closest_path_point,distances)
fprintf(1,'\n\nStation \t Location X \t Location Y \t Distance \n');
for i_station =1:length(stations)
    fprintf(1,'%.2f \t\t %.2f \t\t\t %.2f \t\t\t %.2f\n',stations(i_station),closest_path_point(i_station,1),closest_path_point(i_station,2),distances(i_station));
end
end
