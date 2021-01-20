% script_test_fcn_Path_FindOrthogonalHitFromPathToPath.m
% This is a script to exercise the function: fcn_Path_FindOrthogonalHitFromPathToPath.m
% This function was written on 2020_11_14 by S. Brennan
%     Modified on 2020_11_14 to prep for Path class
%     Modified on 2020_12_25 to include situation where central path and nearby path are the same
% Questions or comments? sbrennan@psu.edu

close all;

%% BASIC example 1 - parallel lines, query is in middle area
stations = 1; % Define the station

% Create a dummy central path and convert it to a traversal
central_path = [0 0; 2 0];  
central_traversal = fcn_Path_convertXYtoTraversalStructure(central_path(:,1),central_path(:,2));

% Define a "nearby" path and convert it to a traversal
nearby_path = [0 4; 2 4];
nearby_traversal =  fcn_Path_convertXYtoTraversalStructure(nearby_path(:,1),nearby_path(:,2));

% Set default values
flag_rounding_type = 3;
search_radius = 20;
fig_num = 1;

% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_FindOrthogonalHitFromPathToPath(stations,...
    central_traversal,nearby_traversal,...
    flag_rounding_type,search_radius,fig_num);

print_results(stations,closest_path_point,distances);
%% BASIC example 2 - angled line segment adjacent to endpoint query
stations = 1;
central_path = [0 0; 2 0];
central_traversal = fcn_Path_convertXYtoTraversalStructure(central_path(:,1),central_path(:,2));

nearby_path = [0 4; 2 7];
nearby_traversal =  fcn_Path_convertXYtoTraversalStructure(nearby_path(:,1),nearby_path(:,2));

% Set default values
flag_rounding_type = 3;
search_radius = 20;
fig_num = 2;

% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_FindOrthogonalHitFromPathToPath(stations,...
    central_traversal,nearby_traversal,...
    flag_rounding_type,search_radius,fig_num);

print_results(stations,closest_path_point,distances);
%% BASIC example 3 - angled line segment adjacent to endpoint query 
stations = 10;
central_path = [0 0; 10 0];
central_traversal = ...
    fcn_Path_convertXYtoTraversalStructure(central_path(:,1),central_path(:,2));
nearby_path = [0 4; 10 7];
nearby_traversal = ...
    fcn_Path_convertXYtoTraversalStructure(nearby_path(:,1),nearby_path(:,2));

% Set default values
flag_rounding_type = 3;
search_radius = 20;
fig_num = 3;

% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_FindOrthogonalHitFromPathToPath(stations,...
    central_traversal,nearby_traversal,...
    flag_rounding_type,search_radius,fig_num);

print_results(stations,closest_path_point,distances);
%% BASIC example 4 - angled line segment adjacent to startpoint query
stations = 0;
central_path = [0 0; 10 0];
central_traversal = fcn_Path_convertXYtoTraversalStructure(central_path(:,1),central_path(:,2));

nearby_path = [-1 4; 12 7];
nearby_traversal =  fcn_Path_convertXYtoTraversalStructure(nearby_path(:,1),nearby_path(:,2));

% Set default values
flag_rounding_type = 3;
search_radius = 20;
fig_num = 4;

% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_FindOrthogonalHitFromPathToPath(stations,...
    central_traversal,nearby_traversal,...
    flag_rounding_type,search_radius,fig_num);

print_results(stations,closest_path_point,distances);

%% BASIC example 5 - parallel line segment adjacent to startpoint query
stations = 0;
central_path = [0 0; 10 0];
central_traversal = fcn_Path_convertXYtoTraversalStructure(central_path(:,1),central_path(:,2));

nearby_path = [0 4; 10 4];
nearby_traversal =  fcn_Path_convertXYtoTraversalStructure(nearby_path(:,1),nearby_path(:,2));

% Set default values
flag_rounding_type = 3;
search_radius = 20;
fig_num = 5;

% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_FindOrthogonalHitFromPathToPath(stations,...
    central_traversal,nearby_traversal,...
    flag_rounding_type,search_radius,fig_num);

print_results(stations,closest_path_point,distances);

%% BASIC example 6 - central path and nearby path are the same
stations = 1;
flag_rounding_type = 3;
search_radius = 10;
fig_num = 6;

central_path = [0 0; 2 2];
central_traversal = fcn_Path_convertXYtoTraversalStructure(central_path(:,1),central_path(:,2));

nearby_path = central_path;
nearby_traversal =  fcn_Path_convertXYtoTraversalStructure(nearby_path(:,1),nearby_path(:,2));

% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_FindOrthogonalHitFromPathToPath(stations,...
    central_traversal,nearby_traversal,...
    flag_rounding_type,search_radius,fig_num);

print_results(stations,closest_path_point,distances);


%% AVERAGING examples

% Set up data
close all
central_path = [0 0; 1 1; 2 0];
central_traversal = fcn_Path_convertXYtoTraversalStructure(central_path(:,1),central_path(:,2));
nearby_path = [-1 0.5; 0.5 2; 1.5 2; 3 0.5];
nearby_traversal =  fcn_Path_convertXYtoTraversalStructure(nearby_path(:,1),nearby_path(:,2));
stations = [0; 1; 2^0.5-0.1; 2^0.5; 2^0.5+.1; 2; central_traversal.Station(end)];

% AVERAGING example 1 - default setting
flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_FindOrthogonalHitFromPathToPath(stations,central_traversal,nearby_traversal,flag_rounding_type);
print_results(stations,closest_path_point,distances);
title('Vertex projection via prior segment (default, flag=1)');

% AVERAGING example 2 - use following segment
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_FindOrthogonalHitFromPathToPath(stations,central_traversal,nearby_traversal,flag_rounding_type);
print_results(stations,closest_path_point,distances);
title('Vertex projection via following segment (flag=2)');

% AVERAGING example 3 - use average of both segments
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_FindOrthogonalHitFromPathToPath(stations,central_traversal,nearby_traversal,flag_rounding_type);
print_results(stations,closest_path_point,distances);
title('Vertex projection via averaging prior and following segment at vertex (flag=3)');

% AVERAGING example 4 - use average always
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_FindOrthogonalHitFromPathToPath(stations,central_traversal,nearby_traversal,flag_rounding_type);
print_results(stations,closest_path_point,distances);
title('Vertex projection via averaging everywhere (flag=4)');



%% NEGATIVE examples

% Prep the example and workspace
close all;
central_path = [-2 1; 1 4; 3 2];
central_traversal = fcn_Path_convertXYtoTraversalStructure(central_path(:,1),central_path(:,2));
nearby_path = [-1 0.5; 0.5 2; 1.5 2; 3 0.5];
nearby_traversal =  fcn_Path_convertXYtoTraversalStructure(nearby_path(:,1),nearby_path(:,2));
stations = [0; 1.5; 3; 3.5; 18^0.5-0.1; 18^0.5; 18^0.5+.1; 5; 5.5; 6.5; central_traversal.Station(end)];

% NEGATIVE example 1 - default setting
flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_FindOrthogonalHitFromPathToPath(stations,central_traversal,nearby_traversal,flag_rounding_type);
print_results(stations,closest_path_point,distances);
title('Vertex projection via prior segment (default, flag=1)');

% NEGATIVE example 2 - using following
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_FindOrthogonalHitFromPathToPath(stations,central_traversal,nearby_traversal,flag_rounding_type);
print_results(stations,closest_path_point,distances);
title('Vertex projection via following segment (flag=2)');

% NEGATIVE example 3 - using average at apex only
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_FindOrthogonalHitFromPathToPath(stations,central_traversal,nearby_traversal,flag_rounding_type);
print_results(stations,closest_path_point,distances);
title('Vertex projection via averaging prior and following segment at vertex (flag=3)');

% NEGATIVE example 4 - using average always
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_FindOrthogonalHitFromPathToPath(stations,central_traversal,nearby_traversal,flag_rounding_type);
print_results(stations,closest_path_point,distances);
title('Vertex projection via averaging everywhere (flag=4)');


%% MULTICROSS examples
close all;

% Setup
central_path = [-2 1; 1 4; 3 2; 5 2; 6 3; 7 2];
central_traversal = fcn_Path_convertXYtoTraversalStructure(central_path(:,1),central_path(:,2));
nearby_path = [-1 0.5; 0.5 2; 1.5 2; 3 0.5; 4 3; 5 4; 6 3; 7 1];
nearby_traversal =  fcn_Path_convertXYtoTraversalStructure(nearby_path(:,1),nearby_path(:,2));
step_size = 0.2;
stations = sort([[0:step_size:central_traversal.Station(end)]'; central_traversal.Station]);
stations = unique(stations);

% MULTICROSS example 1 - default setting
flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_FindOrthogonalHitFromPathToPath(stations,central_traversal,nearby_traversal,flag_rounding_type);
print_results(stations,closest_path_point,distances);
title('Multicross example using projection via prior segment (default)');

% MULTICROSS example 2 - using following
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_FindOrthogonalHitFromPathToPath(stations,central_traversal,nearby_traversal,flag_rounding_type);
print_results(stations,closest_path_point,distances);
title('Multicross example using projection via following segment');

% MULTICROSS example 3 - using average at apex only
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_FindOrthogonalHitFromPathToPath(stations,central_traversal,nearby_traversal,flag_rounding_type);
print_results(stations,closest_path_point,distances);
title('Multicross example using projection via averaging of prior and following segment only at apex');


% MULTICROSS example 4 - using average always
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation
[closest_path_point,distances] = ...
    fcn_Path_FindOrthogonalHitFromPathToPath(stations,central_traversal,nearby_traversal,flag_rounding_type);
print_results(stations,closest_path_point,distances);
title('Multicross example using projection via averaging of prior and following segment always');


%% Real path examples
close all;

% Fill in some dummy data
paths = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
for i_Path = 1:length(paths)
    traversal = fcn_Path_convertXYtoTraversalStructure(paths{i_Path}(:,1),paths{i_Path}(:,2));
    data.traversal{i_Path} = traversal;
end

% Call the plot command to show results in XY
fig_num = 12;
fcn_Path_plotPathXY(data,fig_num);

fig_num = 13;
flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints

central_traversal = data.traversal{1};
step_size = 10;
% stations = sort([[0:step_size:central_traversal.Station(end)]'; central_traversal.Station]);
stations = [0:step_size:central_traversal.Station(end)]';
stations = unique(stations);

for i_Path = 1:length(paths)
    nearby_traversal  = data.traversal{i_Path};

    [closest_path_point,distances] = ...
        fcn_Path_FindOrthogonalHitFromPathToPath(stations,central_traversal,nearby_traversal,flag_rounding_type,30,fig_num);
end




function print_results(stations,closest_path_point,distances)
fprintf(1,'\n\nStation \t Location X \t Location Y \t Distance \n');
for i_station =1:length(stations)
    fprintf(1,'%.2f \t\t %.2f \t\t\t %.2f \t\t\t %.2f\n',stations(i_station),closest_path_point(i_station,1),closest_path_point(i_station,2),distances(i_station));
end
end
