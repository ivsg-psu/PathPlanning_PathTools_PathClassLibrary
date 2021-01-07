% script_test_fcn_Path_FindOrthogonalPathVectorsAtStations.m.m
% This is a script to exercise the function: fcn_Path_FindOrthogonalPathVectorsAtStations.m.m
% This function was written on 2020_12_31 by S. Brennan
%     Modified on 2020_12_31 using script_test_fcn_Path_FindOrthogonalHitFromPathToPath
% Questions or comments? sbrennan@psu.edu

close all;

%% BASIC example 1 - simple horizontal line
stations = 1; % Define the station
flag_rounding_type = 1; % Define the rounding type
fig_num = 1;  % Define the figure

% Create a dummy central path and convert it to a traversal
central_path = [0 0; 2 0];  
central_traversal = fcn_Path_convertXYtoTraversalStructure(central_path(:,1),central_path(:,2));

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_FindOrthogonalPathVectorsAtStations(...
    stations,central_traversal,flag_rounding_type,fig_num); %#ok<*ASGLU>

%% BASIC example 2 - angled line segment 
stations = 1;
central_path = [0 0; 2 2];
central_traversal = fcn_Path_convertXYtoTraversalStructure(central_path(:,1),central_path(:,2));

fig_num = 2;  % Define the figure

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_FindOrthogonalPathVectorsAtStations(...
    stations,central_traversal,flag_rounding_type,fig_num);

%% BASIC example 3 - angled line segment adjacent to endpoint 
stations = 2*2^0.5;
central_path = [0 0; 2 2];
central_traversal = fcn_Path_convertXYtoTraversalStructure(central_path(:,1),central_path(:,2));

fig_num = 3;  % Define the figure

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_FindOrthogonalPathVectorsAtStations(...
    stations,central_traversal,flag_rounding_type,fig_num);

%% BASIC example 4 - angled line segment adjacent to startpoint
stations = 0;
central_path = [0 0; 2 2];
central_traversal = fcn_Path_convertXYtoTraversalStructure(central_path(:,1),central_path(:,2));

fig_num = 4;  % Define the figure

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_FindOrthogonalPathVectorsAtStations(...
    stations,central_traversal,flag_rounding_type,fig_num);


%% AVERAGING examples

% Set up data
close all
central_path = [0 0; 1 1; 2 0];
central_traversal = fcn_Path_convertXYtoTraversalStructure(central_path(:,1),central_path(:,2));
stations = [0; 1; 2^0.5-0.1; 2^0.5; 2^0.5+.1; 2; central_traversal.Station(end)];

% AVERAGING example 1 - default setting
flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

fig_num = 11;  % Define the figure

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_FindOrthogonalPathVectorsAtStations(...
    stations,central_traversal,flag_rounding_type,fig_num);
title('Vertex projection via prior segment (default, flag=1)');

% AVERAGING example 2 - use following segment
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use averagae projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

fig_num = 12;  % Define the figure

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_FindOrthogonalPathVectorsAtStations(...
    stations,central_traversal,flag_rounding_type,fig_num);
title('Vertex projection via following segment (flag=2)');

% AVERAGING example 3 - use average of both segments
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

fig_num = 13;  % Define the figure

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_FindOrthogonalPathVectorsAtStations(...
    stations,central_traversal,flag_rounding_type,fig_num);title('Vertex projection via averaging prior and following segment at vertex (flag=3)');

% AVERAGING example 4 - use average always
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

fig_num = 14;  % Define the figure

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_FindOrthogonalPathVectorsAtStations(...
    stations,central_traversal,flag_rounding_type,fig_num);title('Vertex projection via averaging everywhere (flag=4)');



%% NEGATIVE examples

% Prep the example and workspace
close all;
central_path = [-2 1; 1 4; 3 2];
central_traversal = fcn_Path_convertXYtoTraversalStructure(central_path(:,1),central_path(:,2));
stations = [0; 1.5; 3; 3.5; 18^0.5-0.1; 18^0.5; 18^0.5+.1; 5; 5.5; 6.5; central_traversal.Station(end)];

% NEGATIVE example 1 - default setting
flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

fig_num = 21;  % Define the figure

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_FindOrthogonalPathVectorsAtStations(...
    stations,central_traversal,flag_rounding_type,fig_num);
title('Vertex projection via prior segment (default, flag=1)');

% NEGATIVE example 2 - using following
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

fig_num = 22;  % Define the figure

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_FindOrthogonalPathVectorsAtStations(...
    stations,central_traversal,flag_rounding_type,fig_num);
title('Vertex projection via following segment (flag=2)');

% NEGATIVE example 3 - using average at apex only
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

fig_num = 23;  % Define the figure

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_FindOrthogonalPathVectorsAtStations(...
    stations,central_traversal,flag_rounding_type,fig_num);
title('Vertex projection via averaging prior and following segment at vertex (flag=3)');

% NEGATIVE example 4 - using average always
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

fig_num = 24;  % Define the figure

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_FindOrthogonalPathVectorsAtStations(...
    stations,central_traversal,flag_rounding_type,fig_num);
title('Vertex projection via averaging everywhere (flag=4)');


%% MULTICROSS examples
close all;

% Setup
central_path = [-2 1; 1 4; 3 2; 5 2; 6 3; 7 2];
central_traversal = fcn_Path_convertXYtoTraversalStructure(central_path(:,1),central_path(:,2));
step_size = 0.2;
stations = sort([(0:step_size:central_traversal.Station(end))'; central_traversal.Station]);
stations = unique(stations);

% MULTICROSS example 1 - default setting
flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

fig_num = 31;  % Define the figure

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_FindOrthogonalPathVectorsAtStations(...
    stations,central_traversal,flag_rounding_type,fig_num);
title('Multicross example using projection via prior segment (default)');

% MULTICROSS example 2 - using following
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

fig_num = 32;  % Define the figure

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_FindOrthogonalPathVectorsAtStations(...
    stations,central_traversal,flag_rounding_type,fig_num);
title('Multicross example using projection via following segment');

% MULTICROSS example 3 - using average at apex only
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
% flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

fig_num = 33;  % Define the figure

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_FindOrthogonalPathVectorsAtStations(...
    stations,central_traversal,flag_rounding_type,fig_num);
title('Multicross example using projection via averaging of prior and following segment only at apex');


% MULTICROSS example 4 - using average always
% flag_rounding_type = 1;  % use orthogonal projection of prior segment
% flag_rounding_type = 2;  % use orthogonal projection of following segment
% flag_rounding_type = 3;  % use average projection of prior and following segment, only at endpoints
flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

fig_num = 44;  % Define the figure

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_FindOrthogonalPathVectorsAtStations(...
    stations,central_traversal,flag_rounding_type,fig_num);
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
stations = (0:step_size:central_traversal.Station(end))';
stations = unique(stations);

for i_Path = 1:length(paths)
    nearby_traversal  = data.traversal{i_Path};


    fig_num = 40+i_Path;  % Define the figure
    
    % Calculate the unit normal vectors at given stations and put results into
    % the figure.
    [unit_normal_vector_start, unit_normal_vector_end] = ...
        fcn_Path_FindOrthogonalPathVectorsAtStations(...
        stations,central_traversal,flag_rounding_type,fig_num);
end



% 
% function print_results(stations,closest_path_point,distances)
% fprintf(1,'\n\nStation \t Location X \t Location Y \t Distance \n');
% for i_station =1:length(stations)
%     fprintf(1,'%.2f \t\t %.2f \t\t\t %.2f \t\t\t %.2f\n',stations(i_station),closest_path_point(i_station,1),closest_path_point(i_station,2),distances(i_station));
% end
% end
