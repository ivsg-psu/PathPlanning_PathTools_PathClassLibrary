% script_test_fcn_Path_fillRandomTraversalsAboutTraversal.m
% Tests fcn_Path_fillRandomTraversalsAboutTraversal
       
% Revision history:
% 2021_01_03
% -- first write of the code
% 2021_01_07
% -- cleared up function calls for traversals vs paths
% 2021_01_09
% -- cleared up a few function names in the script
close all





%% Test case: basic call for one trajectory
fig_num = 10001;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});
emtpy_value = [];

random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    emtpy_value,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    emtpy_value,... % (flag_generate_random_stations),...
    emtpy_value,... % (spatial_smoothness),...
    fig_num);

% Check variable types
assert(isstruct(random_traversals));
assert(isfield(random_traversals,'traversal'))
assert(iscell(random_traversals.traversal))
assert(isfield(random_traversals.traversal{1},'X'))
assert(isfield(random_traversals.traversal{1},'Y'))
assert(isfield(random_traversals.traversal{1},'Z'))
assert(isfield(random_traversals.traversal{1},'Diff'))
assert(isfield(random_traversals.traversal{1},'Station'))
assert(isfield(random_traversals.traversal{1},'Yaw'))

% Check variable sizes
NreferencePoints = length(reference_traversal.X(:,1));
assert(length(random_traversals.traversal)==1) %#ok<ISCL>
assert(isequal(size(random_traversals.traversal{1}.X),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Y),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Z),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Diff),[NreferencePoints 2]));
assert(isequal(size(random_traversals.traversal{1}.Station),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Yaw),[NreferencePoints-1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Test case: advanced call for one trajectory - specify figure
fig_num = 10002;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

%      random_traversals = ...
%      fcn_Path_fillRandomTraversalsAboutTraversal(...
%            reference_traversal,...
%            (std_deviation),...
%            (num_trajectories),...
%            (num_points),...
%            (flag_generate_random_stations),...
%            (spatial_smoothness),...
%            (fig_num));
emtpy_value = [];

random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    emtpy_value,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    emtpy_value,... % (flag_generate_random_stations),...
    emtpy_value,... % (spatial_smoothness),...
    fig_num);

% Check variable types
assert(isstruct(random_traversals));
assert(isfield(random_traversals,'traversal'))
assert(iscell(random_traversals.traversal))
assert(isfield(random_traversals.traversal{1},'X'))
assert(isfield(random_traversals.traversal{1},'Y'))
assert(isfield(random_traversals.traversal{1},'Z'))
assert(isfield(random_traversals.traversal{1},'Diff'))
assert(isfield(random_traversals.traversal{1},'Station'))
assert(isfield(random_traversals.traversal{1},'Yaw'))

% Check variable sizes
NreferencePoints = length(reference_traversal.X(:,1));
assert(length(random_traversals.traversal)==1) %#ok<ISCL>
assert(isequal(size(random_traversals.traversal{1}.X),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Y),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Z),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Diff),[NreferencePoints 2]));
assert(isequal(size(random_traversals.traversal{1}.Station),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Yaw),[NreferencePoints-1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Test case: use same station points
fig_num = 10003;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;


% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

%      random_traversals = ...
%      fcn_Path_fillRandomTraversalsAboutTraversal(...
%            reference_traversal,...
%            (std_deviation),...
%            (num_trajectories),...
%            (num_points),...
%            (flag_generate_random_stations),...
%            (spatial_smoothness),...
%            (fig_num));
emtpy_value = [];

flag_generate_random_stations = 0;
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    emtpy_value,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    emtpy_value,... % (spatial_smoothness),...
    fig_num);

% Check variable types
assert(isstruct(random_traversals));
assert(isfield(random_traversals,'traversal'))
assert(iscell(random_traversals.traversal))
assert(isfield(random_traversals.traversal{1},'X'))
assert(isfield(random_traversals.traversal{1},'Y'))
assert(isfield(random_traversals.traversal{1},'Z'))
assert(isfield(random_traversals.traversal{1},'Diff'))
assert(isfield(random_traversals.traversal{1},'Station'))
assert(isfield(random_traversals.traversal{1},'Yaw'))

% Check variable sizes
NreferencePoints = length(reference_traversal.X(:,1));
assert(length(random_traversals.traversal)==1) %#ok<ISCL>
assert(isequal(size(random_traversals.traversal{1}.X),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Y),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Z),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Diff),[NreferencePoints 2]));
assert(isequal(size(random_traversals.traversal{1}.Station),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Yaw),[NreferencePoints-1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Test case: show effects of spatial smoothness with many trajectories
fig_num = 10004;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

%      random_traversals = ...
%      fcn_Path_fillRandomTraversalsAboutTraversal(...
%            reference_traversal,...
%            (std_deviation),...
%            (num_trajectories),...
%            (num_points),...
%            (flag_generate_random_stations),...
%            (spatial_smoothness),...
%            (fig_num));
emtpy_value = [];
flag_generate_random_stations = 0;
num_trajectories = 5;

subplot(2,2,1);
spatial_smoothness = 4;  % Units are meters
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title('Spatial smoothness: 4 meters (below 3 generates warning)');

% Check variable types
assert(isstruct(random_traversals));
assert(isfield(random_traversals,'traversal'))
assert(iscell(random_traversals.traversal))
assert(isfield(random_traversals.traversal{1},'X'))
assert(isfield(random_traversals.traversal{1},'Y'))
assert(isfield(random_traversals.traversal{1},'Z'))
assert(isfield(random_traversals.traversal{1},'Diff'))
assert(isfield(random_traversals.traversal{1},'Station'))
assert(isfield(random_traversals.traversal{1},'Yaw'))

% Check variable sizes
NreferencePoints = length(reference_traversal.X(:,1));
assert(length(random_traversals.traversal)==num_trajectories) 
assert(isequal(size(random_traversals.traversal{1}.X),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Y),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Z),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Diff),[NreferencePoints 2]));
assert(isequal(size(random_traversals.traversal{1}.Station),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Yaw),[NreferencePoints-1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

subplot(2,2,2);
spatial_smoothness = 8;  % Units are meters
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('Spatial smoothness: %.0d meters',spatial_smoothness));

% Check variable types
assert(isstruct(random_traversals));
assert(isfield(random_traversals,'traversal'))
assert(iscell(random_traversals.traversal))
assert(isfield(random_traversals.traversal{1},'X'))
assert(isfield(random_traversals.traversal{1},'Y'))
assert(isfield(random_traversals.traversal{1},'Z'))
assert(isfield(random_traversals.traversal{1},'Diff'))
assert(isfield(random_traversals.traversal{1},'Station'))
assert(isfield(random_traversals.traversal{1},'Yaw'))

% Check variable sizes
NreferencePoints = length(reference_traversal.X(:,1));
assert(length(random_traversals.traversal)==num_trajectories)
assert(isequal(size(random_traversals.traversal{1}.X),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Y),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Z),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Diff),[NreferencePoints 2]));
assert(isequal(size(random_traversals.traversal{1}.Station),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Yaw),[NreferencePoints-1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

subplot(2,2,3);
spatial_smoothness = 20;  % Units are meters
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('Spatial smoothness: %.0d meters',spatial_smoothness));

% Check variable types
assert(isstruct(random_traversals));
assert(isfield(random_traversals,'traversal'))
assert(iscell(random_traversals.traversal))
assert(isfield(random_traversals.traversal{1},'X'))
assert(isfield(random_traversals.traversal{1},'Y'))
assert(isfield(random_traversals.traversal{1},'Z'))
assert(isfield(random_traversals.traversal{1},'Diff'))
assert(isfield(random_traversals.traversal{1},'Station'))
assert(isfield(random_traversals.traversal{1},'Yaw'))

% Check variable sizes
NreferencePoints = length(reference_traversal.X(:,1));
assert(length(random_traversals.traversal)==num_trajectories)
assert(isequal(size(random_traversals.traversal{1}.X),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Y),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Z),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Diff),[NreferencePoints 2]));
assert(isequal(size(random_traversals.traversal{1}.Station),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Yaw),[NreferencePoints-1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

subplot(2,2,4);
spatial_smoothness = 40;  % Units are meters
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('Spatial smoothness: %.0d meters',spatial_smoothness));

% Check variable types
assert(isstruct(random_traversals));
assert(isfield(random_traversals,'traversal'))
assert(iscell(random_traversals.traversal))
assert(isfield(random_traversals.traversal{1},'X'))
assert(isfield(random_traversals.traversal{1},'Y'))
assert(isfield(random_traversals.traversal{1},'Z'))
assert(isfield(random_traversals.traversal{1},'Diff'))
assert(isfield(random_traversals.traversal{1},'Station'))
assert(isfield(random_traversals.traversal{1},'Yaw'))

% Check variable sizes
NreferencePoints = length(reference_traversal.X(:,1));
assert(length(random_traversals.traversal)==num_trajectories)
assert(isequal(size(random_traversals.traversal{1}.X),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Y),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Z),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Diff),[NreferencePoints 2]));
assert(isequal(size(random_traversals.traversal{1}.Station),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Yaw),[NreferencePoints-1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Test case: show effects of standard deviation
fig_num = 10005;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;


% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});


%      random_traversals = ...
%      fcn_Path_fillRandomTraversalsAboutTraversal(...
%            reference_traversal,...
%            (std_deviation),...
%            (num_trajectories),...
%            (num_points),...
%            (flag_generate_random_stations),...
%            (spatial_smoothness),...
%            (fig_num));
emtpy_value = [];
flag_generate_random_stations = 0;
num_trajectories = 5;
spatial_smoothness = 7;  % Units are meters

subplot(3,1,1)
std_deviation = 1;  % Units are meters
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    std_deviation,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('Standard deviation: %.0d meters',std_deviation));

% Check variable types
assert(isstruct(random_traversals));
assert(isfield(random_traversals,'traversal'))
assert(iscell(random_traversals.traversal))
assert(isfield(random_traversals.traversal{1},'X'))
assert(isfield(random_traversals.traversal{1},'Y'))
assert(isfield(random_traversals.traversal{1},'Z'))
assert(isfield(random_traversals.traversal{1},'Diff'))
assert(isfield(random_traversals.traversal{1},'Station'))
assert(isfield(random_traversals.traversal{1},'Yaw'))

% Check variable sizes
NreferencePoints = length(reference_traversal.X(:,1));
assert(length(random_traversals.traversal)==num_trajectories)
assert(isequal(size(random_traversals.traversal{1}.X),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Y),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Z),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Diff),[NreferencePoints 2]));
assert(isequal(size(random_traversals.traversal{1}.Station),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Yaw),[NreferencePoints-1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

subplot(3,1,2)
std_deviation = 2;  % Units are meters
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    std_deviation,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('Standard deviation: %.0d meters',std_deviation));

% Check variable types
assert(isstruct(random_traversals));
assert(isfield(random_traversals,'traversal'))
assert(iscell(random_traversals.traversal))
assert(isfield(random_traversals.traversal{1},'X'))
assert(isfield(random_traversals.traversal{1},'Y'))
assert(isfield(random_traversals.traversal{1},'Z'))
assert(isfield(random_traversals.traversal{1},'Diff'))
assert(isfield(random_traversals.traversal{1},'Station'))
assert(isfield(random_traversals.traversal{1},'Yaw'))

% Check variable sizes
NreferencePoints = length(reference_traversal.X(:,1));
assert(length(random_traversals.traversal)==num_trajectories)
assert(isequal(size(random_traversals.traversal{1}.X),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Y),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Z),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Diff),[NreferencePoints 2]));
assert(isequal(size(random_traversals.traversal{1}.Station),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Yaw),[NreferencePoints-1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

subplot(3,1,3)
std_deviation = 5;  % Units are meters
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    std_deviation,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('Standard deviation: %.0d meters',std_deviation));

% Check variable types
assert(isstruct(random_traversals));
assert(isfield(random_traversals,'traversal'))
assert(iscell(random_traversals.traversal))
assert(isfield(random_traversals.traversal{1},'X'))
assert(isfield(random_traversals.traversal{1},'Y'))
assert(isfield(random_traversals.traversal{1},'Z'))
assert(isfield(random_traversals.traversal{1},'Diff'))
assert(isfield(random_traversals.traversal{1},'Station'))
assert(isfield(random_traversals.traversal{1},'Yaw'))

% Check variable sizes
NreferencePoints = length(reference_traversal.X(:,1));
assert(length(random_traversals.traversal)==num_trajectories)
assert(isequal(size(random_traversals.traversal{1}.X),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Y),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Z),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Diff),[NreferencePoints 2]));
assert(isequal(size(random_traversals.traversal{1}.Station),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Yaw),[NreferencePoints-1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Test case: show effects of num_points
fig_num = 10006;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

%      random_traversals = ...
%      fcn_Path_fillRandomTraversalsAboutTraversal(...
%            reference_traversal,...
%            (std_deviation),...
%            (num_trajectories),...
%            (num_points),...
%            (flag_generate_random_stations),...
%            (spatial_smoothness),...
%            (fig_num));
emtpy_value = [];
flag_generate_random_stations = 1;
num_trajectories = 5;
spatial_smoothness = 10;  % Units are meters

subplot(3,1,1);
num_points = 10;
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    num_points,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('num points: %.0d',num_points));

% Check variable types
assert(isstruct(random_traversals));
assert(isfield(random_traversals,'traversal'))
assert(iscell(random_traversals.traversal))
assert(isfield(random_traversals.traversal{1},'X'))
assert(isfield(random_traversals.traversal{1},'Y'))
assert(isfield(random_traversals.traversal{1},'Z'))
assert(isfield(random_traversals.traversal{1},'Diff'))
assert(isfield(random_traversals.traversal{1},'Station'))
assert(isfield(random_traversals.traversal{1},'Yaw'))

% Check variable sizes
NreferencePoints = num_points;
assert(length(random_traversals.traversal)==num_trajectories) 
assert(isequal(size(random_traversals.traversal{1}.X),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Y),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Z),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Diff),[NreferencePoints 2]));
assert(isequal(size(random_traversals.traversal{1}.Station),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Yaw),[NreferencePoints-1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

subplot(3,1,2);
num_points = 50;
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    num_points,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('num points: %.0d',num_points));

% Check variable types
assert(isstruct(random_traversals));
assert(isfield(random_traversals,'traversal'))
assert(iscell(random_traversals.traversal))
assert(isfield(random_traversals.traversal{1},'X'))
assert(isfield(random_traversals.traversal{1},'Y'))
assert(isfield(random_traversals.traversal{1},'Z'))
assert(isfield(random_traversals.traversal{1},'Diff'))
assert(isfield(random_traversals.traversal{1},'Station'))
assert(isfield(random_traversals.traversal{1},'Yaw'))

% Check variable sizes
NreferencePoints = num_points;
assert(length(random_traversals.traversal)==num_trajectories) 
assert(isequal(size(random_traversals.traversal{1}.X),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Y),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Z),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Diff),[NreferencePoints 2]));
assert(isequal(size(random_traversals.traversal{1}.Station),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Yaw),[NreferencePoints-1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

subplot(3,1,3);
num_points = 200;
random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    num_trajectories,... % (num_trajectories),...
    num_points,... % (num_points),...
    flag_generate_random_stations,... % (flag_generate_random_stations),...
    spatial_smoothness,... % (spatial_smoothness),...
    fig_num);
title(sprintf('num points: %.0d',num_points));

% Check variable types
assert(isstruct(random_traversals));
assert(isfield(random_traversals,'traversal'))
assert(iscell(random_traversals.traversal))
assert(isfield(random_traversals.traversal{1},'X'))
assert(isfield(random_traversals.traversal{1},'Y'))
assert(isfield(random_traversals.traversal{1},'Z'))
assert(isfield(random_traversals.traversal{1},'Diff'))
assert(isfield(random_traversals.traversal{1},'Station'))
assert(isfield(random_traversals.traversal{1},'Yaw'))

% Check variable sizes
NreferencePoints = num_points;
assert(length(random_traversals.traversal)==num_trajectories)
assert(isequal(size(random_traversals.traversal{1}.X),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Y),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Z),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Diff),[NreferencePoints 2]));
assert(isequal(size(random_traversals.traversal{1}.Station),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Yaw),[NreferencePoints-1 1]));

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

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});
emtpy_value = [];

random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    emtpy_value,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    emtpy_value,... % (flag_generate_random_stations),...
    emtpy_value,... % (spatial_smoothness),...
    []);

% Check variable types
assert(isstruct(random_traversals));
assert(isfield(random_traversals,'traversal'))
assert(iscell(random_traversals.traversal))
assert(isfield(random_traversals.traversal{1},'X'))
assert(isfield(random_traversals.traversal{1},'Y'))
assert(isfield(random_traversals.traversal{1},'Z'))
assert(isfield(random_traversals.traversal{1},'Diff'))
assert(isfield(random_traversals.traversal{1},'Station'))
assert(isfield(random_traversals.traversal{1},'Yaw'))

% Check variable sizes
NreferencePoints = length(reference_traversal.X(:,1));
assert(length(random_traversals.traversal)==1) %#ok<ISCL>
assert(isequal(size(random_traversals.traversal{1}.X),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Y),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Z),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Diff),[NreferencePoints 2]));
assert(isequal(size(random_traversals.traversal{1}.Station),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Yaw),[NreferencePoints-1 1]));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});
emtpy_value = [];

random_traversals = ...
    fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
    emtpy_value,... % (std_deviation),...
    emtpy_value,... % (num_trajectories),...
    emtpy_value,... % (num_points),...
    emtpy_value,... % (flag_generate_random_stations),...
    emtpy_value,... % (spatial_smoothness),...
    -1);

% Check variable types
assert(isstruct(random_traversals));
assert(isfield(random_traversals,'traversal'))
assert(iscell(random_traversals.traversal))
assert(isfield(random_traversals.traversal{1},'X'))
assert(isfield(random_traversals.traversal{1},'Y'))
assert(isfield(random_traversals.traversal{1},'Z'))
assert(isfield(random_traversals.traversal{1},'Diff'))
assert(isfield(random_traversals.traversal{1},'Station'))
assert(isfield(random_traversals.traversal{1},'Yaw'))

% Check variable sizes
NreferencePoints = length(reference_traversal.X(:,1));
assert(length(random_traversals.traversal)==1) %#ok<ISCL>
assert(isequal(size(random_traversals.traversal{1}.X),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Y),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Z),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Diff),[NreferencePoints 2]));
assert(isequal(size(random_traversals.traversal{1}.Station),[NreferencePoints 1]));
assert(isequal(size(random_traversals.traversal{1}.Yaw),[NreferencePoints-1 1]));

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

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});
emtpy_value = [];


Niterations = 20;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    random_traversals = ...
        fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
        emtpy_value,... % (std_deviation),...
        emtpy_value,... % (num_trajectories),...
        emtpy_value,... % (num_points),...
        emtpy_value,... % (flag_generate_random_stations),...
        emtpy_value,... % (spatial_smoothness),...
        []);
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    random_traversals = ...
        fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal,...
        emtpy_value,... % (std_deviation),...
        emtpy_value,... % (num_trajectories),...
        emtpy_value,... % (num_points),...
        emtpy_value,... % (flag_generate_random_stations),...
        emtpy_value,... % (spatial_smoothness),...
        -1);
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
