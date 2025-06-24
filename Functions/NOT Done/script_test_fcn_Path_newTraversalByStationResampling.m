% script_test_fcn_fcn_Path_convertPathToTraversalStructure
% Tests the function: fcn_Path_convertPathToTraversalStructure

% Revision history
%     2022_01_05: 
%     -- wrote the code originally - in support of
%     fcn_Path_findAverageTraversalViaOrthoProjection

close all;


%% Basic example 1 - start at zero
% Fill in sample paths (as a starter)
basic_path = [0 0; 10 0; 20 0];
input_traversal = fcn_Path_convertPathToTraversalStructure(basic_path);

fig_num = 2344;

% Redecimate the traversal at 1-meter increments
interval = 5;
new_stations    = (0:interval:5)';
new_traversal = fcn_Path_newTraversalByStationResampling(input_traversal, new_stations, fig_num);

% Make sure function worked!
assert(isequal(round(new_traversal.X,4),[0; 5]));

%% Basic example 2 - start at one
% Fill in sample paths (as a starter)
basic_path = [0 0; 10 0; 20 0];
input_traversal = fcn_Path_convertPathToTraversalStructure(basic_path);

fig_num = 2345;

% Redecimate the traversal at 1-meter increments
interval = 5;
new_stations    = (1:interval:6)';
new_traversal = fcn_Path_newTraversalByStationResampling(input_traversal, new_stations, fig_num);

% Make sure function worked!
assert(isequal(round(new_traversal.X,4),[1; 6]));

%% Basic example 3 - first point before start
% Fill in sample paths (as a starter)
basic_path = [0 0; 10 0; 20 0];
input_traversal = fcn_Path_convertPathToTraversalStructure(basic_path);

fig_num = 23456;

% Redecimate the traversal
interval = 5;
new_stations    = (-5:interval:0)';
new_traversal = fcn_Path_newTraversalByStationResampling(input_traversal, new_stations, fig_num);

% Make sure function worked!
assert(isequal(round(new_traversal.X,4),new_stations));

%% Basic example 4 - perfect overlap
% Fill in sample paths (as a starter)
basic_path = [0 0; 10 0; 20 0];
input_traversal = fcn_Path_convertPathToTraversalStructure(basic_path);

fig_num = 23457;

% Redecimate the traversal at 1-meter increments
interval = 10;
new_stations    = (0:interval:20)';
new_traversal = fcn_Path_newTraversalByStationResampling(input_traversal, new_stations, fig_num);

% Make sure function worked!
assert(isequal(round(new_traversal.X,4),new_stations));

%% Basic example 5 - none on traversal
% Fill in sample paths (as a starter)
basic_path = [0 0; 10 0; 20 0];
input_traversal = fcn_Path_convertPathToTraversalStructure(basic_path);

fig_num = 23457;

% Redecimate the traversal at 1-meter increments
interval = 5;
new_stations    = (-10:interval:-5)';
new_traversal = fcn_Path_newTraversalByStationResampling(input_traversal, new_stations, fig_num);

% Make sure function worked!
assert(isequal(round(new_traversal.X,4),new_stations));


%% Advanced test: show this works with a real traversal

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;
input_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

fig_num = 2345;

% Redecimate the traversal at 1-meter increments
interval = 10;
new_stations    = (0:interval:input_traversal.Station(end))';
new_traversal = fcn_Path_newTraversalByStationResampling(input_traversal, new_stations, fig_num);
