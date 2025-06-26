% script_test_fcn_Path_findTraversalStationSegment.m
% This is a script to exercise the function: 
% fcn_Path_findTraversalStationSegment.m
% This function was written on 2020_11_16 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

%      [traversal_trimmed,flag_outside_start, flag_outside_end] = ...
%      fcn_Path_findTraversalStationSegment(...
%      long_traversal, s_coord_start,s_coord_end, 
%      (fig_num))

% Revision history:
%     2021_01_09
%     -- updated name and types to take traversal inputs
%     -- added input checking
%     -- added flag_do_plots    

close all;
   
%% BASIC example: demo of fcn_Path_findTraversalStationSegment
fig_num = 10001;
titleString = sprintf('BASIC example: demo of fcn_Path_findTraversalStationSegment');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

s_coord_start = 10;
s_coord_end   = 100;

[traversal_trimmed,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, (fig_num));

title(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(traversal_trimmed));
assert(isfield(traversal_trimmed,'X'))
assert(isfield(traversal_trimmed,'Y'))
assert(isfield(traversal_trimmed,'Z'))
assert(isfield(traversal_trimmed,'Diff'))
assert(isfield(traversal_trimmed,'Station'))
assert(isfield(traversal_trimmed,'Yaw'))
assert(isnumeric(flag_outside_start));
assert(isnumeric(flag_outside_end));

% Check variable sizes
NreferencePoints = length(traversal_trimmed.X(:,1));
assert(isequal(size(traversal_trimmed.X),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Y),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Z),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Diff),[NreferencePoints 2]));
assert(isequal(size(traversal_trimmed.Station),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Yaw),[NreferencePoints-1 1]));
assert(isequal(size(flag_outside_start),[1 1]));
assert(isequal(size(flag_outside_end),[1 1]));

% Check variable values. Note: the start point should start at 0. The end
% point will round up to end of whatever segment is cut by s_coord_end
assert(traversal_trimmed.Station(1,1)==0);
assert(traversal_trimmed.Station(end,1)<=(abs(s_coord_end-s_coord_start)+20));
assert(flag_outside_start == 0);
assert(flag_outside_end == 0);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example: changing start and stop
fig_num = 10002;
titleString = sprintf('BASIC example: changing start and stop');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

s_coord_start = 40;
s_coord_end   = 100;

[traversal_trimmed,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, (fig_num));

title(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(traversal_trimmed));
assert(isfield(traversal_trimmed,'X'))
assert(isfield(traversal_trimmed,'Y'))
assert(isfield(traversal_trimmed,'Z'))
assert(isfield(traversal_trimmed,'Diff'))
assert(isfield(traversal_trimmed,'Station'))
assert(isfield(traversal_trimmed,'Yaw'))
assert(isnumeric(flag_outside_start));
assert(isnumeric(flag_outside_end));

% Check variable sizes
NreferencePoints = length(traversal_trimmed.X(:,1));
assert(isequal(size(traversal_trimmed.X),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Y),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Z),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Diff),[NreferencePoints 2]));
assert(isequal(size(traversal_trimmed.Station),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Yaw),[NreferencePoints-1 1]));
assert(isequal(size(flag_outside_start),[1 1]));
assert(isequal(size(flag_outside_end),[1 1]));

% Check variable values. Note: the start point should start at 0. The end
% point will round up to end of whatever segment is cut by s_coord_end
assert(traversal_trimmed.Station(1,1)==0);
assert(traversal_trimmed.Station(end,1)<=(abs(s_coord_end-s_coord_start)+20));
assert(flag_outside_start == 0);
assert(flag_outside_end == 0);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example: tighter start and stop
fig_num = 10003;
titleString = sprintf('BASIC example: tighter start and stop');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

s_coord_start = 70;
s_coord_end   = 80;

[traversal_trimmed,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, (fig_num));

title(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(traversal_trimmed));
assert(isfield(traversal_trimmed,'X'))
assert(isfield(traversal_trimmed,'Y'))
assert(isfield(traversal_trimmed,'Z'))
assert(isfield(traversal_trimmed,'Diff'))
assert(isfield(traversal_trimmed,'Station'))
assert(isfield(traversal_trimmed,'Yaw'))
assert(isnumeric(flag_outside_start));
assert(isnumeric(flag_outside_end));

% Check variable sizes
NreferencePoints = length(traversal_trimmed.X(:,1));
assert(isequal(size(traversal_trimmed.X),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Y),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Z),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Diff),[NreferencePoints 2]));
assert(isequal(size(traversal_trimmed.Station),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Yaw),[NreferencePoints-1 1]));
assert(isequal(size(flag_outside_start),[1 1]));
assert(isequal(size(flag_outside_end),[1 1]));

% Check variable values. Note: the start point should start at 0. The end
% point will round up to end of whatever segment is cut by s_coord_end
assert(traversal_trimmed.Station(1,1)==0);
assert(traversal_trimmed.Station(end,1)<=(abs(s_coord_end-s_coord_start)+20));
assert(flag_outside_start == 0);
assert(flag_outside_end == 0);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example: entire traversal 
% because start is less than 0, end greater than max station
fig_num = 10004;
titleString = sprintf('BASIC example: entire traversal ');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

s_coord_start = -5;
s_coord_end   = 7000;

[traversal_trimmed,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, (fig_num));

title(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(traversal_trimmed));
assert(isfield(traversal_trimmed,'X'))
assert(isfield(traversal_trimmed,'Y'))
assert(isfield(traversal_trimmed,'Z'))
assert(isfield(traversal_trimmed,'Diff'))
assert(isfield(traversal_trimmed,'Station'))
assert(isfield(traversal_trimmed,'Yaw'))
assert(isnumeric(flag_outside_start));
assert(isnumeric(flag_outside_end));

% Check variable sizes
NreferencePoints = length(traversal_trimmed.X(:,1));
assert(isequal(size(traversal_trimmed.X),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Y),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Z),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Diff),[NreferencePoints 2]));
assert(isequal(size(traversal_trimmed.Station),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Yaw),[NreferencePoints-1 1]));
assert(isequal(size(flag_outside_start),[1 1]));
assert(isequal(size(flag_outside_end),[1 1]));

% Check variable values. Note: the start point should start at 0. The end
% point will round up to end of whatever segment is cut by s_coord_end
assert(traversal_trimmed.Station(1,1)==0);
assert(traversal_trimmed.Station(end,1)<=(abs(s_coord_end-s_coord_start)+20));
assert(flag_outside_start == 1);
assert(flag_outside_end == 1);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Degenerate cases start here
% all start with the number 2

close all;

%% BASIC example: degenerate case where start = stop
fig_num = 20001;
titleString = sprintf('BASIC example: degenerate case where start = stop');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

s_coord_start = 70;
s_coord_end   = 70;

[traversal_trimmed,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, (fig_num));

title(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(traversal_trimmed));
assert(isfield(traversal_trimmed,'X'))
assert(isfield(traversal_trimmed,'Y'))
assert(isfield(traversal_trimmed,'Z'))
assert(isfield(traversal_trimmed,'Diff'))
assert(isfield(traversal_trimmed,'Station'))
assert(isfield(traversal_trimmed,'Yaw'))
assert(isnumeric(flag_outside_start));
assert(isnumeric(flag_outside_end));

% Check variable sizes
NreferencePoints = length(traversal_trimmed.X(:,1));
assert(isequal(size(traversal_trimmed.X),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Y),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Z),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Diff),[NreferencePoints 2]));
assert(isequal(size(traversal_trimmed.Station),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Yaw),[NreferencePoints-1 1]));
assert(isequal(size(flag_outside_start),[1 1]));
assert(isequal(size(flag_outside_end),[1 1]));

% Check variable values. Note: the start point should start at 0. The end
% point will round up to end of whatever segment is cut by s_coord_end
assert(traversal_trimmed.Station(1,1)==0);
assert(traversal_trimmed.Station(end,1)<=(abs(s_coord_end-s_coord_start)+20));
assert(flag_outside_start == 0);
assert(flag_outside_end == 0);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% BASIC example: beyond end of traversal
fig_num = 20002;
titleString = sprintf('BASIC example: beyond end of traversal');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

s_coord_start = 200;
s_coord_end   = 7000;

[traversal_trimmed,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, (fig_num));

title(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(traversal_trimmed));
assert(isfield(traversal_trimmed,'X'))
assert(isfield(traversal_trimmed,'Y'))
assert(isfield(traversal_trimmed,'Z'))
assert(isfield(traversal_trimmed,'Diff'))
assert(isfield(traversal_trimmed,'Station'))
assert(isfield(traversal_trimmed,'Yaw'))
assert(isnumeric(flag_outside_start));
assert(isnumeric(flag_outside_end));

% Check variable sizes
NreferencePoints = length(traversal_trimmed.X(:,1));
assert(isequal(size(traversal_trimmed.X),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Y),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Z),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Diff),[NreferencePoints 2]));
assert(isequal(size(traversal_trimmed.Station),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Yaw),[NreferencePoints-1 1]));
assert(isequal(size(flag_outside_start),[1 1]));
assert(isequal(size(flag_outside_end),[1 1]));

% Check variable values. Note: the start point should start at 0. The end
% point will round up to end of whatever segment is cut by s_coord_end
assert(traversal_trimmed.Station(1,1)==0);
assert(traversal_trimmed.Station(end,1)<=(abs(s_coord_end-s_coord_start)+20));
assert(flag_outside_start == 0);
assert(flag_outside_end == 1);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example: before start of traversal
fig_num = 20003;
titleString = sprintf('BASIC example: before start of traversal');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

s_coord_start = -200;
s_coord_end   = 100;

[traversal_trimmed,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, (fig_num));

title(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(traversal_trimmed));
assert(isfield(traversal_trimmed,'X'))
assert(isfield(traversal_trimmed,'Y'))
assert(isfield(traversal_trimmed,'Z'))
assert(isfield(traversal_trimmed,'Diff'))
assert(isfield(traversal_trimmed,'Station'))
assert(isfield(traversal_trimmed,'Yaw'))
assert(isnumeric(flag_outside_start));
assert(isnumeric(flag_outside_end));

% Check variable sizes
NreferencePoints = length(traversal_trimmed.X(:,1));
assert(isequal(size(traversal_trimmed.X),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Y),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Z),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Diff),[NreferencePoints 2]));
assert(isequal(size(traversal_trimmed.Station),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Yaw),[NreferencePoints-1 1]));
assert(isequal(size(flag_outside_start),[1 1]));
assert(isequal(size(flag_outside_end),[1 1]));

% Check variable values. Note: the start point should start at 0. The end
% point will round up to end of whatever segment is cut by s_coord_end
assert(traversal_trimmed.Station(1,1)==0);
assert(traversal_trimmed.Station(end,1)<=(abs(s_coord_end-s_coord_start)+20));
assert(flag_outside_start == 1);
assert(flag_outside_end == 0);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% BASIC example: all stations before start
fig_num = 20004;
titleString = sprintf('BASIC example: all stations before start');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

s_coord_start = -200;
s_coord_end   = -100;

[traversal_trimmed,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, (fig_num));

title(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(traversal_trimmed));
assert(isfield(traversal_trimmed,'X'))
assert(isfield(traversal_trimmed,'Y'))
assert(isfield(traversal_trimmed,'Z'))
assert(isfield(traversal_trimmed,'Diff'))
assert(isfield(traversal_trimmed,'Station'))
assert(isfield(traversal_trimmed,'Yaw'))
assert(isnumeric(flag_outside_start));
assert(isnumeric(flag_outside_end));

% Check variable sizes
NreferencePoints = length(traversal_trimmed.X(:,1));
assert(isequal(size(traversal_trimmed.X),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Y),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Z),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Diff),[NreferencePoints 2]));
assert(isequal(size(traversal_trimmed.Station),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Yaw),[NreferencePoints-1 1]));
assert(isequal(size(flag_outside_start),[1 1]));
assert(isequal(size(flag_outside_end),[1 1]));

% Check variable values. Note: the start point should start at 0. The end
% point will round up to end of whatever segment is cut by s_coord_end
assert(traversal_trimmed.Station(1,1)==0);
assert(traversal_trimmed.Station(end,1)<=(abs(s_coord_end-s_coord_start)+20));
assert(flag_outside_start == 1);
assert(flag_outside_end == 0);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example: all stations beyond end
fig_num = 20005;
titleString = sprintf('BASIC example: all stations beyond end');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

s_coord_start = 2000;
s_coord_end   = 3000;

[traversal_trimmed,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, (fig_num));

title(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(traversal_trimmed));
assert(isfield(traversal_trimmed,'X'))
assert(isfield(traversal_trimmed,'Y'))
assert(isfield(traversal_trimmed,'Z'))
assert(isfield(traversal_trimmed,'Diff'))
assert(isfield(traversal_trimmed,'Station'))
assert(isfield(traversal_trimmed,'Yaw'))
assert(isnumeric(flag_outside_start));
assert(isnumeric(flag_outside_end));

% Check variable sizes
NreferencePoints = length(traversal_trimmed.X(:,1));
assert(isequal(size(traversal_trimmed.X),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Y),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Z),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Diff),[NreferencePoints 2]));
assert(isequal(size(traversal_trimmed.Station),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Yaw),[NreferencePoints-1 1]));
assert(isequal(size(flag_outside_start),[1 1]));
assert(isequal(size(flag_outside_end),[1 1]));

% Check variable values. Note: the start point should start at 0. The end
% point will round up to end of whatever segment is cut by s_coord_end
assert(traversal_trimmed.Station(1,1)==0);
assert(traversal_trimmed.Station(end,1)<=(abs(s_coord_end-s_coord_start)+20));
assert(flag_outside_start == 0);
assert(flag_outside_end == 1);

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

% Convert paths to traversal structures
traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

s_coord_start = 10;
s_coord_end   = 100;

[traversal_trimmed,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, ([]));

title(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(traversal_trimmed));
assert(isfield(traversal_trimmed,'X'))
assert(isfield(traversal_trimmed,'Y'))
assert(isfield(traversal_trimmed,'Z'))
assert(isfield(traversal_trimmed,'Diff'))
assert(isfield(traversal_trimmed,'Station'))
assert(isfield(traversal_trimmed,'Yaw'))
assert(isnumeric(flag_outside_start));
assert(isnumeric(flag_outside_end));

% Check variable sizes
NreferencePoints = length(traversal_trimmed.X(:,1));
assert(isequal(size(traversal_trimmed.X),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Y),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Z),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Diff),[NreferencePoints 2]));
assert(isequal(size(traversal_trimmed.Station),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Yaw),[NreferencePoints-1 1]));
assert(isequal(size(flag_outside_start),[1 1]));
assert(isequal(size(flag_outside_end),[1 1]));

% Check variable values. Note: the start point should start at 0. The end
% point will round up to end of whatever segment is cut by s_coord_end
assert(traversal_trimmed.Station(1,1)==0);
assert(traversal_trimmed.Station(end,1)<=(abs(s_coord_end-s_coord_start)+20));
assert(flag_outside_start == 0);
assert(flag_outside_end == 0);

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

s_coord_start = 10;
s_coord_end   = 100;

[traversal_trimmed,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, (-1));

% Check variable types
assert(isstruct(traversal_trimmed));
assert(isfield(traversal_trimmed,'X'))
assert(isfield(traversal_trimmed,'Y'))
assert(isfield(traversal_trimmed,'Z'))
assert(isfield(traversal_trimmed,'Diff'))
assert(isfield(traversal_trimmed,'Station'))
assert(isfield(traversal_trimmed,'Yaw'))
assert(isnumeric(flag_outside_start));
assert(isnumeric(flag_outside_end));

% Check variable sizes
NreferencePoints = length(traversal_trimmed.X(:,1));
assert(isequal(size(traversal_trimmed.X),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Y),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Z),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Diff),[NreferencePoints 2]));
assert(isequal(size(traversal_trimmed.Station),[NreferencePoints 1]));
assert(isequal(size(traversal_trimmed.Yaw),[NreferencePoints-1 1]));
assert(isequal(size(flag_outside_start),[1 1]));
assert(isequal(size(flag_outside_end),[1 1]));

% Check variable values. Note: the start point should start at 0. The end
% point will round up to end of whatever segment is cut by s_coord_end
assert(traversal_trimmed.Station(1,1)==0);
assert(traversal_trimmed.Station(end,1)<=(abs(s_coord_end-s_coord_start)+20));
assert(flag_outside_start == 0);
assert(flag_outside_end == 0);

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

% Convert paths to traversal structures
traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

s_coord_start = 10;
s_coord_end   = 100;

Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [traversal_trimmed,flag_outside_start, flag_outside_end] = ...
        fcn_Path_findTraversalStationSegment(...
        traversal, s_coord_start,s_coord_end, ([]));

end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [traversal_trimmed,flag_outside_start, flag_outside_end] = ...
        fcn_Path_findTraversalStationSegment(...
        traversal, s_coord_start,s_coord_end, (-1));

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

if 1==0
    %% Intentional warnings
    fprintf(1,'\n\nTHE FOLLOWING WILL INTENTIONALLY PRODUCE ERRORS OR WARNINGS\n');

    %% BASIC example: warning thrown because start and end are out of order
    fig_num = 90001;
    titleString = sprintf('BASIC example: warning thrown because start and end are out of order');
    fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
    figure(fig_num); clf;

    % Fill in sample paths (as a starter)
    paths_array = fcn_Path_fillSamplePaths;

    % Convert paths to traversal structures
    traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});

    s_coord_start = 200;
    s_coord_end   = 100;

    [traversal_trimmed,flag_outside_start, flag_outside_end] = ...
        fcn_Path_findTraversalStationSegment(...
        traversal, s_coord_start,s_coord_end, (fig_num));

    title(titleString, 'Interpreter','none');

    % Check variable types
    assert(isstruct(traversal_trimmed));
    assert(isfield(traversal_trimmed,'X'))
    assert(isfield(traversal_trimmed,'Y'))
    assert(isfield(traversal_trimmed,'Z'))
    assert(isfield(traversal_trimmed,'Diff'))
    assert(isfield(traversal_trimmed,'Station'))
    assert(isfield(traversal_trimmed,'Yaw'))
    assert(isnumeric(flag_outside_start));
    assert(isnumeric(flag_outside_end));

    % Check variable sizes
    NreferencePoints = length(traversal_trimmed.X(:,1));
    assert(isequal(size(traversal_trimmed.X),[NreferencePoints 1]));
    assert(isequal(size(traversal_trimmed.Y),[NreferencePoints 1]));
    assert(isequal(size(traversal_trimmed.Z),[NreferencePoints 1]));
    assert(isequal(size(traversal_trimmed.Diff),[NreferencePoints 2]));
    assert(isequal(size(traversal_trimmed.Station),[NreferencePoints 1]));
    assert(isequal(size(traversal_trimmed.Yaw),[NreferencePoints-1 1]));
    assert(isequal(size(flag_outside_start),[1 1]));
    assert(isequal(size(flag_outside_end),[1 1]));

    % Check variable values. Note: the start point should start at 0. The end
    % point will round up to end of whatever segment is cut by s_coord_end
    assert(traversal_trimmed.Station(1,1)==0);
    assert(traversal_trimmed.Station(end,1)<=(abs(s_coord_end-s_coord_start)+20));
    assert(flag_outside_start == 0);
    assert(flag_outside_end == 0);

    % Make sure plot opened up
    assert(isequal(get(gcf,'Number'),fig_num));
end

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
