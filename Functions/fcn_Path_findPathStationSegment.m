function [path_trimmed, flag_outside_start, flag_outside_end] = ...
    fcn_Path_findPathStationSegment(...
    long_path, s_coord_start, s_coord_end, varargin)
% fcn_Path_findPathStationSegment
% Finds portion of a long_path that contains the given s_coordinates,
% starting from s_coord_start to s_coord_end, and returns that portion as
% another path "trimmed" out of the original by finding the path
% segments within the long path closest to the queried s-coordinates.
% If the s-coordinates are outside those of the long_path, then flags
% are set to 1 for either flag_outside_start, flag_outside_end; otherwise,
% these flags are zero. The end station of the trimemd path will round
% up to end of whatever segment is cut by s_coord_end.
%
% Note that this function always returns at least 2 points representing the
% closest path segment, even if both s-point queries are outside the given
% path.
%
% 
% FORMAT: 
%
%      [path_trimmed,flag_outside_start, flag_outside_end] = ...
%      fcn_Path_findPathStationSegment(...
%      long_path, s_coord_start,s_coord_end, 
%      (fig_num))
%
% INPUTS:
%
%      long_path:  a N x 2 or N x 3 set of coordinates
%      representing the [X Y] or [X Y Z] coordinates, in sequence, of a
%      path that is being used for trimming by the given station coordinates.
%
%      s_coord_start: a 1x1 (scalar) indicating the s-coordinate location
%      at which the query starts. The path segment output will start at
%      previous s-value to this station.
%
%      s_coord_end: a 1x1 (scalar) indicating the s-coordinate location
%      at which the query ends. The path segment output will end at
%      subsequent s-value to this station.
%
%      (OPTIONAL INPUTS)
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%      path_trimmed: a path output trimmed out of the original
%      path. 
%
%      flag_outside_start, flag_outside_end: flags that are set equal to 1
%      if the query is outside the s-distance within the given path at
%      either the start, the end, or both
%
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_Path_calcPathStation
%
% EXAMPLES:
%      
% See the script: 
% script_test_fcn_Path_findPathStationSegment.m
% for a full test suite.
%
% This function was written on 2020_10_14 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2020_10_14:
% - first write of the code
% 2020_11_15:
% - changed the name to prep for Paths class
% 2021_01_08
% -- started updating for new class
% 2021_01_09
% -- updated name and types to take traversal inputs
% -- added input checking
% -- added flag_do_plots
% 2025_06_23 - S. Brennan
% -- Updated debugging and input checks
% 2025_07_01 - S. Brennan
% -- Removed traversal input type and replaced with path types
% -- Renamed function from fcn_Path_findTraversalStationSegment

% TO-DO
% (none)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 4; % The largest Number of argument inputs to the function
flag_max_speed = 0;
if (nargin==MAX_NARGIN && isequal(varargin{end},-1))
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_PATHCLASS_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_PATHCLASS_FLAG_CHECK_INPUTS");
    MATLABFLAG_PATHCLASS_FLAG_DO_DEBUG = getenv("MATLABFLAG_PATHCLASS_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_PATHCLASS_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_PATHCLASS_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_PATHCLASS_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_PATHCLASS_FLAG_CHECK_INPUTS);
    end
end

% flag_do_debug = 1;

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 999978; %#ok<NASGU>
else
    debug_fig_num = []; %#ok<NASGU>
end

%% check input arguments?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0==flag_max_speed
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(3,MAX_NARGIN);

        % Check the long_path input
        fcn_DebugTools_checkInputsToFunctions(long_path, 'path2or3D');

        % Check the s_coord_start input
        fcn_DebugTools_checkInputsToFunctions(s_coord_start, 'station');

        % Check the s_coord_end input
        fcn_DebugTools_checkInputsToFunctions(s_coord_end, 'station');

    end
end

% Check that the start is strictly before the end
if s_coord_start>s_coord_end
    warning('S coordinates of start and end are out of order. These will be automatically corrected, but the results may be incorrect.');
    temp = s_coord_end;
    s_coord_end = s_coord_start;
    s_coord_start = temp;
end

% Does user want to show the plots?
flag_do_plots = 0; % Default is to NOT show plots
if (0==flag_max_speed) && (MAX_NARGIN == nargin) 
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        fig_num = temp;
        figure(fig_num);
        flag_do_plots = 1;
    end
else
    if flag_do_debug
        fig_debug = 272763; %#ok<NASGU>
    end
end


%% Find the closest point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Are the input vectors the right shape?
Npoints_in_path = length(long_path(:,1));

% Set default outputs
flag_outside_start = 0;
flag_outside_end   = 0;

% Use find function to grab values
[long_path_Stations, ~] = fcn_Path_calcPathStation(long_path,-1);
path_segment_start_index = find(long_path_Stations < s_coord_start,1,'last');
path_segment_end_index   = find(long_path_Stations >s_coord_end,1,'first');

% Check if we've gone past the ends of the path
if isempty(path_segment_start_index) % There is no s-coordinate in the path smaller than the start
    path_segment_start_index = 1;
    flag_outside_start = 1;
end
if isempty(path_segment_end_index) % There is no s-coordinate in the path larger than the end
    path_segment_end_index = Npoints_in_path;
    flag_outside_end   = 1;
end

% Check to see if start index is the end of path
% Are all path s-coordinates are smaller than the start?
if Npoints_in_path == path_segment_start_index  
    path_segment_start_index = Npoints_in_path - 1;
    path_segment_end_index = Npoints_in_path;
end

% Check to see if end index is the start of path?
if 1 == path_segment_end_index  % All path s-coordinates are larger than the end
    path_segment_start_index = 1;
    path_segment_end_index = 2;
end

% Check to see if start and end of segment are the same (degenerate case)
if path_segment_start_index == path_segment_end_index
    path_segment_start_index = path_segment_end_index - 1;
end

% Grab the path segment that is closest to the given s-coordinates
path_trimmed = long_path(path_segment_start_index:path_segment_end_index,:); 

%% Plot the results (for debugging)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _                 
%  |  __ \     | |                
%  | |  | | ___| |__  _   _  __ _ 
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_do_plots
 % Prep the figure for plotting
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end
    
    % Is this 2D or 3D?
    dimension_of_points = 2;

    % Find size of plotting domain
    allPointsBeingPlotted = [long_path; path_trimmed];
    max_plotValues = max(allPointsBeingPlotted);
    min_plotValues = min(allPointsBeingPlotted);
    sizePlot = max(max_plotValues) - min(min_plotValues);
    nudge = sizePlot*0.006; %#ok<NASGU>

    % Find size of plotting domain
    if flag_rescale_axis
        percent_larger = 0.3;
        axis_range = max_plotValues - min_plotValues;
        if (0==axis_range(1,1))
            axis_range(1,1) = 2/percent_larger;
        end
        if (0==axis_range(1,2))
            axis_range(1,2) = 2/percent_larger;
        end
        if dimension_of_points==3 && (0==axis_range(1,3))
            axis_range(1,3) = 2/percent_larger;
        end

        % Force the axis to be equal?
        if 1==1
            min_valuesInPlot = min(min_plotValues);
            max_valuesInPlot = max(max_plotValues);
        else
            min_valuesInPlot = min_plotValues;
            max_valuesInPlot = max_plotValues;
        end

        % Stretch the axes
        stretched_min_vertexValues = min_valuesInPlot - percent_larger.*axis_range;
        stretched_max_vertexValues = max_valuesInPlot + percent_larger.*axis_range;
        axesTogether = [stretched_min_vertexValues; stretched_max_vertexValues];
        newAxis = reshape(axesTogether, 1, []);
        axis(newAxis);

    end
    % goodAxis = axis;

    hold on;
    grid on;
    axis equal;

    xlabel('X [m]');
    ylabel('Y [m]');


    % Plot the path
    plot(long_path(:,1),long_path(:,2),'r.-','Linewidth',5, 'MarkerSize',30,'DisplayName','Input path');       
    
    % Plot the results
    plot(path_trimmed(:,1),path_trimmed(:,2),'b.-','Linewidth',3,'MarkerSize',15,'DisplayName','Trimmed path');   

    legend;

    text(path_trimmed(1,1),path_trimmed(1,2),'Start');   
    text(path_trimmed(end,1),path_trimmed(end,2),'End');     
    
end % Ends the flag_do_debug if statement


if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends the function

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
