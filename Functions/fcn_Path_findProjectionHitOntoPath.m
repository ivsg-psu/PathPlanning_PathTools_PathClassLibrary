function [distance,location,path_segment, t, u] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path, ...
    sensor_vector_start, ...
    sensor_vector_end, ...
    varargin)   
%% fcn_Path_findProjectionHitOntoPath 
% calculates hits between a sensor projection and a path, returning the
% distance and location of the hit. Also returns as optional outputs which
% path segment number was hit (e.g. segment 1, 2, 3, etc.) and the distance
% along both that specific path segment (t) and the distance in the sensor
% direction (u).
%
% FORMAT: 
%
%      [distance, location, path_segment, t, u] = ...
%         fcn_Path_findProjectionHitOntoPath(path,...
%         sensor_vector_start,sensor_vector_end,...
%         (flag_search_type),(fig_num))  
%
% INPUTS:
%
%      path: an N x 2 vector containing the X,Y points of the path to be
%      checked for intersections
%
%      sensor_vector_start: a 1 x 2 vector containing the X,Y points of the
%      sensor's start location
%
%      sensor_vector_end: a 1 x 2 vector containing the X,Y points of the
%      sensor's end location
%
%      (OPTIONAL INPUTS)
%      flag_search_type: an integer specifying the type of search.
%
%            0: returns distance and location of first intersection only if
%            the given sensor_vector overlaps the path (this is the
%            default)
%
%            1: returns distance and location of FIRST intersection if ANY
%            projection of the sensor vector, in any direction, hits the
%            path (in other words, if there is any intersection). Note that
%            distance returned could be negative if the nearest
%            intersection is in the opposite direction of the given sensor
%            vector.
%
%            2: returns distances and locations of ALL the detected
%            intersections of where the given sensor_vector overlaps the
%            path (e.g., this gives ALL the results of the flag=0 case).
%            Outputs results as M x 1 and M x 2 vectors respectively, where
%            the M rows represent the ordered intersections.  In cases
%            where the sensor vector completely overlaps a path segment and
%            thus there are infinite points, only the start and end of
%            overlap are given. Note that distance returned will always be
%            positive because only the given sensor vector is checked.
%
%            3: returns distance and location of the FIRST intersection if
%            ANY projection of the path segments, in any direction, hits
%            the sensor (in other words, if there is any intersection).
%            This is the opposite behavior of the flag=1 case. Note that
%            distance returned will always be positive because only the
%            given sensor vector is checked.
%
%            4: returns distance and location of the FIRST intersection of
%            any projection of the path segment vector, in any direction,
%            hits the sensor or if any projection of the sensor vector, in
%            any direction, hits the path segment (in other words, if there
%            is any intersection of the lines that fit the sensor and every
%            path segment). Note that distance returned will be negative if
%            the nearest intersection is in the opposite direction of the
%            given sensor vector.  If multple segments hit at the
%            same intersection, the segment number of the first segment is
%            returned.
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%      distance: a 1 x 1 scalar representing the distance to the closest
%      intersection of the sensor with the path
%
%      location: a 1 x 2 vector of the X,Y location of intersection point
%
%      path_segment: the segment number of the path that was hit (1 is the
%      first segment, 2 is the second, etc)
%
%       t and u: t is distance along the path, and u is distance
%       along the sensor, as fractions of the input vector lengths.    
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%      
%       See the script: script_test_fcn_Path_findProjectionHitOntoPath.m
%       for a full test suite. 
%
% Adopted from https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
% This function was written on 2020_11_14 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2020_11_14 - S. Brennan
% - wrote the code
% 2020_12_29 - S. Brennan
% - fixed bug with missing sensor vector if starts on a path
% - added better comments
% 2020_12_31 - S. Brennan
% - fixed the input arguments so that they are more clear (start/end)
% - fixed the incorrect location of the debug echo at top of the code
% - fixed flag usage to decouple plotting with debugging
% 2021_01_08
% -- Added input check on path type
% 2021_01_23 - S. Brennan
% -- Added flag_search_type = 2 option, to allow multiple cross points
% to be returned
% -- Fixed bug with partially overlapping vectors not returning a
% result
% -- Added path segment output so that we can ID which segment was hit
% 2021_01_24 - S. Brennan
% -- Fixed bug with overlapping colinear where two path segments
% identified when there is only one
% 2021_12_27 - S. Brennan
% -- Added better comments on flags
% 2024_03_14 - S. Brennan
% -- Added better comments
% -- Fixed bug where the figure plotting breaks if someone gives an
% empty figure number
% -- Added flag 3 and 4 cases
% 2024_05_14 - Aneesh Batchu
% -- Added max speed options
% 2025_06_23 - S. Brennan
% -- changed dependency to only use findSensorHitOnWall

%% Set up for debugging

flag_max_speed = 0;
if (nargin==5 && isequal(varargin{end},-1))
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS");
    MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG = getenv("MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS);
    end
end


if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 34838; %#ok<NASGU>
else
    debug_fig_num = []; %#ok<NASGU>
end

% 
% flag_do_debug = 0; % Flag to plot the results for debugging
% flag_check_inputs = 1; % Flag to perform input checking
% 
% if flag_do_debug
%     st = dbstack; %#ok<*UNRCH>
%     fprintf(1,'Starting function: %s, in file: %s\n',st(1).name,st(1).file);
% end

%% check input arguments
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
        narginchk(3,5);

        % Check path input
        fcn_DebugTools_checkInputsToFunctions(path, 'path');

    end
end

% Does user wish to specify flag_search_type type?
flag_search_type = 0;
if 4 <= nargin
    flag_search_type = varargin{1};
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (5<= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end



%% Calculations begin here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define each path as a set of walls that can be hit
wall_start = [path(1:end-1,1) path(1:end-1,2)];
wall_end   = [path(2:end,1) path(2:end,2)];

switch flag_search_type
    case {0}
        % 0: returns distance and location of first intersection only if
        % the given sensor_vector overlaps the path (this is the
        % default)
        flag_search_return_type = 0;
        flag_search_range_type = 0;

    case {1}
        % 1: returns distance and location of FIRST intersection if ANY
        % projection of the sensor vector, in any direction, hits the
        % path (in other words, if there is any intersection). Note that
        % distance returned could be negative if the nearest
        % intersection is in the opposite direction of the given sensor
        % vector.
        flag_search_return_type = 0;
        flag_search_range_type = 1;

    case {2}       
        % 2: returns distances and locations of ALL the detected
        % intersections of where the given sensor_vector overlaps the
        % path (e.g., this gives ALL the results of the flag=0 case).
        % Outputs results as M x 1 and M x 2 vectors respectively, where
        % the M rows represent the ordered intersections.  In cases
        % where the sensor vector completely overlaps a path segment and
        % thus there are infinite points, only the start and end of
        % overlap are given. Note that distance returned will always be
        % positive because only the given sensor vector is checked.
        flag_search_return_type = 1;
        flag_search_range_type = 0;

    case {3}
        % 3: returns distance and location of the FIRST intersection if
        % ANY projection of the path segments, in any direction, hits
        % the sensor (in other words, if there is any intersection).
        % This is the opposite behavior of the flag=1 case. Note that
        % distance returned will always be positive because only the
        % given sensor vector is checked.
        flag_search_return_type = 0;
        flag_search_range_type = 2;
        
    case {4}        
        % 4: returns distance and location of the FIRST intersection of
        % any projection of the path segment vector, in any direction,
        % hits the sensor or if any projection of the sensor vector, in
        % any direction, hits the path segment (in other words, if there
        % is any intersection of the lines that fit the sensor and every
        % path segment). Note that distance returned will be negative if
        % the nearest intersection is in the opposite direction of the
        % given sensor vector.  If multple segments hit at the
        % same intersection, the segment number of the first segment is
        % returned.
        flag_search_return_type = 0;
        flag_search_range_type  = 3;

    otherwise
        error('unrecognized value');
end

tolerance = [];

[distance, location, path_segment, t, u] = ...
    fcn_Path_findSensorHitOnWall(...
    wall_start, wall_end,...
    sensor_vector_start,sensor_vector_end,...
    (flag_search_return_type), (flag_search_range_type), ...
    (tolerance), -1);


%% Any debugging?
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
    
    % Set up the figure
    figure(fig_num);
    clf;
    hold on;
    axis equal;
    grid on; grid minor;
    
    % Plot the path in black
    plot(path(:,1),path(:,2),'k.-','Linewidth',5);
    handle_text = text(path(1,1),path(1,2),'Path');
    set(handle_text,'Color',[0 0 0]);
       
    % Plot the sensor vector
    quiver(sensor_vector_start(:,1),sensor_vector_start(:,2),sensor_vector_end(:,1)-sensor_vector_start(:,1),sensor_vector_end(:,2)-sensor_vector_start(:,2),0,'r','Linewidth',3);
    plot(sensor_vector_end(:,1),sensor_vector_end(:,2),'r.','Markersize',10);
    
    handle_text = text(sensor_vector_start(:,1),sensor_vector_start(:,2),'Sensor');
    set(handle_text,'Color',[1 0 0]);
    
    axis_size = axis;
    y_range = axis_size(4)-axis_size(3);
        
    % Plot any hits in blue
    for i_result = 1:length(distance)
        plot(location(i_result,1),location(i_result,2),'bo','Markersize',30);
        handle_text = text(location(i_result,1),location(i_result,2)-0.05*y_range,sprintf('Hit at distance: %.2f',distance(i_result)));
        set(handle_text,'Color',[0 0 1]);
    end
    
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end
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




