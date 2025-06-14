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
%            default). Returns [nan nan] if no overlap found.
%
%            1: returns distance and location of FIRST intersection if ANY
%            projection of the sensor vector, in any direction, hits the
%            path (in other words, if there is any intersection). Note that
%            distance returned could be negative if the nearest
%            intersection is in the opposite direction of the given sensor
%            vector. Returns [nan nan] if no overlap found.
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
%      fig_num: a figure number to plot results. Turns debugging on.
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
%      fcn_Path_checkInputsToFunctions
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
%      2020_11_14 - S. Brennan
%      - wrote the code
%      2020_12_29 - S. Brennan
%      - fixed bug with missing sensor vector if starts on a path
%      - added better comments
%      2020_12_31 - S. Brennan 
%      - fixed the input arguments so that they are more clear (start/end)
%      - fixed the incorrect location of the debug echo at top of the code
%      - fixed flag usage to decouple plotting with debugging
%      2021_01_08 
%      -- Added input check on path type
%      2021_01_23 - S. Brennan
%      -- Added flag_search_type = 2 option, to allow multiple cross points
%      to be returned
%      -- Fixed bug with partially overlapping vectors not returning a
%      result
%      -- Added path segment output so that we can ID which segment was hit
%      2021_01_24 - S. Brennan
%      -- Fixed bug with overlapping colinear where two path segments
%      identified when there is only one
%      2021_12_27 - S. Brennan
%      -- Added better comments on flags
%      2024_03_14 - S. Brennan
%      -- Added better comments
%      -- Fixed bug where the figure plotting breaks if someone gives an
%      empty figure number
%      -- Added flag 3 and 4 cases
%      2024_06_19 - S. Brennan
%      -- fixed tolerance issue with overlapping vertical lines
%      2025_06_14 - S. Brennan
%      -- added expanded plotting so not tight
%      -- added fast mode, global flags, better input checking

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==5 && isequal(varargin{end},-1))
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_PATH_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_PATH_FLAG_CHECK_INPUTS");
    MATLABFLAG_PATH_FLAG_DO_DEBUG = getenv("MATLABFLAG_PATH_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_PATH_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_PATH_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_PATH_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_PATH_FLAG_CHECK_INPUTS);
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
        narginchk(3,5);

        % Check the path input
        fcn_DebugTools_checkInputsToFunctions(path, 'path');
    end
end

% Is the user entering a flag_search_type?
flag_search_type = 0;
if 4 <= nargin
    flag_search_type = varargin{1};
end


% Does user want to show the plots?
flag_do_plot = 0; % Default is to NOT show plots
if (0==flag_max_speed) && (5 == nargin) 
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plot = 1;
    end
else
    if flag_do_debug
        fig = figure; 
        fig_num = fig.Number;
        flag_do_plot = 1;
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

% Define meaning of near zero
near_zero = (eps*1000);

% Define each path as a set of walls that can be hit
wall_start = [path(1:end-1,1) path(1:end-1,2)];
wall_end   = [path(2:end,1) path(2:end,2)];
path_segments = (1:length(wall_start(:,1)))';

% Define p, q, r and s vectors
p = wall_start;
q = sensor_vector_start;
r = wall_end - wall_start;
s = sensor_vector_end - sensor_vector_start;

% Define useful intermediate terms
r_cross_s = crossProduct(r,s);
q_minus_p =  q - p;
q_minus_p_cross_s = crossProduct(q_minus_p,s);
q_minus_p_cross_r = crossProduct(q_minus_p,r);

% Are any of these parallel?
parallel_non_intersecting_indices = find((near_zero>=abs(r_cross_s)).*(near_zero<abs(q_minus_p_cross_r)));
if any(parallel_non_intersecting_indices)
    r_cross_s(parallel_non_intersecting_indices) = 1; % They are colinear or parallel, so make dummy length
end

% Are any of these colinear?
colinear_indices = find((near_zero>=abs(r_cross_s)).*(near_zero>=abs(q_minus_p_cross_r)));
if any(colinear_indices)
    r_cross_s(colinear_indices) = 1; % They are colinear or parallel, so make dummy length
    r_dot_r = sum(r.*r,2);
    q_minus_p_dot_r = sum(q_minus_p.*r,2);
    s_dot_r = sum(s.*r,2);
    t0 = q_minus_p_dot_r./r_dot_r;
    t1 = t0 + s_dot_r./r_dot_r;
    
    % Keep only the good indices
    %     % For debugging:
    %     conditions = [-0.5 -0.4; -0.5 0; -0.5 .2; 0 0.2; 0.2 0.4; 0.2 1; 0.2 1.2; 1 1.2; 1.2 1.3; -0.5 1.2]
    %     t0 = conditions(:,1);
    %     t1 = conditions(:,2);
    %     t0_inside = (t0>=0)&(t0<=1);
    %     t1_inside = (t1>=0)&(t1<=1);
    %     t0_t1_surround = (t0<0)&(t1>1) | (t1<0)&(t0>1);
    %     any_within = t0_inside | t1_inside | t0_t1_surround;
    %     fprintf(1,'t0 inside flag:\n');
    %     [conditions t0_inside]
    %     fprintf(1,'t1 inside flag:\n');
    %     [conditions t1_inside]
    %     fprintf(1,'surround flag:\n');
    %     [conditions t0_t1_surround]
    %     fprintf(1,'any_wighin flag:\n');
    %     [conditions any_within]
    
    % Check whether there is overlap by seeing of the t0 and t1 values are
    % within the interval of [0 1], endpoint inclusive
    t0_inside = (t0>=0)&(t0<=1);
    t1_inside = (t1>=0)&(t1<=1);
    t0_t1_surround = (t0<0)&(t1>1) | (t1<0)&(t0>1);
    any_within = t0_inside | t1_inside | t0_t1_surround;
    good_indices = find(any_within);
    good_colinear_indices = intersect(colinear_indices,good_indices);
    
    % Fix the ranges to be within 0 and 1
    t0(good_colinear_indices) = max(0,t0(good_colinear_indices));
    t0(good_colinear_indices) = min(1,t0(good_colinear_indices));

    t1(good_colinear_indices) = max(0,t1(good_colinear_indices));
    t1(good_colinear_indices) = min(1,t1(good_colinear_indices));

end

% Calculate t and u, where t is distance along the path, and u is distance
% along the sensor.

t = q_minus_p_cross_s./r_cross_s; % Distance along the path
u = q_minus_p_cross_r./r_cross_s; % Distance along the sensor

% Fix any situations that are parallel and non-intersecting, as these will
% give wrong calculation results from the t and u calculations above
t(parallel_non_intersecting_indices) = inf;
u(parallel_non_intersecting_indices) = inf;

% Fix any situations that are colinear. For these, we save the start point
% as the point where overlap starts, and the end point where overlap ends.
% Note that this creates NEW intersections beyond the number of segements
if any(colinear_indices)
    % Shut off colinear ones to start
    t(colinear_indices) = inf;
    u(colinear_indices) = inf;
    
    % Correct the t values
    u(good_colinear_indices) = 1;
    t(good_colinear_indices) = t0(good_colinear_indices);

    % Do we need to add more hit points?
    indices_hit_different_point = find(t0~=t1);
    more_indices = intersect(indices_hit_different_point,good_colinear_indices);
         
    % Make p and r, t and u longer so that additional hit points are
    % calculated in special case of overlaps
    p = [p; p(more_indices,:)];
    r = [r; r(more_indices,:)];
    u = [u; u(more_indices)];
    t = [t; t1(more_indices)];   
    path_segments = [path_segments; path_segments(more_indices)];

end

% Initialize all intersections to infinity
intersections = NaN*ones(length(p(:,1)),2);

% Note: Since doing many segments at once, need to use vector form (e.g.
% the .* format of dot products).
% 
% Note: Tolerance added as numerical errors can cause points to be missed
% for some segments that right next to or through points. This biases -
% very slightly - the data to include intersections along segments that
% "graze" next to each other. For example, the segment from (0,0) to
% (1,0) barely grazes the segment from (0.5,0) to (0.5,1).

tolerance = eps*1000;
zero_threshold = 0 - tolerance;
one_threshold  = 1 + tolerance;
if 0 == flag_search_type
    good_vector = ((zero_threshold<=t).*(one_threshold>=t).*(zero_threshold<=u).*(one_threshold>=u));
elseif 1 == flag_search_type
    good_vector = ((zero_threshold<=t).*(one_threshold>=t));
elseif 2 == flag_search_type
    good_vector = ((zero_threshold<=t).*(one_threshold>=t).*(zero_threshold<=u).*(one_threshold>=u));
elseif 3 == flag_search_type           % Changed on 2024_03_14
    good_vector = ((zero_threshold<=u).*(one_threshold>=u));
elseif 4 == flag_search_type           % Changed on 2024_03_14
    good_vector = ((Inf>u).*(Inf>t));
else
    error('Incorrect flag_search_type entered');
end

% Keep only the indices that work
good_indices = find(good_vector>0);

if ~isempty(good_indices)
    result = p + t.*r; 
    intersections(good_indices,:) = result(good_indices,:);    
end

% Find the distances via Euclidian distance to the sensor's origin
% note: a faster way to do this might be to just
% calculate t*r as a length
distances_squared = sum((intersections - sensor_vector_start).^2,2);

if (flag_search_type ==0) || (flag_search_type==1) || (flag_search_type ==3) || (flag_search_type ==4)
    % Keep just the minimum distance
    [closest_distance_squared,closest_index] = min(distances_squared);
    
    distance = closest_distance_squared^0.5*sign(u(closest_index));
    location = intersections(closest_index,:);
    path_segment = path_segments(closest_index);

    if isnan(distance)
        path_segment = nan;
    end
elseif flag_search_type ==2
    % Return all the results
    good_indices = find(~isnan(distances_squared));
    distance = distances_squared(good_indices).^0.5.*sign(u(good_indices));
    location = intersections(good_indices,:);
    path_segment = path_segments(good_indices);

else
    error('unexpected flag_search_type')
end

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
if flag_do_plot
     % check whether the figure already has data
    h_fig = figure(fig_num);
    flag_rescale_axis = 0; 
    if isempty(get(h_fig,'Children')) 
        flag_rescale_axis = 1; 
    else
        child_handle = get(h_fig,'Children');
        if isfield(child_handle,'TileArrangement') && strcmp(get(child_handle,'TileArrangement'),'flow')
            flag_rescale_axis = 1;
        end
    end

    hold on;
    axis equal;
    grid on; grid minor;
    

    % Find size of plotting domain
    allPoints = [path; sensor_vector_start; sensor_vector_end];

    max_plotValues = max(allPoints);
    min_plotValues = min(allPoints);
    sizePlot = max(max_plotValues) - min(min_plotValues);
    nudge = sizePlot*0.006; %#ok<NASGU>

    % Set size of plotting domain
    if flag_rescale_axis
        
        percent_larger = 0.3;
        axis_range = max_plotValues - min_plotValues;
        if (0==axis_range(1,1))
            axis_range(1,1) = 2/percent_larger;
        end
        if (0==axis_range(1,2))
            axis_range(1,2) = 2/percent_larger;
        end
        
        % Force the axis to be equal
        min_vertexValuesInPlot = min(min_plotValues);
        max_vertexValuesInPlot = max(max_plotValues);

        % Stretch the axes
        stretched_min_vertexValues = min_vertexValuesInPlot - percent_larger.*axis_range;
        stretched_max_vertexValues = max_vertexValuesInPlot + percent_larger.*axis_range;
        axesTogether = [stretched_min_vertexValues; stretched_max_vertexValues];
        newAxis = reshape(axesTogether, 1, []);
        axis(newAxis);

    end
    % goodAxis = axis;

    % Plot the path in black
    plot(path(:,1),path(:,2),'k.-','Linewidth',5);
    handle_text = text(path(1,1),path(1,2),'Path');
    set(handle_text,'Color',[0 0 0]);

    % Plot the sensor vector
    quiver(q(:,1),q(:,2),s(:,1),s(:,2),'r','Linewidth',3);
    plot(sensor_vector_start(:,1),sensor_vector_start(:,2),'r.','Markersize',20);
    plot(sensor_vector_end(:,1),sensor_vector_end(:,2),'r.','Markersize',10);

    handle_text = text(q(:,1),q(:,2),'Sensor');
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



%% Calculate cross products
function result = crossProduct(v,w)
result = v(:,1).*w(:,2)-v(:,2).*w(:,1);
end



