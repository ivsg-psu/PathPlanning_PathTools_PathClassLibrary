function [distance,location,wall_segment, t, u] = ...
    fcn_Path_findProjectionHitOntoPathImproved(...
    wall_start, wall_end, ...
    sensor_vector_start, sensor_vector_end, ...
    varargin)   
%% fcn_Path_findProjectionHitOntoPathImproved 
% calculates hits between a sensor projection and a wall, returning the
% distance and location of the hit. Also returns as optional outputs which
% wall segment number was hit (e.g. segment 1, 2, 3, etc.) and the distance
% along both that specific wall segment (t) and the distance in the sensor
% direction (u).
%
% FORMAT: 
%
%      [distance, location, wall_segment, t, u] = ...
%         fcn_Path_findProjectionHitOntoPathImproved(...
%         wall_start, wall_end,...
%         sensor_vector_start,sensor_vector_end,...
%         (flag_search_return_type), (flag_search_range_type), ...
%         (tolerance), (fig_num))  
%
% INPUTS:
%
%      wall_start: an N x 2 vector containing the X,Y points of the
%      starting points of segments to be tested
%
%      wall_end: an N x 2 vector containing the X,Y points of the
%      end points of segments to be tested
%
%      sensor_vector_start: a 1 x 2 vector containing the X,Y points of the
%      sensor's start location
%
%      sensor_vector_end: a 1 x 2 vector containing the X,Y points of the
%      sensor's end location
%
%      (OPTIONAL INPUTS)
%
%      flag_search_return_type: an integer specifying how many points to
%      return
%
%            0: (default) returns distance and location of FIRST
%            intersection only if the given sensor_vector overlaps the
%            wall. Returns [nan nan] if no overlap found.  In cases where
%            the sensor vector completely overlaps a wall segment and thus
%            there are infinite points, only the starting point of overlap
%            is given. If a sensor starts on top of a wall segment, the
%            sensor's starting point is returned. If a sensor overlaps
%            multiple wall segments at the start, only the first wall
%            segment is given.
%
%            1: returns distance and location of ALL intersections where
%            the given sensor_vector overlaps the wall. Returns [nan nan]
%            if no overlap found.  In cases where the sensor vector
%            completely overlaps a wall segment and thus there are infinite
%            points, the starting point of overlap is given. If a sensor
%            starts on top of a wall segment, the sensor's starting point
%            and ending points are returned. If a sensor overlaps multiple
%            wall segments, the first and last point of overlap for each
%            wall segment is given. Outputs results as M x 1 and M x 2
%            vectors respectively, where the M rows represent the ordered
%            intersections. Some intersections may be repeated locations
%            because if they occur on different wall segments.
%
%      flag_search_range_type: an integer specifying if the range of the
%      sensor and/or wall segments to be evaluated. When projections are
%      used, the range extends in the vector direction of the sensor and/or
%      wall, in both the negative and positive directions. Distances in the
%      negative sensor direction are reported as negative values.
%
%            0: (default) the GIVEN sensor and GIVEN wall used.
%            1: ANY projection of the sensor used with the GIVEN wall
%            2: ANY projection of the wall used with the GIVEN sensor
%            3: ANY projection of BOTH the wall and sensor
%
%      tolerance: (default is 1000*eps) How close points must be to a
%      segment to be counted as intersecting. Positive values are
%      inclusive, negative values are exclusive. For example, consider 3
%      wall segments:
%            A: from the origin [0 0] to [1 0], 
%            B: from [0.5 0] to [0.5 1],
%            C: from [1 0] to [2 0]
%      if the tolerance is 0, whether or not A and B intersect can be
%      unclear and is strongly affected by the numerical accuray of number
%      representations. If the tolerance is higher than numberical accuracy
%      amounts (say, 1000 times larger than this value or 1000*eps), then A
%      and B definitely intersect at 1 point that is -almost- exactly in
%      the middle of A, and -almost- exactly at the start of B. However,
%      for the same tolerance, A and C overlap at infinite points.
%      Similarly, for negative tolerances that are sufficiently large, A
%      and B do NOT overlap, nor do A and C. The use of tolerance allows
%      the desirable definition of "near miss" situations.
%
%      fig_num: a figure number to plot results. Turns debugging on.
%
% OUTPUTS:
%
%      distance: a 1 x 1 scalar representing the distance to the closest
%      intersection of the sensor with the wall
%
%      location: a 1 x 2 vector of the X,Y location of intersection point
%
%      wall_segment: the segment number of the wall that was hit (1 is the
%      first segment, 2 is the second, etc)
%
%       t and u: t is distance along the wall, and u is distance
%       along the sensor, as fractions of the input vector lengths.    
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%      
%       See the script: script_test_fcn_Path_findProjectionHitOntoPathImproved.m
%       for a full test suite. 
%
% Adopted from https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
% This function was written in original form on 2020_11_14 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%      2020_11_14 - S. Brennan
%      - wrote the code
%      2020_12_29 - S. Brennan
%      - fixed bug with missing sensor vector if starts on a wall
%      - added better comments
%      2020_12_31 - S. Brennan 
%      - fixed the input arguments so that they are more clear (start/end)
%      - fixed the incorrect location of the debug echo at top of the code
%      - fixed flag usage to decouple plotting with debugging
%      2021_01_08 
%      -- Added input check on wall type
%      2021_01_23 - S. Brennan
%      -- Added flag_search_type = 2 option, to allow multiple cross points
%      to be returned
%      -- Fixed bug with partially overlapping vectors not returning a
%      result
%      -- Added wall segment output so that we can ID which segment was hit
%      2021_01_24 - S. Brennan
%      -- Fixed bug with overlapping colinear where two wall segments
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
%      -- changed flag inputs to give more options and clearer usage
%      -- expanded out wall inputs and outputs to allow more general use
%      -- added tolerances

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==8 && isequal(varargin{end},-1))
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
        narginchk(4,8);

        % Check the wall_start input
        fcn_DebugTools_checkInputsToFunctions(wall_start, '2column_of_numbers');

        % Check the wall_end input
        fcn_DebugTools_checkInputsToFunctions(wall_end, '2column_of_numbers');

        % Check the wall_start input
        fcn_DebugTools_checkInputsToFunctions(sensor_vector_start, '2column_of_numbers');

        % Check the wall_end input
        fcn_DebugTools_checkInputsToFunctions(sensor_vector_end, '2column_of_numbers');

    end
end

% Is the user entering a flag_search_return_type?
flag_search_return_type = 0;
if 5 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        flag_search_return_type = temp;
    end
end

% Is the user entering a flag_search_range_type?
flag_search_range_type = 0;
if 6 <= nargin
    temp = varargin{2};
    if ~isempty(temp)
        flag_search_range_type = temp;
    end
end

% Is the user entering a tolerance?
tolerance = eps*1000;
if 7 <= nargin
    temp = varargin{3};
    if ~isempty(temp)
        tolerance = temp;
    end
end

% Does user want to show the plots?
flag_do_plot = 0; % Default is to NOT show plots
if (0==flag_max_speed) && (8 == nargin) 
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
near_zero = tolerance;

% Define wall indexes
wall_indexes = (1:length(wall_start(:,1)))';

% Define p, q, r and s vectors. Per the reference, these have the following
% meanings:
% p: the [x y] location of the start of each "wall". Can be Nx2
% q: the [x y] location of the start of the "sensor". Must be 1x2
% r: the [x y] vector denoting the transition from start of a wall to end
% s: the [x y] vector denoting the transition from start of sensor to end

% To solve for:
% t: the scalar denoting the percentage along the wall vector of intersect
% u: the scalar denoting the percentage along the sensor vector of intersect
% Solution: setup up p + t*r = q + u*s

p = wall_start;
q = sensor_vector_start;
r = wall_end - wall_start;
s = sensor_vector_end - sensor_vector_start;

% Define useful intermediate terms. CrossProduct if function defined below.
r_cross_s = crossProduct(r,s);
q_minus_p =  q - p;
q_minus_p_cross_s = crossProduct(q_minus_p,s);
q_minus_p_cross_r = crossProduct(q_minus_p,r);

%% Step 1: classify all the solutions

% Find which walls are parallel to the sensor.
% If r × s = 0, then the two lines are parallel
flag_isParallel = (near_zero>=abs(r_cross_s));

% Are any of these parallel and non-intersecting?
% If r × s = 0 and (q − p) × r ≠ 0, then the two lines are parallel and non-intersecting
parallel_non_intersecting_indices = find(flag_isParallel.*(near_zero<abs(q_minus_p_cross_r)));
if any(parallel_non_intersecting_indices)
    % If they are colinear or parallel, make dummy length that is removed
    % later, to avoid divide by zero in intermediate calculations.
    r_cross_s(parallel_non_intersecting_indices) = 1; 
end


% Are any of these collinear (e.g. on same line)? If they are collinear, then
% need to find whether there is overlapping, e.g. there may be infinite
% solutions. In other words, we need to know when the t-value, the
% percentage along each wall, starts (t0) and stops (t1). For a given wall,
% these values can never be larger than 1 or less than 0, physcially. So we
% need to "cap" t0 and t1 at these values.
%
% If r × s = 0 and (q − p) × r = 0, then the two lines are collinear. 
% NOTE: later steps require division by r x s, so we need to be careful to
% set this value to 1 to avoid division by zero.
% 
% In this case, express the endpoints of the second segment (q and q + s)
% in terms of the equation of the first line segment (p + t r):
% 
% t0 = (q − p) · r / (r · r)
% 
% t1 = (q + s − p) · r / (r · r) = t0 + s · r / (r · r)
% 
% If the interval between t0 and t1 intersects the interval [0, 1] then the
% line segments are collinear and overlapping; otherwise they are collinear
% and disjoint.
% 
% Note that if s and r point in opposite directions, then s · r < 0 and so
% the interval to be checked is [t1, t0] rather than [t0, t1].

colinear_indices = find(flag_isParallel.*(near_zero>=abs(q_minus_p_cross_r)));

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
    %     fprintf(1,'any_within flag:\n');
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



% Calculate t and u, where t is distance along the wall, and u is distance
% along the sensor.

t = q_minus_p_cross_s./r_cross_s; % Distance along the wall
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
    wall_indexes = [wall_indexes; wall_indexes(more_indices)];

end


%% Apply flag_search_range_type
% 0: (default) the GIVEN sensor and GIVEN wall used.
% 1: ANY projection of the sensor used with the GIVEN wall
% 2: ANY projection of the wall used with the GIVEN sensor
% 3: ANY projection of BOTH the wall and sensor
% 
% To apply these, we impose restrictions on:
% t: the scalar denoting the percentage along the wall vector of intersect
% u: the scalar denoting the percentage along the sensor vector of intersect
%
% Note: Since doing many segments at once, need to use vector form (e.g.
% the .* format of dot products).
% 
% Note: Tolerance added as numerical errors can cause points to be missed
% for some segments that right next to or through points. This biases -
% very slightly - the data to include intersections along segments that
% "graze" next to each other. For example, the segment from (0,0) to
% (1,0) barely grazes the segment from (0.5,0) to (0.5,1).

zero_threshold = 0 - tolerance;  % Positive tolerance assumed
one_threshold  = 1 + tolerance;
if 0 == flag_search_range_type
    % 0: (default) the GIVEN sensor and GIVEN wall used.
    % This constrains both t (wall extent) and u (sensor extent)
    good_vector = ((zero_threshold<=t).*(one_threshold>=t).*(zero_threshold<=u).*(one_threshold>=u));
elseif 1 == flag_search_range_type
    % 1: ANY projection of the sensor used with the GIVEN wall
    % This constrains ONLY t (wall extent) 
    good_vector = ((zero_threshold<=t).*(one_threshold>=t));
elseif 2 == flag_search_range_type
    % 2: ANY projection of the wall used with the GIVEN sensor
    % This constrains ONLY u (sensor extent) 
    good_vector = ((zero_threshold<=u).*(one_threshold>=u));
elseif 3 == flag_search_range_type
    % 3: ANY projection of BOTH the wall and sensor
    % NO constraints
    good_vector = ((Inf>u).*(Inf>t));
else
    warning('on','backtrace');
    warning('Expecting a flag_search_range_type as integer in range of 0 to 3, but found: %.3f',flag_search_range_type);
    error('Bad flag_search_range_type encountered');
end

% Keep only the indices that work
good_indices = find(good_vector>0);

% Initialize all intersections to infinity
intersections = NaN*ones(length(p(:,1)),2);
if ~isempty(good_indices)
    % Calculate the intersection point (finally)
    result = p + t.*r; 
    intersections(good_indices,:) = result(good_indices,:);    
end


%% Apply flag_search_return_type
% 0: (default) returns results of FIRST intersection
% 1: returns distance and location of ALL intersections


% Find the distances via Euclidian distance to the sensor's origin
% note: a faster way to do this might be to just
% calculate t*r as a length
distances_squared = sum((intersections - sensor_vector_start).^2,2);

if 0==flag_search_return_type
    % Keep only the minimum distance
    [closest_distance_squared,closest_index] = min(distances_squared);
    
    distance = closest_distance_squared^0.5*sign(u(closest_index));
    location = intersections(closest_index,:);
    wall_segment = wall_indexes(closest_index);

    if isnan(distance)
        wall_segment = nan;
    end
elseif 1==flag_search_return_type
    % Return all the results
    if all(isnan(distances_squared))
        distance = nan;
        location = [nan nan];
        wall_segment = nan;
    else
        good_indices = find(~isnan(distances_squared));
        distance = distances_squared(good_indices).^0.5.*sign(u(good_indices));
        location = intersections(good_indices,:);
        wall_segment = wall_indexes(good_indices);
    end
else
    warning('on','backtrace');
    warning('Expecting a flag_search_return_type as integer with values of 0 or 1, but found: %.3f',flag_search_return_type);
    error('Bad flag_search_range_type encountered');
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
    
    % Find all the walls in one plottable format
    Nwalls = length(wall_start(:,1));
    allWallsX = [wall_start(:,1) wall_end(:,1) nan(Nwalls,1)];
    allWallsX = reshape(allWallsX',1,[]);
    allWallsY = [wall_start(:,2) wall_end(:,2) nan(Nwalls,1)];
    allWallsY = reshape(allWallsY',1,[]);
    allWalls = [allWallsX' allWallsY'];


    % Find size of plotting domain
    allPoints = [allWalls; sensor_vector_start; sensor_vector_end];

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

    % Plot the walls in black   
    plot(allWalls(:,1),allWalls(:,2),'k.-','Linewidth',5);
    handle_text = text(allWalls(1,1),allWalls(1,2),'Walls');
    set(handle_text,'Color',[0 0 0]);

    % Plot the sensor vector
    quiver(q(:,1),q(:,2),s(:,1),s(:,2),'r','Linewidth',2,'MaxHeadSize',1);
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



