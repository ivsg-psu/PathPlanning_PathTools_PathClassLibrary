function [distance,location,wall_segment, t, u] = ...
    fcn_Path_findSensorHitOnWall(...
    wall_start, wall_end, ...
    sensor_vector_start, sensor_vector_end, ...
    varargin)
%% fcn_Path_findSensorHitOnWall
% calculates hits between a sensor projection and a wall, returning the
% distance and location of the hit. Also returns as optional outputs which
% wall segment number was hit (e.g. segment 1, 2, 3, etc.) and the distance
% along both that specific wall segment (t) and the distance in the sensor
% direction (u).
%
% FORMAT:
%      [distance, location, wall_segment, t, u] = ...
%         fcn_Path_findSensorHitOnWall(...
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
%            0: GIVEN projection of the sensor used with GIVEN projection of wall (default)
%            1: ANY   projection of the sensor used with GIVEN projection of wall
%            2: GIVEN projection of the sensor used with ANY   projection of wall
%            3: ANY   projection of the sensor used with ANY   projection of wall
%
%     If any walls are entered as a single point (wall_start = wall_end),
%     then the flag_search_range_type means the following:
%            0: GIVEN projection of the sensor hits GIVEN wall point (default)
%            1: ANY   projection of the sensor hits GIVEN wall point
%            2: GIVEN projection of the sensor hits SENSOR ORTHO PROJECTION of wall onto sensor
%            2: ANY   projection of the sensor hits SENSOR ORTHO PROJECTION of wall onto sensor
%     and the t vector refers to the ortho projection distance (as
%     percentage of the sensor), and the u vector has the typical meaning.
%
%     If the sensor is entered as a single point (sensor_vector_start = sensor_vector_end),
%     then the flag_search_range_type means the following:
%            0: GIVEN point of sensor                     is inside GIVEN projection of walls (default)
%            1: ORTHO PROJECTION of sensor onto each wall is inside GIVEN projection of walls
%            2: GIVEN point of sensor                     is inside ANY   projection of walls
%            3: ORTHO PROJECTION of sensor onto each wall is inside ANY   projection of walls
%     and the t vector has the typical meaning, and the u vector is the
%     orthogonal projection distance along each wall
%
%      tolerance: (default is 1000*eps) How close points must be to a
%      segment to be counted as intersecting. Positive values are
%      inclusive, negative values are exclusive. For example, consider 3
%      wall segments:
%            A: from the origin [0 0] to [1 0],
%            B: from [0.5 0] to [0.5 1] (head of B sits on middle of A)
%            C: from [1 0] to [2 0]  (head of C sits on tail of A)
%      if the tolerance is 0, whether or not A and B intersect can be
%      unclear as this is strongly affected by the numerical accuray of number
%      representations. If the tolerance is higher than numberical accuracy
%      amounts (say, 1000 times larger than this value or 1000*eps), then A
%      and B definitely intersect at 1 point that is -almost- exactly in
%      the middle of A, and -almost- exactly at the start of B. However,
%      for the same tolerance, A and C technically overlap at infinite points.
%      Similarly, for negative tolerances that are sufficiently large, A
%      and B do NOT overlap, nor do A and C. The use of tolerance allows
%      the desirable definition of "near miss" situations. The primary
%      effect of tolerance is for when the sensor and wall are collinear -
%      a negative tolerance will ALWAYS have the sensor "miss" the wall.
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
%       See the script: script_test_fcn_Path_findSensorHitOnWall.m
%       for a full test suite.
%
% Adopted from https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
% This function was written in original form on 2020_11_14 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2020_11_14 - S. Brennan
% - wrote the code
% 2020_12_29 - S. Brennan
% - fixed bug with missing sensor vector if starts on a wall
% - added better comments
% 2020_12_31 - S. Brennan
% - fixed the input arguments so that they are more clear (start/end)
% - fixed the incorrect location of the debug echo at top of the code
% - fixed flag usage to decouple plotting with debugging
% 2021_01_08
% -- Added input check on wall type
% 2021_01_23 - S. Brennan
% -- Added flag_search_type = 2 option, to allow multiple cross points
% to be returned
% -- Fixed bug with partially overlapping vectors not returning a
% result
% -- Added wall segment output so that we can ID which segment was hit
% 2021_01_24 - S. Brennan
% -- Fixed bug with overlapping colinear where two wall segments
% identified when there is only one
% 2021_12_27 - S. Brennan
% -- Added better comments on flags
% 2024_03_14 - S. Brennan
% -- Added better comments
% -- Fixed bug where the figure plotting breaks if someone gives an
% empty figure number
% -- Added flag 3 and 4 cases
% 2024_06_19 - S. Brennan
% -- fixed tolerance issue with overlapping vertical lines
% 2025_06_14 - S. Brennan
% -- added expanded plotting so not tight
% -- added fast mode, global flags, better input checking
% -- changed flag inputs to give more options and clearer usage
% -- expanded out wall inputs and outputs to allow more general use
% -- added tolerances
% -- added "point wall" option, where if distance of wall_start to
% wall_end less than tolerance, wall is treated as point.
% -- added "point sensor" option, where if distance of
% sensor_vector_start to sensor_vector_end less than tolerance,
% sensor is treated as point.
% 2025_06_23 - S. Brennan
% -- Updated debugging and input checks

% TO-DO
% (nothing yet)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 8; % The largest Number of argument inputs to the function
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
        narginchk(4,MAX_NARGIN);

        % Check the wall_start input
        fcn_DebugTools_checkInputsToFunctions(wall_start, '2column_of_numbers');

        % Check the wall_end input
        fcn_DebugTools_checkInputsToFunctions(wall_end, '2column_of_numbers');

        % Check the sensor_vector_start input
        fcn_DebugTools_checkInputsToFunctions(sensor_vector_start, '2column_of_numbers',[1 1]);

        % Check the sensor_vector_end input
        fcn_DebugTools_checkInputsToFunctions(sensor_vector_end, '2column_of_numbers',[1 1]);

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
if (0==flag_max_speed) && (MAX_NARGIN == nargin)
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
sensor_vector = sensor_vector_end - sensor_vector_start;


%%%%
% Identify point cases. This will happen with r dot r is (near) zero.
% For these, artificially connect the point to the
% end of the sensor to create a virtual "wall" that can be used hereafter.
flag_isPointWall = sum(r.^2,2)<tolerance.^2;
pointWall_indices = find(flag_isPointWall);
r(pointWall_indices,:) = sensor_vector_end - wall_start(pointWall_indices,:);

% Is the sensor a point?
flag_isPointSensor = sum(sensor_vector.^2,2)<tolerance.^2;

if flag_isPointSensor
    % If a flag is a point sensor, then set the s value to the orthogonal
    % projection of each wall

    % First, make sure there are no point walls
    if any(flag_isPointWall)
        error('Sensor and individual walls cannot both be single points. Exiting.');
    end

    % Take orthogonal projection of the wall vectors. Use negative value
    % because distances measured from sensor to wall
    q = ones(length(r(:,1)),1)*q;
    sensor_vector = -r*[0 1; -1 0];
end



% Define useful intermediate terms. CrossProduct if function defined below.
r_cross_s = fcn_INTERNAL_crossProduct(r,sensor_vector);
q_minus_p =  q - p;
q_minus_p_cross_s = fcn_INTERNAL_crossProduct(q_minus_p,sensor_vector);
q_minus_p_cross_r = fcn_INTERNAL_crossProduct(q_minus_p,r);

%% Step 1: find t and u values for all cases
% Avoid division by zero for parallel cases
% If r × s = 0, then the two lines are parallel
flag_isParallel = (near_zero>=abs(r_cross_s));

if any(flag_isParallel)
    % If any walls are parallel to sensor, make dummy length that is
    % removed later, to avoid divide by zero in intermediate calculations.
    r_cross_s(flag_isParallel) = 1;
end


%%%%%
% General case for t and u calculations
% Calculate t and u, where t is distance along the wall, and u is distance
% along the sensor.

tValues = q_minus_p_cross_s./r_cross_s; % Distance along the wall
uValues = q_minus_p_cross_r./r_cross_s; % Distance along the sensor

%%%%%
% Parallel cases are nan
% NOTE: this includes non-intersecting cases
% If r × s = 0 and (q − p) × r ≠ 0, then the two lines are parallel and non-intersecting
% NOTE: the way to calculate this is as follows. The code that follows does
% not need this variable, but it is included anyway in case future
% modifications are needed:
% parallel_non_intersecting_indices = find(flag_isParallel.*(near_zero<abs(q_minus_p_cross_r)));
%
% Fix any situations that are parallel as these will give wrong calculation
% results from the t and u calculations above
tValues(flag_isParallel) = nan;
uValues(flag_isParallel) = nan;

%%%%
% Fix these collinear cases (e.g. on same line)
% IF wall segments and sensor are collinear, then there may be overlap
% producing infinite intersections.
% If r × s = 0 (e.g., it is parallel) and (q − p) × r = 0, then the two
% lines are collinear.
flag_isCollinear = flag_isParallel.*(near_zero>=abs(q_minus_p_cross_r)).*(~flag_isPointWall);
collinear_indices = find(flag_isCollinear);
if any(collinear_indices)
    [t0_alongwall, t1_alongwall, indices_hit_different_point] = ...
        fcn_INTERNAL_scaleTforCollinear(r(collinear_indices,:), q_minus_p(collinear_indices,:), sensor_vector(1,:), flag_search_range_type);

    % Make p and r, t and u longer so that additional hit points are
    % calculated in special case of overlaps
    tValues(collinear_indices) = t0_alongwall;
    uValues(collinear_indices) = 0.5; % Forces the t values to be accepted in later step
    tValues = [tValues; t1_alongwall(indices_hit_different_point)];
    uValues = [uValues; 0.5*ones(length(indices_hit_different_point),1)];


    % Check if there are 2 hits at different locations. Sometimes the hit
    % will land at the same location over/over, and we only want the first
    % one.
    indices_twoHits = collinear_indices(indices_hit_different_point);

    if ~isempty(indices_twoHits)
        % Duplicate selected values. Note that q (sensor start) and s (sensor
        % vector) stay 1x2 vectors.
        p = [p; p(indices_twoHits,:)];
        r = [r; r(indices_twoHits,:)];
        wall_indexes = [wall_indexes; wall_indexes(indices_twoHits,:)];
        NindicesWithTwoHits = length(indices_twoHits(:,1));
        flag_isCollinear = [flag_isCollinear; ones(NindicesWithTwoHits,1)];
    end
end

%%%%
% Fix point wall cases
if any(pointWall_indices)
    Npointwalls = length(pointWall_indices);

    % For point wall cases, the t value will be 0 if there's no "wall"
    % projection, since there is no "distance" along a point. If there is a
    % wall projection, the code returns the orthogonal distance from the
    % sensor to the point as t. The orthogonal distance, as a ratio of the
    % sensor length, is found by projecting from the sensor's start (q) to
    % the point (p). This is (p-q)
    %
    % The u value will the distance from the sensor's start (q) to the
    % point (p). This is given by (p - q), since p is the wall start and q
    % is the sensor_vector_start. Note: we've already calculated q-p
    % (vector from point to sensor start), so just use the negative here: p
    % = wall_start;
    % q = sensor_vector_start;
    % r = wall_end - wall_start;
    % s = sensor_vector_end - sensor_vector_start;
    q_minus_p_singlePoint = q_minus_p(pointWall_indices,:);
    s_singlePoint = ones(Npointwalls, 1)*sensor_vector(1,:);
    s_dot_s_singlePoint = sum(s_singlePoint.^2,2);
    u_singlePoint = sum(s_singlePoint.*(-1*q_minus_p_singlePoint),2)./s_dot_s_singlePoint;
    ortho_s_singlePoint = s_singlePoint*[0 1; -1 0];
    t_singlePoint = sum(ortho_s_singlePoint.*(-1*q_minus_p_singlePoint),2)./s_dot_s_singlePoint;
    if 0 == flag_search_range_type || 1 == flag_search_range_type
        pointsNotHit = abs(t_singlePoint)>near_zero;
        t_singlePoint(pointsNotHit) = nan;
    end

    % Save t's and u's
    tValues(pointWall_indices) = t_singlePoint;
    uValues(pointWall_indices) = u_singlePoint;

    % Set r values to ortho vector. Note: use the negative value so that
    % the wall points to the sensor, which is how t is typically defined.
    r(pointWall_indices,:) = -ortho_s_singlePoint;

end

% At this point, all the u and t vectors should be filled with "stacked"
% points

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


% u values must be zero to have sensor on wall
if any(flag_isPointSensor)
    if 0 == flag_search_range_type || 2 == flag_search_range_type
        uValues(uValues~=0) = inf;
    end
end


if 0 == flag_search_range_type
    % 0: (default) the GIVEN sensor and GIVEN wall used.
    % This constrains both t (wall extent) and u (sensor extent)
    good_vector = ((zero_threshold<=tValues).*(one_threshold>=tValues).*(zero_threshold<=uValues).*(one_threshold>=uValues));
elseif 1 == flag_search_range_type
    % 1: ANY projection of the sensor used with the GIVEN wall
    good_vector = ((zero_threshold<=tValues).*(one_threshold>=tValues));
elseif 2 == flag_search_range_type
    % 2: ANY projection of the wall used with the GIVEN sensor
    % This constrains ONLY u (sensor extent)
    good_vector = ((zero_threshold<=uValues).*(one_threshold>=uValues));
elseif 3 == flag_search_range_type
    % 3: ANY projection of BOTH the wall and sensor
    % NO constraints
    good_vector = ((Inf>=uValues).*(-Inf<=tValues));
else
    warning('on','backtrace');
    warning('Expecting a flag_search_range_type as integer in range of 0 to 3, but found: %.3f',flag_search_range_type);
    error('Bad flag_search_range_type encountered');
end

% Keep only the indices that work
within_indices = find(good_vector>0);

% Keep only intersections that are within good_vector
intersections = NaN*ones(length(p(:,1)),2);
if ~isempty(within_indices)
    % Calculate the intersection point (finally)
    result = p + tValues.*r;

    % Fix values that are infinite
    if any(isinf(tValues))
        % Will get a NaN value if multiply infinity by zero. Force NaN
        % values to be zero
        infiniteIndices = isinf(tValues);
        fixedResults = result(infiniteIndices,:);
        fixedResults(isnan(fixedResults))=0;
        result(infiniteIndices,:) = fixedResults;
    end

    % Save result
    intersections(within_indices,:) = result(within_indices,:);

end

%% Apply flag_search_return_type
% 0: (default) returns results of FIRST intersection
% 1: returns distance and location of ALL intersections
[within_indices, distances_squared] = fcn_INTERNAL_selectClosestPoint(sensor_vector_start, intersections, flag_search_return_type);
if any(flag_isPointSensor)

end

if isempty(within_indices)
    distance = nan;
    location = [nan nan];
    wall_segment = nan;
    t = nan;
    u = nan;
else
    location = intersections(within_indices,:);
    wall_segment = wall_indexes(within_indices);
    p = p(within_indices,:);
    r = r(within_indices,:);
    t = tValues(within_indices);
    u = uValues(within_indices);

    % The following converts t-values into u-values. Only need to do this for colinear indices
    within_and_colinear = find(flag_isCollinear(within_indices));
    if ~isempty(within_and_colinear)
        temp = fcn_Path_convertPerA2PerB(p, ones(length(p(:,1)),1)*q, r, ones(length(p(:,1)),1)*sensor_vector(1,:), t, -1);
        u(within_and_colinear) = temp(within_and_colinear);
    end
    validDistancesSquared = distances_squared(within_indices);
    distance = validDistancesSquared.^0.5.*sign(u);
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
    nudge = sizePlot*0.006;

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

    colorsWalls = zeros(Nwalls,3);
    if 1==1
        aveWalls = (wall_start + wall_end)/2;
        % Plot all walls in different colors
        for ith_wall = 1:Nwalls
            walls = [wall_start(ith_wall,:); wall_end(ith_wall,:)];
            h_plot = plot(walls(:,1),walls(:,2),'.-','Linewidth',3,'MarkerSize',15);
            colorsWalls(ith_wall,:) = get(h_plot,'Color');

            % Label the wall
            text(aveWalls(ith_wall,1)+nudge,aveWalls(ith_wall,2),sprintf('%.0f',ith_wall),'Color',0.5*colorsWalls(ith_wall,:));

        end

    else
        % Plot all walls as black
        plot(allWalls(:,1),allWalls(:,2),'k.-','Linewidth',5);
        handle_text = text(allWalls(1,1),allWalls(1,2),'Walls');
        set(handle_text,'Color',[0 0 0]);
    end

    % Plot the sensor vector
    if ~any(flag_isPointSensor)
        quiver(sensor_vector_start(1,1),sensor_vector_start(1,2),sensor_vector(1,1),sensor_vector(1,2),0,'r','Linewidth',2,'MaxHeadSize',1);
    end
    plot(sensor_vector_start(:,1),sensor_vector_start(:,2),'r.','Markersize',20);
    plot(sensor_vector_end(:,1),sensor_vector_end(:,2),'r.','Markersize',10);

    handle_text = text(q(:,1),q(:,2),'Sensor');
    set(handle_text,'Color',[1 0 0]);

    axis_size = axis;
    y_range = axis_size(4)-axis_size(3);

    % Plot any hits in blue
    for i_result = 1:length(distance)
        plot(location(i_result,1),location(i_result,2),'bo','Markersize',10);
        handle_text = text(location(i_result,1),location(i_result,2)-0.05*y_range,sprintf('D: %.1f',distance(i_result)));
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

%% fcn_INTERNAL_crossProduct
% Calculate cross products
function result = fcn_INTERNAL_crossProduct(v,w)
result = v(:,1).*w(:,2)-v(:,2).*w(:,1);
end % Ends fcn_INTERNAL_crossProduct

%% fcn_INTERNAL_saturateRange
function valueOut = fcn_INTERNAL_saturateRange(valueIn)
valueOut = max(0,valueIn);
valueOut = min(1,valueOut);
end % Ends fcn_INTERNAL_saturateRange

%% fcn_INTERNAL_scaleTforCollinear
function [t0_alongwall, t1_alongwall, indices_hit_different_point] = fcn_INTERNAL_scaleTforCollinear(r, q_minus_p, s, flag_search_range_type)

Nwalls = length(r(:,1));

% If they are collinear, then need to find whether there is overlapping,
% e.g. there may be infinite solutions. In other words, we need to know
% when the t-value, the percentage along each wall, starts (t0) and stops
% (t1), and where the sensor projection starts (u0) and stops (u1). For a
% given wall, these values can never be larger than 1 or less than 0,
% physcially. So, depending on flag_search_range_type, we need to "cap" t0
% and t1 at these values.
%
%
% In this case, express the endpoints of the second segment (q and q + s)
% in terms of the equation of the first line segment (p + t r):
% t = (q − p) × s / (r × s) becomes:
% t0 = (q − p) · r / (r · r)
% t1 = (q + s − p) · r / (r · r) = t0 + s · r / (r · r)
%
% If the interval between t0 and t1 intersects the interval [0, 1] then the
% line segments are collinear and overlapping; otherwise they are collinear
% and disjoint.
%
% Note that if s and r point in opposite directions, then s · r < 0 and so
% the interval to be checked is [t1, t0] rather than [t0, t1].

% Calculate raw t's. As a reminder, the t coordinates say where
% the sensor "hits" in wall coordinates. When calculating the
% intersections, we ONLY use the t-values.
r_dot_r = sum(r.*r,2);
q_minus_p_dot_r = sum(q_minus_p.*r,2);
s_dot_r = sum(s.*r,2);
t0_alongwall = q_minus_p_dot_r./r_dot_r;
t1_alongwall = t0_alongwall + s_dot_r./r_dot_r;


%%%%%
% The following steps "saturate" the t values to the interval of
% 0 and 1. As a reminder, flag_search_range_type means:
%
% 0: (default) the GIVEN sensor and GIVEN wall used.
% 1: ANY projection of the sensor used with the GIVEN wall
% 2: ANY projection of the wall used with the GIVEN sensor
% 3: ANY projection of BOTH the wall and sensor

% Are we not in "any projection of the wall mode"?
if flag_search_range_type==0
    % In this case, we keep only the t values that are both in the
    % [0,1] interval and overlap the sensor

    % Check whether there is overlap by seeing of the t0 and t1 values are
    % within the interval of [0 1], endpoint inclusive. The meaning of the
    % variables is as follows:
    % t0_inside means the start of the sensor is "inside" a wall.
    % t1_inside means the end of the sensor is "inside" a wall.
    % t0_t1_surround means the sensor starts before a wall and ends after a wall.
    t0_inside = (t0_alongwall>=0)&(t0_alongwall<=1);
    t1_inside = (t1_alongwall>=0)&(t1_alongwall<=1);
    t0_t1_surround = (t0_alongwall<0)&(t1_alongwall>1) | (t1_alongwall<0)&(t0_alongwall>1);
    any_within = t0_inside | t1_inside | t0_t1_surround;
    within_indices = find(any_within);

    % Fix the ranges to be within 0 and 1
    t0_alongwall(within_indices) = fcn_INTERNAL_saturateRange(t0_alongwall(within_indices));
    t1_alongwall(within_indices) = fcn_INTERNAL_saturateRange(t1_alongwall(within_indices));

elseif flag_search_range_type==1
    % This is the case where the sensor projects, but the wall stays
    % where it is at. So the sensor ALWAYS hits the wall at t=0 and
    % t=1, because it always overlaps the wall
    t0_alongwall = zeros(Nwalls,1);
    t1_alongwall = ones(Nwalls,1);

elseif flag_search_range_type==2
    % This is the case where the wall projects, but the sensor stays
    % where it is at. So the sensor ALWAYS hits the wall the locations
    % where the sensor is located, e.g. t0_alongwall

    % Do nothing. The previously calculated t0 and t1 values are the
    % correct ones!

elseif flag_search_range_type==3
    % This is the case where the wall and the sensor both project. The
    % intersection points are at -inf and inf;
    t0_alongwall = -inf(Nwalls,1);
    t1_alongwall = inf(Nwalls,1);

else
    warning('on','backtrace');
    warning('Expecting a flag_search_range_type as integer in range of 0 to 3, but found: %.3f',flag_search_range_type);
    error('Bad flag_search_range_type encountered');
end


%%%%%
% Allow multiple intersections because of overlap. At this point, we do
% not know which intersection to keep since there may be many, so we
% have to generate all possibilities and search through them later. Do
% we need to add more hit points
indices_hit_different_point = find(t0_alongwall~=t1_alongwall);

end % Ends fcn_INTERNAL_scaleTforCollinear


%% fcn_INTERNAL_selectClosestPoint
function [within_indices, distances_squared] = fcn_INTERNAL_selectClosestPoint(sensor_vector_start, intersections, flag_search_return_type)

% Find the distances via Euclidian distance to the sensor's origin
% note: a faster way to do this might be to just
% calculate t*r as a length
distances_squared = sum((intersections - sensor_vector_start).^2,2);
within_indices = find(~isnan(distances_squared));
if ~isempty(within_indices)
    if 0==flag_search_return_type
        % Keep only the minimum distance result
        [~,within_indices] = min(distances_squared);

    elseif 1==flag_search_return_type
        % Return all the results by default
    else
        warning('on','backtrace');
        warning('Expecting a flag_search_return_type as integer with values of 0 or 1, but found: %.3f',flag_search_return_type);
        error('Bad flag_search_range_type encountered');
    end
end

end % Ends fcn_INTERNAL_selectClosestPoint