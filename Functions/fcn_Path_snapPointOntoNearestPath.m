function [closest_path_points,...
    s_coordinates,...
    first_path_point_index,...
    second_path_point_index,...
    percent_along_length,...
    distance_real,...
    distance_imaginary] = ...
    fcn_Path_snapPointOntoNearestPath(points, path, varargin)
% fcn_Path_snapPointOntoNearestPath
% Finds location on a path that is closest to a given query point, e.g.
% snapping the point onto the path.
%
% The solution method is as follows:
%  1. Find the closest point on the path to the query point
%  2. Find vectors behind and ahead of the closest point from step 1
%   -> Check for end cases. For start/end, use the adjacent points.
%  3. Find percentage of travel on both the segments using dot products
%  4. Find projected point on the segment that has percentage of travel
%  between 0 and 1
%  5. If percentage of travel is not between 0 and 1 for both the segments,
%  then choose the snap point as closest point, e.g. the result from step.
%  For this special case, if the point is not clearly orthogonal to the
%  path (which can happen at endpoints for flag types other than 4), it
%  calculates the imaginary portion which is the orthogonal distance to the
%  point as measured from the unit orthogonal projection from the path.
% 
% FORMAT: 
%
%    [closest_path_points,...
%     s_coordinate,...
%     first_path_point_index,...
%     second_path_point_index,...
%     percent_along_length,...
%     distance_real,...
%     distance_imaginary] = ...
%      fcn_Path_snapPointOntoNearestPath(points, path,{fig_num})
%
% INPUTS:
%
%      points: a Mx2 vector containing the [X Y] location of the point or
%              a Mx3 vector containing the [X Y Z] location of the point
%              where M is the number of rows of different points
%
%      path: a Nx2 or Nx3 vector of [X Y (Z)] path points, where N is the
%      number of points the points on the path, with N >= 2.
%
%      (OPTIONAL INPUTS)
%      flag_rounding_type: a flag to indicate which type of projection is
%      used, especially when stations are located at the end-points of
%      segments within the nearby_traversal. When stations are at the
%      end-points of segments, the normal vector is undefined as it depends
%      on whether to use the prior or subsequent segment, or some
%      combination of these.
%
%      Note that the very first point always uses projections from the
%      following segement, and the very last point always uses the prior.
%      Otherwise, the flag determines behaviors for endpoints of internal
%      segments. The options include:
%
%          flag_rounding_type = 1;  % This is the default, and indicates
%          that the orthogonal projection of an endpoint is created by the
%          PRIOR segment leading up to each station query point.
% 
%          flag_rounding_type = 2;  % This indicates that the orthogonal
%          projection of an endpoint is created by the FOLLOWING segment
%          after each station query point.
% 
%          flag_rounding_type = 3;  % This indicates that the orthogonal
%          projection, ONLY if the station query falls at the joining point
%          between two segments (e.g. is on the "joint"), then the
%          projection is created by averaging the vector projections
%          created from the PRIOR segment and FOLLOWING segment.
% 
%          flag_rounding_type = 4;  % This indicates that the orthogonal
%          projections along segments should be calculated at the midpoints
%          of each segment, and then for each station qeuary, the vector
%          projections are interpolated from their prior and subsequent
%          vectors. For the endpoints - the start and end - the vectors are
%          aligned with the endpoints in the negative and positive
%          directions, respectively.
%
%      figure_number: figure number where results are plotted
%
% OUTPUTS:
%
%      closest_path_points: a Mx2 or Mx3 vector containing the [X Y] or [X Y
%      Z] location of the nearest point on the path
%
%      s_coordinate: a scalar (Mx1) representing the s-coordinate distance
%      along the path. Negative values indicate distances "before" the path
%      starts, measured from the projection of the first path segment.
%      Similarly, values larger than the total path length represent
%      distances "after" the path ends, measured from the projection of the
%      last path segment.
%
%      first_path_point_index: a scalar (Mx1) representing the index of the
%      starting "row" of the path where the closest point landed.
%
%      second_path_point_index: a scalar (Mx1) representing the index of the
%      ending "row" of the path where the closest point landed.
%
%      percent_along_length: a scalar (Mx1) of the percent along the
%      segment of the path where the closest point landed. Values less than
%      zero represent locations before a path starts, and values above 1
%      represent locations after the path ends.
%
% EXAMPLES:
%      
% See the script: script_test_fcn_Path_snapPointOntoNearestPath
% for a full test suite.
%
% DEPENDENCIES:
%
%     fcn_Path_checkInputsToFunctions
%     fcn_Path_convertPathToTraversalStructure
%
% This function was written on 2020_10_10 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%     2020_10_10 - first write of the code
%     2020_11_10 - changed function names in prep for DataClean class
%     2020_12_03 - updated some of the plotting/debug details to improve
%     2021_01_08 
%     -- cleaned up comments
%     -- added argument checking
%     2021_01_09
%     -- cleaned up header
%     2021_03_21
%     -- modified to allow 3D snapping (!!!)
%     -- changed input checks to include 3D paths
%     -- updated plotting to allow 3D
%     2021_12_10
%     -- updated header for clarity
%     2023_04_24
%     -- better comments
%     -- fixed a bug where it was snapping to wrong locations for acute
%     angles
%     2023_04_30
%     -- added real and imaginary distance outputs
%     2023_05_02
%     -- added flag options for snap options
%     2023_06_01 to 2023_06_05 by sbrennan@psu.edu
%     -- vectorized code
%     2023_08_08 by sbrennan@psu.edu
%     -- cleaned up end case of vector cacluations, bug fix
%     -- cleaned up plotting flags
%     -- improved the comments


flag_do_debug = 1; % Flag to plot the results for debugging
flag_do_plots = 0;
flag_check_inputs = 1; % Flag to perform input checking

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

if flag_check_inputs == 1
    % Are there the right number of inputs?
    narginchk(2,4);
    
    % Check the data input
    fcn_Path_checkInputsToFunctions(path, 'path2or3D');
    
    % Check that the dimension of the point and path match
    if length(points(1,:)) ~= length(path(1,:))
        error('The dimension of the query point, in number of columns, must match the dimension of the path, in number of columnts');
    end
end


% Does user want to specify the rounding type?
flag_rounding_type = 1;
if 3 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        flag_rounding_type=temp;
    end
end

% Does user want to show the plots?
if 4 == nargin
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end

if flag_do_debug
    fig_debug = 888;
end


% % Does user want distances? This is computationally heavy, so skip this if
% % not needed
% flag_calculate_distances = 0;
% flag_calculate_imaginary_distances = 0;
% if nargout>5
%     flag_calculate_distances = 1;
%     if nargout>6
%         flag_calculate_imaginary_distances = 1;
%     end
% end


%% Find the closest point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% How many points are on the path? How many points are we testing?
Npath = length(path(:,1));
Npoints = length(points(:,1));

% The solution method is as follows:
%  1. Find the closest point on the path to the query point
%  2. Find vectors behind and ahead of the closest point from step 1
%     2a. Check for end cases. For start/end, use the adjacent points.
%     2b. For other cases, just use ends
%  3. Find percentage of travel on both the segments using dot products
%  4. Find projected point on the segment that has percentage of travel
%  between 0 and 1
%  5. If percentage of travel is not between 0 and 1 for both the segments,
%  then choose the projected point as closest point from step 1

% Find the differences, then square X, Y, (Z) terms, and take square root
% to find lengths
path_segment_lengths = sum(diff(path,1,1).^2,2).^0.5;
% Add them up to generate stations
path_stations = [0; cumsum(path_segment_lengths)];
% Make sure last point has a "length"
path_segment_lengths = [path_segment_lengths; path_segment_lengths(end,:)];


%% Step 1: find distance from a point to every point on the path
closest_path_point_indicies = fcn_Path_findNearestPathPoints(points, path);

% Plot everything up to now?
if flag_do_debug
    figure(fig_debug)
    clf;
    hold on;
    grid on;
    axis equal

    % Plot the path
    plot(path(:,1),path(:,2),'k.-','LineWidth',3,'MarkerSize',10);

    % Plot the test points
    plot(points(:,1),points(:,2),'g.','LineWidth',3,'MarkerSize',20);
    
    temp = axis;
    axis_range_x = temp(2)-temp(1);
    axis_range_y = temp(4)-temp(3);
    percent_larger = 0.3
    axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);

    % Plot the closest points
     plot(...
         [path(closest_path_point_indicies,1)],...
         [path(closest_path_point_indicies,2)],...
         'go','MarkerSize',20);
    for ith_point = 1:Npoints
        plot(...
            [path(closest_path_point_indicies(ith_point),1) points(ith_point,1)],...
            [path(closest_path_point_indicies(ith_point),2) points(ith_point,2)],...
            'g--');
    end
end

%% Step 2.
% Find vectors behind and ahead of each of the index points in the
% path. Notes on this: Each point has associated with it two vectors:
% ahead, and behind The ahead vector for a point is defined by the line
% segment of the point, and the point after it. The behind vector for a
% point is defined by the line segment of the prior point leading up to the
% point. For the first, the behind vector is not defined and so is just a
% repeat of the "ahead" vector. For the last point, the
% ahead vector is not defined, and so is just a repeat of the "behind"
% vector. Nearly all the vectors overlap, so we only calculate them once
% from the path. The first path vector is from points 1 to 2, the second is
% points 2 to 3, etc. Tangential vectors are vectors along the directions
% of the segments Orthogonal vectors are ones that are 90 degrees, in
% positive cross-product direction, from the tangent vectors.
%
% These normal and tangential unit vectors are used to perform dot products
% to determine distances along a segment, or orthogonal to the segment.
path_vectors  = path(2:end,:) - path(1:end-1,:);

% Fill in the path vector at the very end as a repeat of the previous one
path_vectors(end+1,:)  = path_vectors(end,:);

% Calculate the unit vectors
unit_tangential_vectors = path_vectors./path_segment_lengths;
unit_orthogonal_vectors = unit_tangential_vectors*[0 1; -1 0]; % Rotate by 90 degrees



% Fill in arrays for the back and front. For example, the "back" tangent
% vector for point 3 will be the unit vector connecting 2 to 3, e.g. the
% 2nd unit vector. Similary, the "front" tangent vector at point 3 will be
% the vector connecting 3 to 4, or the 3rd unit vector.
back_indices = (0:Npath-1)';
front_indices = (1:Npath)';
back_indices(1) = 1;

% Fill in the back and front orthogonal vectors
behind_unit_tangential_vectors = unit_tangential_vectors(back_indices,:);
ahead_unit_tangential_vectors  = unit_tangential_vectors(front_indices,:);

back_unit_orthogonal_vector   = unit_orthogonal_vectors(back_indices,:);
front_unit_orthogonal_vector  = unit_orthogonal_vectors(front_indices,:);

% Show the vectors?
if flag_do_debug

    figure(fig_debug)
    
    % Plot the vectors for ahead and behind
    for ith_point = 1:Npath
        quiver(...
            path(ith_point,1)-behind_unit_tangential_vectors(ith_point,1),...
            path(ith_point,2)-behind_unit_tangential_vectors(ith_point,2),...
            behind_unit_tangential_vectors(ith_point,1),...
            behind_unit_tangential_vectors(ith_point,2),...
            'off','filled',...
            'r','LineWidth',3);
        quiver(...
            path(ith_point,1),...
            path(ith_point,2),...
            ahead_unit_tangential_vectors(ith_point,1),...
            ahead_unit_tangential_vectors(ith_point,2),...
            'off','filled',...
            'g','LineWidth',3);
    end
    plot(path(:,1),path(:,2),'k.','MarkerSize',20);
end

%% STEP 3. Find percentage of travel on both the segments using dot products


% Find the measurements to the back segment
back_start_to_query_vectors = points-path(back_indices,:);
back_projection_distances  = dot(behind_unit_tangential_vectors,back_start_to_query_vectors); % Do dot product
back_percent_along_lengths = back_projection_distances./path_segment_lengths(back_indices);
back_closest_path_points = path(back_indices,:) + path_vectors(back_indices,:).*back_percent_along_lengths;
back_s_coordinates       = path_stations(back_indices) + path_segment_lengths(back_indices).*back_percent_along_lengths;

% Find the measurements to the front segment
front_start_to_query_vectors = points-path(front_indices,:);
front_projection_distances  = dot(ahead_unit_tangential_vectors,front_start_to_query_vectors); % Do dot product
front_percent_along_lengths = front_projection_distances./path_segment_lengths(front_indices);
front_closest_path_points = path(front_indices,:) + path_vectors(front_indices,:).*front_percent_along_lengths;
front_s_coordinates       = path_stations(front_indices) + path_segment_lengths(front_indices).*front_percent_along_lengths;




%%  STEPS 4 and 5:
%  4. Find projected point on the segment that has percentage of travel
%  between 0 and 1
%  5. If percentage of travel is not between 0 and 1 for both the segments,
%  then choose the snap point as closest point, e.g. the result from step.
%  For this special case, if the point is not clearly orthogonal to the
%  path (which can happen at endpoints for flag types other than 4), it
%  calculates the imaginary portion which is the orthogonal distance to the
%  point as measured from the unit orthogonal projection from the path.

% Fill in output defaults
first_path_point_index  = closest_path_point_indicies;
second_path_point_index = closest_path_point_indicies;
percent_along_length = 0*closest_path_point_indicies;
closest_path_points = path(closest_path_point_indicies,:);
s_coordinates = 0*closest_path_point_indicies;

%% Check for the case of points before the start
indicies_to_check = find(back_percent_along_lengths<0, 1);
if ~isempty(indicies_to_check)
    % Point must be located BEHIND the vector that is the rear-most vector -
    % This typically happens for points before the start point. 

    % Check that the nearest point is indeed the first point
    if ~all(closest_path_point_indicies(indicies_to_check)==1)
        error('Negative percent distances detected when closest point NOT first point. This needs to be debugged!');
    end

    % Calculate the outputs
    first_path_point_index(indicies_to_check)  = 1; % closest_path_point_indicies(indicies_to_check);
    second_path_point_index(indicies_to_check) = 2; % closest_path_point_indicies(indicies_to_check);
    percent_along_length(indicies_to_check)    = back_percent_along_lengths(indicies_to_check);
    closest_path_points(indicies_to_check,:)   = path(1,:);
    s_coordinates(indicies_to_check)           = back_projection_distances(1).*percent_along_length(indicies_to_check);

    % The conversion into imaginary distance depends on the snap
    % type
    %          flag_rounding_type = 1;  % This is the default, and indicates
    %          that the orthogonal projection of an endpoint is created by the
    %          PRIOR segment leading up to each station query point.
    %
    %          flag_rounding_type = 2;  % This indicates that the orthogonal
    %          projection of an endpoint is created by the FOLLOWING segment
    %          after each station query point.
    %
    %          flag_rounding_type = 3;  % This indicates that the orthogonal
    %          projection, ONLY if the station query falls at the joining point
    %          between two segments (e.g. is on the "joint"), then the
    %          projection is created by averaging the vector projections
    %          created from the PRIOR segment and FOLLOWING segment.
    %
    %          flag_rounding_type = 4;  % This indicates that the orthogonal
    %          projections along segments should be calculated at the midpoints
    %          of each segment, and then for each station qeuary, the vector
    %          projections are interpolated from their prior and subsequent
    %          vectors.
    %
    switch flag_rounding_type
        case 1 % Use the rear vector
            unit_orthogonal_projection_vector(indicies_to_check,:) = back_unit_orthogonal_vector(indicies_to_check,:);
        case 2 % Use the front vector (WHICH IN THIS CASE IS BACK)
            unit_orthogonal_projection_vector(indicies_to_check,:) = back_unit_orthogonal_vector(indicies_to_check,:);
        case 3 % Use the average of the front and rear vectors
            unit_orthogonal_projection_vector(indicies_to_check,:) = back_unit_orthogonal_vector(indicies_to_check,:);
        case 4
            % This is a very special case, handled separately
            error('Not coded yet!');
        otherwise
            error('Unknown flag_rounding_type');
    end
end

%% Check for points past the end
indicies_to_check = find(front_percent_along_lengths>1, 1); 
if ~isempty(indicies_to_check)
    % Point must be located AHEAD the vector that is the front-most vector -
    % This typically happens for points after the end point. 

    % Check that the nearest point is indeed the end point
    if ~all(closest_path_point_indicies(indicies_to_check)==Npath)
        error('Greater than 100 percent distances detected when closest point NOT end point. This needs to be debugged!');
    end

    % Calculate the outputs
    first_path_point_index(indicies_to_check)  = Npath-1; % closest_path_point_indicies(indicies_to_check);
    second_path_point_index(indicies_to_check) = Npath; % closest_path_point_indicies(indicies_to_check);
    percent_along_length(indicies_to_check)    = front_percent_along_lengths(indicies_to_check);
    closest_path_points(indicies_to_check,:)   = path(Npath,:);
    s_coordinates(indicies_to_check)           = front_projection_distances(end).*percent_along_length(indicies_to_check);

    % The conversion into imaginary distance depends on the snap
    % type
    %          flag_rounding_type = 1;  % This is the default, and indicates
    %          that the orthogonal projection of an endpoint is created by the
    %          PRIOR segment leading up to each station query point.
    %
    %          flag_rounding_type = 2;  % This indicates that the orthogonal
    %          projection of an endpoint is created by the FOLLOWING segment
    %          after each station query point.
    %
    %          flag_rounding_type = 3;  % This indicates that the orthogonal
    %          projection, ONLY if the station query falls at the joining point
    %          between two segments (e.g. is on the "joint"), then the
    %          projection is created by averaging the vector projections
    %          created from the PRIOR segment and FOLLOWING segment.
    %
    %          flag_rounding_type = 4;  % This indicates that the orthogonal
    %          projections along segments should be calculated at the midpoints
    %          of each segment, and then for each station qeuary, the vector
    %          projections are interpolated from their prior and subsequent
    %          vectors.
    %
    switch flag_rounding_type
        case 1 % Use the rear vector
            unit_orthogonal_projection_vector(indicies_to_check,:) = back_unit_orthogonal_vector(indicies_to_check,:);
        case 2 % Use the front vector (WHICH IN THIS CASE IS BACK)
            unit_orthogonal_projection_vector(indicies_to_check,:) = back_unit_orthogonal_vector(indicies_to_check,:);
        case 3 % Use the average of the front and rear vectors
            unit_orthogonal_projection_vector(indicies_to_check,:) = back_unit_orthogonal_vector(indicies_to_check,:);
        case 4
            % This is a very special case, handled separately
            error('Not coded yet!');
        otherwise
            error('Unknown flag_rounding_type');
    end
end

%% Check the case that the point is between back and front segments 
% e.g. back percentage along length is positive, but before the front
% percentage is negative. This only happens when the point is on the
% outside angle of an joint between segments, where the segments are
% bending away from the point.
indicies_to_check = find((back_percent_along_lengths>1).*(front_percent_along_lengths<0)); 
if ~isempty(indicies_to_check)

    % point is BEFORE start of front yet also AFTER end of back!
    % This is the special situation when points are before the
    % "front" segment, e.g. the next path segment AND points are
    % after the previous "back" path segment. In this special case,
    % the distance calculation just uses the snap point, and both the first
    % and second indicies are set to the same values. Because the point is
    % not actually on either segment, the percent along length is zero.
    
    % Calculate the outputs
    first_path_point_index(indicies_to_check)  = closest_path_point_indicies(indicies_to_check);
    second_path_point_index(indicies_to_check) = closest_path_point_indicies(indicies_to_check);
    percent_along_length(indicies_to_check)    = 0;
    closest_path_points(indicies_to_check,:)   = path(closest_path_point_indicies((indicies_to_check),:));
    s_coordinates(indicies_to_check)           = path_stations(closest_path_point_indicies(indicies_to_check),1);
    
    % The conversion into imaginary distance depends on the snap
    % type
    %          flag_rounding_type = 1;  % This is the default, and indicates
    %          that the orthogonal projection of an endpoint is created by the
    %          PRIOR segment leading up to each station query point.
    %
    %          flag_rounding_type = 2;  % This indicates that the orthogonal
    %          projection of an endpoint is created by the FOLLOWING segment
    %          after each station query point.
    %
    %          flag_rounding_type = 3;  % This indicates that the orthogonal
    %          projection, ONLY if the station query falls at the joining point
    %          between two segments (e.g. is on the "joint"), then the
    %          projection is created by averaging the vector projections
    %          created from the PRIOR segment and FOLLOWING segment.
    %
    %          flag_rounding_type = 4;  % This indicates that the orthogonal
    %          projections along segments should be calculated at the midpoints
    %          of each segment, and then for each station qeuary, the vector
    %          projections are interpolated from their prior and subsequent
    %          vectors.
    %
    switch flag_rounding_type
        case 1 % Use the rear vector
            unit_orthogonal_projection_vector(indicies_to_check,:) = back_unit_orthogonal_vector(indicies_to_check,:);
        case 2 % Use the front vector
            unit_orthogonal_projection_vector(indicies_to_check,:) = front_unit_orthogonal_vector(indicies_to_check,:);
        case 3 % Use the average of the front and rear vectors
            mixed_vector = back_unit_orthogonal_vector(indicies_to_check,:)+front_unit_orthogonal_vector(indicies_to_check,:);
            magnitude_mixed_vector = sum(mixed_vector.^2,2).^0.5;
            unit_mixed_vector = mixed_vector/magnitude_mixed_vector;
            unit_orthogonal_projection_vector(indicies_to_check,:) = unit_mixed_vector;
        case 4
            % This is a very special case, handled separately
            error('Not coded yet!');
        otherwise
            error('Unknown flag_rounding_type');
    end
end

%% Check the case that the point is within the "front" segment, 
% e.g. back percentage along length is positive, and front percentage is
% positive.
indicies_to_check = find((back_percent_along_lengths>1).*(front_percent_along_lengths>=0)); 
if ~isempty(indicies_to_check)
    % Only way to enter here is if point is after the back segement,
    % and exclusively in the front vector area. This occurs when the point
    % is on the front vectors "main" area. We just use the front segment
    % result then.

    % Calculate the outputs
    first_path_point_index(indicies_to_check,:)   = closest_path_point_indicies(indicies_to_check,:) ;
    second_path_point_index(indicies_to_check,:)  = closest_path_point_indicies(indicies_to_check,:) +1;
    percent_along_length(indicies_to_check,:)     = front_percent_along_lengths(indicies_to_check,:) ;
    closest_path_points(indicies_to_check,:)      = front_closest_path_points(indicies_to_check,:) ;
    s_coordinates(indicies_to_check,:)            = front_s_coordinates(indicies_to_check,:) ;
    unit_orthogonal_projection_vector(indicies_to_check,:)  = front_unit_orthogonal_vector(indicies_to_check,:) ;
end


%% Check the case that the point is purely within the "back" segment, 
% e.g. back percentage along length is positive, and front percentage is
% negative.
indicies_to_check = find((back_percent_along_lengths>=0).*(back_percent_along_lengths<=1).*(front_percent_along_lengths<0)); 
if ~isempty(indicies_to_check)
    % Point is in back segment, and before front segment. Just use
    % the back segment results.

    % Calculate the outputs
    first_path_point_index(indicies_to_check,:)  = closest_path_point_indicies(indicies_to_check,:);
    second_path_point_index(indicies_to_check,:) = closest_path_point_indicies(indicies_to_check,:)+1;
    percent_along_length(indicies_to_check,:)    = back_percent_along_lengths(indicies_to_check,:);
    closest_path_points(indicies_to_check,:)     = back_closest_path_points(indicies_to_check,:);
    s_coordinates(indicies_to_check,:)           = back_s_coordinates(indicies_to_check,:);
    unit_orthogonal_projection_vector(indicies_to_check,:) = back_unit_orthogonal_vector(indicies_to_check,:);
end

%% Check the case that the point is within the "back" and "front" segments
% e.g. back percentage along length is positive, and also within the front
% segment, e.g. the front percentage is positive. This occurs when points
% are within the interior bend of two segments joining at an angle.
indicies_to_check = find((back_percent_along_lengths>=0).*(back_percent_along_lengths<=1).*(front_percent_along_lengths>=0)); 
if ~isempty(indicies_to_check)
    % TESTED:
    % Only way to enter here is if point is in the back segement,
    % and in front segment. This occurs when the point is
    % on the both the back and front vectors "main" area. We just
    % have to use the point that is closest to the test point.
    distance_squared_front_to_point = sum((points-front_closest_path_points).^2,2);
    distance_squared_back_to_point  = sum((points-back_closest_path_points).^2,2);

    % Calculate the outputs
    front_closer_indicies = indicies_to_check(distance_squared_front_to_point<=distance_squared_back_to_point);
    rear_closer_indicies  = indicies_to_check(distance_squared_front_to_point>distance_squared_back_to_point);

    % Front segment is closer
    if ~isempty(front_closer_indicies)
        first_path_point_index(front_closer_indicies,:)  = closest_path_point_indicies(front_closer_indicies,:);
        second_path_point_index(front_closer_indicies,:) = closest_path_point_indicies(front_closer_indicies,:) + 1;
        percent_along_length(front_closer_indicies,:)    = front_percent_along_lengths(front_closer_indicies,:) ;
        closest_path_points(front_closer_indicies,:)     = front_closest_path_points(front_closer_indicies,:) ;
        s_coordinates(front_closer_indicies,:)           = front_s_coordinates(front_closer_indicies,:) ;
        unit_orthogonal_projection_vector(front_closer_indicies,:)  = front_unit_orthogonal_vector(front_closer_indicies,:) ;
    end

    % Back segment is closer
    if ~isempty(rear_closer_indicies)
        first_path_point_index(rear_closer_indicies,:)   = closest_path_point_indicies(rear_closer_indicies,:) - 1;
        second_path_point_index(rear_closer_indicies,:)  = closest_path_point_indicies(rear_closer_indicies,:);
        percent_along_length(rear_closer_indicies,:)     = back_percent_along_lengths(rear_closer_indicies,:) ;
        closest_path_points(rear_closer_indicies,:)      = back_closest_path_points(rear_closer_indicies,:) ;
        s_coordinates(rear_closer_indicies,:)            = back_s_coordinates(rear_closer_indicies,:) ;
        unit_orthogonal_projection_vector(rear_closer_indicies,:)  = back_unit_orthogonal_vector(rear_closer_indicies,:);
    end
end

% Calculate real and imaginary portions
reference_points = closest_path_points;

off_segment_indicies = find((percent_along_length<=0) || (percent_along_length>=1));
if ~isempty(off_segment_indicies)
    reference_points(off_segment_indicies,:) = path(closest_path_point_indicies(off_segment_indicies,1),:);
end
[distance_real, distance_imaginary] = ...
        fcn_INTERNAL_convertOffsetCoordinateToImaginaryNumber(unit_orthogonal_projection_vector,reference_points, points);



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
    figure(fig_num);
    clf;
    hold on;
    grid on;
    
    % Is this a 2-D query
    if length(path(1,:))==2
        % Plot the path
        plot(path(:,1),path(:,2),'r-','Linewidth',5);
        plot(path(:,1),path(:,2),'r.','Markersize',20);
        for ith_point = 1:length(path(:,1))
            text(path(ith_point,1),path(ith_point,2),sprintf('%.0d',ith_point),'Color',[1 0 0],'FontSize',12,'VerticalAlignment','bottom');
        end
        
        axis equal;
        
        % Plot the query points
        plot(points(:,1),points(:,2),'k.','MarkerSize',20);

        
        %         % Plot the closest path points;
        %         plot(...
        %             path(first_path_point_index:second_path_point_index,1),...
        %             path(first_path_point_index:second_path_point_index,2),'m.','MarkerSize',20);
        
        % Plot the reference vector
        for ith_point = 1:length(points(:,1))
            quiver(path(closest_path_point_indicies(ith_point),1),path(closest_path_point_indicies(ith_point),2),...
                unit_orthogonal_projection_vector(ith_point,1),unit_orthogonal_projection_vector(ith_point,2),...
                0,'Color',[0.5 0.5 0.5]);
            % midpoint = (path(closest_path_point_indicies(ith_point),:) + path(closest_path_point_indicies(ith_point),:)+unit_orthogonal_projection_vector(ith_point,:))/2;
            %text(midpoint(1,1),midpoint(1,2),'Reference vector','LineWidth',2,'VerticalAlignment','bottom');
            
            % Label the points with distances
            midpoint = (closest_path_points(ith_point,:) + points(ith_point,:))/2;
            if distance_imaginary==0
                distance_string = sprintf('Between points %.0d %.0d, distance = %.2f',first_path_point_index, second_path_point_index, distance_real(ith_point));
            else
                distance_string = sprintf('Between points %.0d %.0d, distance = %.2f + %.2f i',first_path_point_index, second_path_point_index, distance_real(ith_point),distance_imaginary(ith_point));
            end
            text(midpoint(1,1),midpoint(1,2),distance_string,'VerticalAlignment','top');
            
            % Plot the closest point on path
            plot(closest_path_points(ith_point,1),closest_path_points(ith_point,2),'g.','Markersize',20);
            text(closest_path_points(ith_point,1),closest_path_points(ith_point,2),'Snap Point on Path');
            
            
            % Connect closest point on path to query point
            plot(...
                [points(ith_point,1) closest_path_points(ith_point,1)],...
                [points(ith_point,2) closest_path_points(ith_point,2)],'g-','Linewidth',2);
        end
        
    elseif length(path(1,:))==3
        
        % Plot the path
        plot3(path(:,1),path(:,2),path(:,3),'r-','Linewidth',5);
        plot3(path(:,1),path(:,2),path(:,3),'ro','Markersize',20);
        
        axis equal;
        
        % Plot the query point
        plot3(points(:,1),points(:,2),points(:,3),'ko');
        text( points(:,1),points(:,2),points(:,3),'Query point');
        
        % Plot the closest path points;
        plot3(...
            path(first_path_point_index:second_path_point_index,1),...
            path(first_path_point_index:second_path_point_index,2),...
            path(first_path_point_index:second_path_point_index,3),'r*');
        
        
        % Plot the closest point on path
        plot3(...
            closest_path_points(:,1),...
            closest_path_points(:,2),...
            closest_path_points(:,3),...
            'go','Markersize',20);
        text(...
            closest_path_points(:,1),...
            closest_path_points(:,2),...
            closest_path_points(:,3),...
            'Snap Point on Path');
        
        
        % Connect closest point on path to query point
        plot3(...
            [points(:,1) closest_path_points(:,1)],...
            [points(:,2) closest_path_points(:,2)],...
            [points(:,3) closest_path_points(:,3)],...
            'g-','Linewidth',2);
        

    end
end % Ends the flag_do_debug if statement



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

%% fcn_INTERNAL_calculateVectorMeasures    
function     [closest_path_points,s_coordinates,percent_along_lengths,unit_orthogonal_projection_vectors] = ...
    fcn_INTERNAL_calculateVectorMeasures(points,...
    segment_end_points,...
    segment_start_points,...
    segment_start_stations)
% Do the dot products - define the vectors first
% See: https://mathinsight.org/dot_product for explanation Basically,
% we are seeing what amount the point_vector points in the direction of
% the path_vector.


path_vectors  = segment_end_points - segment_start_points;
path_segment_lengths  = sum(path_vectors.^2,2).^0.5;
unit_orthogonal_projection_vectors = (path_vectors./path_segment_lengths)*[0 1; -1 0]; % Rotate by 90 degrees

start_to_point_vectors = points-segment_start_points;
projection_distances  = dot(path_vectors,start_to_point_vectors)./path_segment_lengths; % Do dot product
% offset_distance  = dot(ortho_path_vector,start_to_point_vector)/path_segment_length; % Do dot product
percent_along_lengths = projection_distances./path_segment_lengths;

% Calculate the remaining outputs
closest_path_points = segment_start_points + path_vectors.*percent_along_lengths;
s_coordinates       = segment_start_stations + path_segment_lengths.*percent_along_lengths;
end % Ends fcn_INTERNAL_calculateVectorMeasures


%% fcn_INTERNAL_convertOffsetCoordinateToImaginaryNumber
function     [real_distance, imag_distance] = ...
    fcn_INTERNAL_convertOffsetCoordinateToImaginaryNumber(...
    unit_projection_vector_of_real_axis,...
    origin_point, ...
    point_to_convert)

% This function is used to convert a point's location into an imaginary
% distance, where the real portion of the distance represents the distance in
% the direction of a given unit vector, and the imaginary portion of the
% distance represents the orthogonal component. The conversion assumes a
% coordinate orientation of the Re/Im plane such that the real axis is
% aligned with the given unit projection vector representing the real axis.
%
% The reason for this function is that, when snapping a point onto a
% reference traversal, there are situations - particularly before the start
% of the traversal, after the end of the traversal, and along sharp corners
% of the traversal, where the point cannot be correctly represented solely
% by the orthogonal projection distance. For example, in sharp corners, the
% orthogonal projection is undefined and must be chosen by the user. The
% conversion of orthogonal distance thus allows correct and invertible
% converstions from path offsets in XY to Sy representations, and vice
% versa.

unit_projection_vector_of_imag_axis = unit_projection_vector_of_real_axis*[0 1; -1 0]; % Rotate by 90 degrees

start_to_point_vector = point_to_convert-origin_point;
real_distance  = dot(unit_projection_vector_of_real_axis,start_to_point_vector); % Do dot product
imag_distance  = dot(unit_projection_vector_of_imag_axis,start_to_point_vector); % Do dot product

% % Measure the distance from the closest_path_point station
% distance_imaginary = distance_imaginary - path_station(closest_path_point_index,1);
% 
% % Need to flip the sign so that the projection is correctly positioned
% % relative to the origin
% distance_imaginary = -1.0*distance_imaginary;
end % Ends fcn_INTERNAL_convertOffsetCoordinateToImaginaryNumber

