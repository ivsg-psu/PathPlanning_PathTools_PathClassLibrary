function [closest_path_points,...
    s_coordinates,...
    first_path_point_indicies,...
    second_path_point_indicies,...
    percent_along_length,...
    distances_real,...
    distances_imaginary] = ...
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
%      fcn_Path_snapPointOntoNearestPath(points, path,(flag_rounding_type), (fig_num))
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
%      first_path_point_indicies: a scalar (Mx1) representing the index of the
%      starting "row" of the path where the closest point landed.
%
%      second_path_point_indicies: a scalar (Mx1) representing the index of the
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
%     2023_08_27 by sbrennan@psu.edu
%     -- fixed bug in case 3 where unit vector calculated wrong!

flag_do_debug = 0; % Flag to plot the results for debugging
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
    fig_debug = 888; %#ok<*UNRCH> 
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
% flag_3D = 0;
point_dimension = length(points(1,:));


% The solution method is as follows:
%  1. Find the closest point on the path to the query point
%  2. Find unit vectors behind and ahead of the closest point from step 1
%  3. Find distance of travel on both the segments using dot products
%  with each direction
%  4. Find projected point on the segment that has positive tangential
%  distance
%  4a. If one is positive, use that segment to calculate values
%  4b. If both are positive, use the one that has the smallest orthogonal distance,
%  4c. If neither are positive, use the flag type to disambiguate

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
    plot(path(:,1),path(:,2),'k.-','LineWidth',5,'MarkerSize',10);

    % Plot the test points
    plot(points(:,1),points(:,2),'b.','LineWidth',5,'MarkerSize',40);
    
    % Make axis slightly larger
    temp = axis;
    axis_range_x = temp(2)-temp(1);
    axis_range_y = temp(4)-temp(3);
    percent_larger = 0.3;
    axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);

    % Plot the closest points
     plot(...
         path(closest_path_point_indicies,1),...
         path(closest_path_point_indicies,2),...
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


% % These normal and tangential unit vectors are used to perform dot products
% % to determine distances along a segment, or orthogonal to the segment.
% path_vectors  = path(2:end,:) - path(1:end-1,:);
% 
% % Fill in the path vector at the very end as a repeat of the previous one
% path_vectors(end+1,:)  = path_vectors(end,:);
% 
% % Calculate the unit vectors
% unit_tangential_vectors = path_vectors./path_segment_lengths;
% if point_dimension == 3
%     rotation_matrix = [0 1 0; -1 0 0; 0 0 1];
% else
%     rotation_matrix = [0 1; -1 0];
% end
% unit_orthogonal_vectors_at_midpoints = unit_tangential_vectors*rotation_matrix; % Rotate by 90 degrees

% Call the function fcn_Path_findPathOrthogonalVectors to find ortho
% vectors
% unit_orthogonal_vectors_at_midpoints is not used, so is commented out
[unit_orthogonal_vectors_at_midpoints, unit_orthogonal_vectors_at_joints] = ...
    fcn_Path_findPathOrthogonalVectors(path,flag_rounding_type);
unit_tangential_vectors_at_midpoints = unit_orthogonal_vectors_at_midpoints*[0 -1; 1 0];
unit_tangential_vectors_at_joints = unit_orthogonal_vectors_at_joints*[0 -1; 1 0];

% Fill in arrays for the back and front. For example, the "back" tangent
% vector for point 3 will be the unit vector connecting 2 to 3, e.g. the
% 2nd unit vector. Similary, the "front" tangent vector at point 3 will be
% the vector connecting 3 to 4, or the 3rd unit vector.
back_indices = [1; (1:Npath-1)'];
front_indices = [(1:Npath-1)'; Npath-1];

% Fill in the back and front orthogonal vectors
behind_unit_tangential_vectors_at_midpoints = -unit_tangential_vectors_at_midpoints(back_indices,:);
ahead_unit_tangential_vectors_at_midpoints  = unit_tangential_vectors_at_midpoints(front_indices,:);

behind_unit_orthogonal_vectors_at_midpoints   = unit_orthogonal_vectors_at_midpoints(back_indices,:);
ahead_unit_orthogonal_vectors_at_midpoints    = unit_orthogonal_vectors_at_midpoints(front_indices,:);

% Show the vectors?
if flag_do_debug

    figure(fig_debug)
    
    % Plot the tangent vectors for ahead and behind
    for ith_point = 1:Npath
        quiver(...
            path(ith_point,1),...
            path(ith_point,2),...
            behind_unit_tangential_vectors_at_midpoints(ith_point,1),...
            behind_unit_tangential_vectors_at_midpoints(ith_point,2),...
            0,'g-','filled',...
            'LineWidth',3);
        quiver(...
            path(ith_point,1),...
            path(ith_point,2),...
            ahead_unit_tangential_vectors_at_midpoints(ith_point,1),...
            ahead_unit_tangential_vectors_at_midpoints(ith_point,2),...
            0,'c-','filled',...
            'LineWidth',3);
    end

    % Plot the orthogonal vectors for ahead and behind
    for ith_point = 1:Npath
        quiver(...
            path(ith_point,1)+behind_unit_tangential_vectors_at_midpoints(ith_point,1),...
            path(ith_point,2)+behind_unit_tangential_vectors_at_midpoints(ith_point,2),...
            behind_unit_orthogonal_vectors_at_midpoints(ith_point,1),...
            behind_unit_orthogonal_vectors_at_midpoints(ith_point,2),...
            0,'r-','filled',...
            'LineWidth',3);
        quiver(...
            path(ith_point,1)+ahead_unit_tangential_vectors_at_midpoints(ith_point,1),...
            path(ith_point,2)+ahead_unit_tangential_vectors_at_midpoints(ith_point,2),...
            ahead_unit_orthogonal_vectors_at_midpoints(ith_point,1),...
            ahead_unit_orthogonal_vectors_at_midpoints(ith_point,2),...
            0,'m-','filled',...
            'LineWidth',3);
    end

    plot(path(:,1),path(:,2),'k.','MarkerSize',20);
end

%% STEP 3. Find distance of travel on both the directions using dot products

% Find the distances to the front and back segment
query_vectors = points-path(closest_path_point_indicies,:);
query_point_back_tangential_distances   = sum(behind_unit_tangential_vectors_at_midpoints(closest_path_point_indicies,:).*query_vectors,2); % Do dot product
query_point_front_tangential_distances  = sum(ahead_unit_tangential_vectors_at_midpoints(closest_path_point_indicies,:) .*query_vectors,2); % Do dot product
query_point_back_orthogonal_distances   = sum(behind_unit_orthogonal_vectors_at_midpoints(closest_path_point_indicies,:).*query_vectors,2); % Do dot product
query_point_front_orthogonal_distances  = sum(ahead_unit_orthogonal_vectors_at_midpoints(closest_path_point_indicies,:) .*query_vectors,2); % Do dot product


% Show the results?
if flag_do_debug

    figure(fig_debug)
    
    % Plot the tangent vectors for ahead and behind
    for ith_point = 1:Npoints
        % Plot any positive measurements in rear
        if query_point_back_tangential_distances>0
            quiver(...
                path(closest_path_point_indicies(ith_point),1),...
                path(closest_path_point_indicies(ith_point),2),...
                query_point_back_tangential_distances*behind_unit_tangential_vectors_at_midpoints(closest_path_point_indicies(ith_point),1),...
                query_point_back_tangential_distances*behind_unit_tangential_vectors_at_midpoints(closest_path_point_indicies(ith_point),2),...
                0,'-','filled',...
                'Color',[1 0 1],...
                'LineWidth',3);

            quiver(...
                path(closest_path_point_indicies(ith_point),1)+query_point_back_tangential_distances*behind_unit_tangential_vectors_at_midpoints(closest_path_point_indicies(ith_point),1),...
                path(closest_path_point_indicies(ith_point),2)+query_point_back_tangential_distances*behind_unit_tangential_vectors_at_midpoints(closest_path_point_indicies(ith_point),2),...
                query_point_back_orthogonal_distances*behind_unit_orthogonal_vectors_at_midpoints(closest_path_point_indicies(ith_point),1),...
                query_point_back_orthogonal_distances*behind_unit_orthogonal_vectors_at_midpoints(closest_path_point_indicies(ith_point),2),...
                0,'-','filled',...
                'Color',[1 0 1],...
                'LineWidth',3);
            
            
        end

        % Plot any positive measurements in front
        if query_point_front_tangential_distances>0
            quiver(...
                path(closest_path_point_indicies(ith_point),1),...
                path(closest_path_point_indicies(ith_point),2),...
                query_point_front_tangential_distances*ahead_unit_tangential_vectors_at_midpoints(closest_path_point_indicies(ith_point),1),...
                query_point_front_tangential_distances*ahead_unit_tangential_vectors_at_midpoints(closest_path_point_indicies(ith_point),2),...
                0,'-','filled',...
                'Color',[1 0 1],...
                'LineWidth',3);


            quiver(...
                path(closest_path_point_indicies(ith_point),1)+query_point_front_tangential_distances*ahead_unit_tangential_vectors_at_midpoints(closest_path_point_indicies(ith_point),1),...
                path(closest_path_point_indicies(ith_point),2)+query_point_front_tangential_distances*ahead_unit_tangential_vectors_at_midpoints(closest_path_point_indicies(ith_point),2),...
                query_point_front_orthogonal_distances*ahead_unit_orthogonal_vectors_at_midpoints(closest_path_point_indicies(ith_point),1),...
                query_point_front_orthogonal_distances*ahead_unit_orthogonal_vectors_at_midpoints(closest_path_point_indicies(ith_point),2),...
                0,'-','filled',...
                'Color',[1 0 1],...
                'LineWidth',3);

        end
    end
end


%%  STEP 4
%  4. Find projected point on the segment that has positive tangential
%  distance
%  4a. If one is strictly positive, use that segment to calculate values
%  4b. If both are strictly positive, use the one that has the smallest orthogonal distance,
%  4c. If neither are strictly positive, use the flag type to disambiguate


% Fill in output defaults
closest_path_vertex_points = path(closest_path_point_indicies,:);
closest_s_coordinates = path_stations(closest_path_point_indicies);
s_coordinate_perturbations = nan*closest_s_coordinates;
unit_orthogonal_projection_vectors = nan(length(closest_s_coordinates),point_dimension);

% Keep flags to check whether we need to calculate an imaginary component
flags_calculate_imaginary = nan(length(closest_s_coordinates),1);

%% 4a. If one is strictly positive, use that segment to calculate values
% Is the point in the back segement, and definitely NOT in the
% front segment?
query_point_indicies_to_check = find((query_point_back_tangential_distances>0).*(query_point_front_tangential_distances<0)); 
if ~isempty(query_point_indicies_to_check)
    % This is the case where the point is squarely in the back section
    % TESTED
    % Point must be located BEHIND the closest point 
    s_coordinate_perturbations(query_point_indicies_to_check) = -query_point_back_tangential_distances(query_point_indicies_to_check);
    unit_orthogonal_projection_vectors(query_point_indicies_to_check,:) = behind_unit_orthogonal_vectors_at_midpoints(closest_path_point_indicies(query_point_indicies_to_check),:);
    flags_calculate_imaginary(query_point_indicies_to_check) = 0;
end

% Is the point definitely NOT in the back segement, and definitely in the
% front segment?
query_point_indicies_to_check = find((query_point_back_tangential_distances<0).*(query_point_front_tangential_distances>0)); 
if ~isempty(query_point_indicies_to_check)
    % This is the case where the point is squarely in the front section

    % TESTED
    % Point must be located AHEAD of the closest point 
    s_coordinate_perturbations(query_point_indicies_to_check) =  query_point_front_tangential_distances(query_point_indicies_to_check);
    unit_orthogonal_projection_vectors(query_point_indicies_to_check,:) = ahead_unit_orthogonal_vectors_at_midpoints(closest_path_point_indicies(query_point_indicies_to_check),:);
    flags_calculate_imaginary(query_point_indicies_to_check) = 0;
end

%%  4b. If both are positive, use the one that has the smallest orthogonal distance,
% Are both strictly positive?
query_point_indicies_to_check = find((query_point_back_tangential_distances>0).*(query_point_front_tangential_distances>0)); 
if ~isempty(query_point_indicies_to_check)
    % This is the case where the point is in the inside of a bend between front and back

    % TESTED
    % Is point closer to rear, orthogonally?
    abs_back_distances = abs(query_point_back_orthogonal_distances(query_point_indicies_to_check));
    abs_front_distances = abs(query_point_front_orthogonal_distances(query_point_indicies_to_check));
    back_is_closer = query_point_indicies_to_check(abs_back_distances<abs_front_distances);
    s_coordinate_perturbations(back_is_closer) =  -query_point_back_tangential_distances(back_is_closer);
    unit_orthogonal_projection_vectors(back_is_closer,:) = behind_unit_orthogonal_vectors_at_midpoints(closest_path_point_indicies(back_is_closer),:);

    % Is point closer to front, orthogonally?
    front_is_closer = query_point_indicies_to_check(abs_back_distances>=abs_front_distances);
    s_coordinate_perturbations(front_is_closer) =  query_point_front_tangential_distances(front_is_closer);
    unit_orthogonal_projection_vectors(front_is_closer,:) = ahead_unit_orthogonal_vectors_at_midpoints(closest_path_point_indicies(front_is_closer),:);

    flags_calculate_imaginary(query_point_indicies_to_check) = 0;
end

%%  4c. If neither are positive, it's a joint. Use the flag type to disambiguate
%query_point_indicies_to_check = find((query_point_back_tangential_distances<=0).*(query_point_front_tangential_distances<=0)); 
query_point_indicies_to_check = find(isnan(s_coordinate_perturbations));
if ~isempty(query_point_indicies_to_check)
    % This is the case where the point is in the outside of a bend between front and back

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
            unit_orthogonal_projection_vectors(query_point_indicies_to_check,:) = behind_unit_orthogonal_vectors_at_midpoints(closest_path_point_indicies(query_point_indicies_to_check),:);

        case 2 % Use the front vector 
            unit_orthogonal_projection_vectors(query_point_indicies_to_check,:) = ahead_unit_orthogonal_vectors_at_midpoints(closest_path_point_indicies(query_point_indicies_to_check),:);

        case 3 % Use the average of the front and rear vectors
            average_vector = (ahead_unit_orthogonal_vectors_at_midpoints(closest_path_point_indicies(query_point_indicies_to_check),:)+behind_unit_orthogonal_vectors_at_midpoints(closest_path_point_indicies(query_point_indicies_to_check),:))/2;
            mag_average_vector = sum(average_vector.^2,2).^0.5;
            unit_orthogonal_projection_vectors(query_point_indicies_to_check,:) = average_vector./mag_average_vector;
        
        case 4
            % This is a very special case, handled separately
            error('This method is not coded yet! Try again later, maybe?');
        otherwise
            error('Unknown flag_rounding_type');
    end
    s_coordinate_perturbations(query_point_indicies_to_check) = 0;
    flags_calculate_imaginary(query_point_indicies_to_check) = 1;

end

% Calculate the remaining outputs
% The s-coordinate will be the station at the closest point plus the
% perturbations
s_coordinates          = closest_s_coordinates + s_coordinate_perturbations;

% Check for s-coordinates less than zero or greater than the path_length
query_point_indicies_to_check = find((s_coordinates<0)+(s_coordinates>path_stations(end))); 
if ~isempty(query_point_indicies_to_check)
    flags_calculate_imaginary(query_point_indicies_to_check) = 1;
end


% Set the defaults assuming s-coordinate perturbations are positive
first_path_point_indicies  = closest_path_point_indicies;
second_path_point_indicies = min(closest_path_point_indicies+1,Npath);
percent_along_length    = s_coordinate_perturbations./path_segment_lengths(closest_path_point_indicies);

% For situations where the 1st and 2nd indicies are at the end, this
% is a weird s-coordinate situation for path percentage
second_indicies_to_check = find((s_coordinate_perturbations>0).*(first_path_point_indicies==second_path_point_indicies));
if ~isempty(second_indicies_to_check)
    percent_along_length(second_indicies_to_check) = (path_segment_lengths(first_path_point_indicies(second_indicies_to_check)) + ...
        s_coordinate_perturbations(second_indicies_to_check))./path_segment_lengths(first_path_point_indicies(second_indicies_to_check));
end

% Replace any joint situations - these are weird
query_point_indicies_to_check = find((flags_calculate_imaginary==1)); 
if ~isempty(query_point_indicies_to_check)
    second_path_point_indicies(query_point_indicies_to_check) = first_path_point_indicies(query_point_indicies_to_check);
end


% Replace any negative values
query_point_indicies_to_check = find((s_coordinate_perturbations<0)); 
if ~isempty(query_point_indicies_to_check)
    first_path_point_indicies(query_point_indicies_to_check)  = max(1,closest_path_point_indicies(query_point_indicies_to_check)-1);
    second_path_point_indicies(query_point_indicies_to_check) = closest_path_point_indicies(query_point_indicies_to_check);
    percent_along_length(query_point_indicies_to_check)    = (path_segment_lengths(first_path_point_indicies(query_point_indicies_to_check)) + ...
        s_coordinate_perturbations(query_point_indicies_to_check))./path_segment_lengths(first_path_point_indicies(query_point_indicies_to_check));

    % For situations where the 1st and 2nd indicies are at the start, this
    % is a weird s-coordinate situation for path percentage
    second_indicies_to_check = find((s_coordinate_perturbations<0).*(first_path_point_indicies==second_path_point_indicies)); 
    if ~isempty(second_indicies_to_check)
        percent_along_length(second_indicies_to_check) = s_coordinate_perturbations(second_indicies_to_check)./path_segment_lengths(closest_path_point_indicies(second_indicies_to_check));
    end
end

% Calculate real and imaginary portions
reference_points = closest_path_vertex_points;

[distances_real, distances_imaginary] = ...
        fcn_INTERNAL_convertOffsetCoordinateToImaginaryNumber(unit_orthogonal_projection_vectors,reference_points, points);
distances_imaginary = distances_imaginary.*(1.00*flags_calculate_imaginary);

% Calculate closest_path_points
closest_path_points = path(first_path_point_indicies,:) + ...
                percent_along_length.*(path(second_path_point_indicies,:)-path(first_path_point_indicies,:));

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
    axis equal;
    
    % Is this a 2-D query
    if length(path(1,:))==2
        % Plot the path
        plot(path(:,1),path(:,2),'r.-','Linewidth',5,'Markersize',20);


        % Label the path points
        for ith_point = 1:Npath
            text(path(ith_point,1),path(ith_point,2),sprintf('%.0d',ith_point),'Color',[1 0 0],'FontSize',12,'VerticalAlignment','bottom');
        end
                
        % Plot the query points
        plot(points(:,1),points(:,2),'k.','MarkerSize',20);

        % Label the query points
        for ith_point = 1:Npoints
            text(points(ith_point,1),points(ith_point,2),sprintf('%.0d',ith_point),'Color',[0 0 0],'FontSize',12,'VerticalAlignment','bottom');
        end

        
        %         % Plot the closest path points;
        %         plot(...
        %             path(first_path_point_index:second_path_point_index,1),...
        %             path(first_path_point_index:second_path_point_index,2),'m.','MarkerSize',20);
        
        % Plot the reference vector
        for ith_point = 1:Npoints
            quiver(path(closest_path_point_indicies(ith_point),1),path(closest_path_point_indicies(ith_point),2),...
                unit_orthogonal_projection_vectors(ith_point,1),unit_orthogonal_projection_vectors(ith_point,2),...
                0,'Color',[0.5 0.5 0.5]);
            % midpoint = (path(closest_path_point_indicies(ith_point),:) + path(closest_path_point_indicies(ith_point),:)+unit_orthogonal_projection_vector(ith_point,:))/2;
            %text(midpoint(1,1),midpoint(1,2),'Reference vector','LineWidth',2,'VerticalAlignment','bottom');
            
            % Label the points with distances
            midpoint = (closest_path_points(ith_point,:) + points(ith_point,:))/2;

            if distances_imaginary(ith_point)==0
                distance_string = sprintf('Between points %.0d %.0d,\n distance = %.2f,\n %% on path = %.2f',...
                    first_path_point_indicies(ith_point), second_path_point_indicies(ith_point), distances_real(ith_point), percent_along_length(ith_point));
            else
                distance_string = sprintf('Between points %.0d %.0d,\n distance = %.2f + %.2f i,\n %% on path = %.2f',...
                    first_path_point_indicies(ith_point), second_path_point_indicies(ith_point), distances_real(ith_point),distances_imaginary(ith_point),percent_along_length(ith_point));
            end

            text(midpoint(1,1),midpoint(1,2),distance_string,'VerticalAlignment','top');            


            % Connect closest point on path to query point
            closest_point = closest_path_points(ith_point,:);

            plot(...
                [points(ith_point,1) closest_point(1,1)],...
                [points(ith_point,2) closest_point(1,2)],'g-','Linewidth',2);

            % Plot the closest point on path
            plot(closest_point(1,1),closest_point(1,2),'g.','Markersize',20);
            text(closest_point(1,1),closest_point(1,2),'Snap Point on Path');
        end

        % Make axis slightly larger
        temp = axis;
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);

        
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
            path(first_path_point_indicies:second_path_point_indicies,1),...
            path(first_path_point_indicies:second_path_point_indicies,2),...
            path(first_path_point_indicies:second_path_point_indicies,3),'r*');
        
        
        % Plot the closest point on path
        plot3(...
            closest_path_vertex_points(:,1),...
            closest_path_vertex_points(:,2),...
            closest_path_vertex_points(:,3),...
            'go','Markersize',20);
        text(...
            closest_path_vertex_points(:,1),...
            closest_path_vertex_points(:,2),...
            closest_path_vertex_points(:,3),...
            'Snap Point on Path');
        
        
        % Connect closest point on path to query point
        plot3(...
            [points(:,1) closest_path_vertex_points(:,1)],...
            [points(:,2) closest_path_vertex_points(:,2)],...
            [points(:,3) closest_path_vertex_points(:,3)],...
            'g-','Linewidth',2);
        
        view(3);

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

% %% fcn_INTERNAL_calculateVectorMeasures    
% function     [closest_path_points,s_coordinates,percent_along_lengths,unit_orthogonal_projection_vectors] = ...
%     fcn_INTERNAL_calculateVectorMeasures(points,...
%     segment_end_points,...
%     segment_start_points,...
%     segment_start_stations)
% % Do the dot products - define the vectors first
% % See: https://mathinsight.org/dot_product for explanation Basically,
% % we are seeing what amount the point_vector points in the direction of
% % the path_vector.
% 
% 
% path_vectors  = segment_end_points - segment_start_points;
% path_segment_lengths  = sum(path_vectors.^2,2).^0.5;
% unit_orthogonal_projection_vectors = (path_vectors./path_segment_lengths)*[0 1; -1 0]; % Rotate by 90 degrees
% 
% start_to_point_vectors = points-segment_start_points;
% projection_distances  = dot(path_vectors,start_to_point_vectors)./path_segment_lengths; % Do dot product
% % offset_distance  = dot(ortho_path_vector,start_to_point_vector)/path_segment_length; % Do dot product
% percent_along_lengths = projection_distances./path_segment_lengths;
% 
% % Calculate the remaining outputs
% closest_path_points = segment_start_points + path_vectors.*percent_along_lengths;
% s_coordinates       = segment_start_stations + path_segment_lengths.*percent_along_lengths;
% end % Ends fcn_INTERNAL_calculateVectorMeasures


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

point_dimension = length(unit_projection_vector_of_real_axis(1,:));
if point_dimension == 3
    rotation_matrix = [0 1 0; -1 0 0; 0 0 1];
else
    rotation_matrix = [0 1; -1 0];
end

unit_projection_vector_of_imag_axis = unit_projection_vector_of_real_axis*rotation_matrix; % Rotate by 90 degrees

start_to_point_vector = point_to_convert-origin_point;
real_distance  = sum(unit_projection_vector_of_real_axis.*start_to_point_vector,2); % Do dot product
imag_distance  = sum(unit_projection_vector_of_imag_axis.*start_to_point_vector,2); % Do dot product

% % Measure the distance from the closest_path_point station
% distance_imaginary = distance_imaginary - path_station(closest_path_point_index,1);
% 
% % Need to flip the sign so that the projection is correctly positioned
% % relative to the origin
% distance_imaginary = -1.0*distance_imaginary;
end % Ends fcn_INTERNAL_convertOffsetCoordinateToImaginaryNumber

