function [closest_path_points,...
    s_coordinates,...
    first_path_point_indicies,...
    second_path_point_indicies,...
    percent_along_length,...
    distances_real,...
    distances_imaginary] = ...
    fcn_Path_snapPointToPathViaVectors(points, path, varargin)
% fcn_Path_snapPointToPathViaVectors
% Finds location on a path that is closest to a given query point, e.g.
% snapping the point onto the path. This function is similar to
% fcn_Path_snapPointOntoNearestPath and gives the same outputs; however, it
% uses a vector calculation method rather than point-distance metric. This
% improves robustness.
%
% The solution method is as follows:
%  1. Find the unit orthogonal vectors of each path segment
%  2. Project from the end of each segment to the query point(s), creating
%  two vectors per segment, one at each end-point of each segment, pointing
%  to each query point.
%  3. Take the cross product of the point-projection vectors at each end of
%  a line segment with the orthogonal projection vector for the segment.
%  4. If the cross-product signs change, then the point must be within the
%  orthogonal projection of the line segment. For these cases, calculate
%  the dot product of the orthogonal projection with one of the point
%  projection vectors to find the distance to the point. Keep the minimimum
%  distance among all segments, save the segment indicies, percent along
%  length, etc.
%  5. There are cases where a point is at the vertex of a segment, and not
%  within the projection of the line segment. To query these cases, find
%  the distance of the query point to all points in the line segment,
%  keeping the minimum. If the minimum distance to a vertex is less than
%  any of the minimum distances to segment projections, use this distance,
%  keep the vertex as the segment indicies, set percent along length to
%  zero, and calculate real/imaginary components depending on the
%  projection type.
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
%      fcn_Path_snapPointToPathViaVectors(points, path,(flag_rounding_type), (fig_num))
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
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
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
% See the script: script_test_fcn_Path_snapPointToPathViaVectors
% for a full test suite.
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%     fcn_Path_findPathOrthogonalVectors
%
% This function was written on 2023_09_28 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2023_09_28 - S. Brennan
% --  first write of the code, using fcn_Path_snapPointOntoNearestPath 
% 2024_03_14 - S. Brennan
% -- fixed bug where snap breaks if path is passed in as a 3D vector
% 2025_06_23 - S. Brennan
% -- Updated debugging and input checks

% TO-DO
% (none)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==4 && isequal(varargin{end},-1))
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
        narginchk(2,4);

        % Check the data input
        fcn_DebugTools_checkInputsToFunctions(path, 'path2or3D');

        % Check that the dimension of the point and path match
        if length(points(1,:)) ~= length(path(1,:))
            error('The dimension of the query point, in number of columns, must match the dimension of the path, in number of columnts');
        end
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
flag_do_plots = 0; % Default is to NOT show plots
if (0==flag_max_speed) && (4 == nargin) 
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        fig_num = temp;
        figure(fig_num);
        flag_do_plots = 1;
    end
else
    if flag_do_debug
        fig = figure;  
        fig_num = fig.Number;
        flag_do_plots = 1;
    end
end


if flag_do_debug
    fig_debug = 888; %#ok<*UNRCH>
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot all the inputs
if flag_do_debug
    figure(fig_debug)
    clf;
    hold on;
    grid on;
    axis equal

    % Plot the path
    plot(path(:,1),path(:,2),'k.-','LineWidth',5,'MarkerSize',20);


    % Plot the test points
    plot(points(:,1),points(:,2),'b.','LineWidth',5,'MarkerSize',40);

    % Make axis slightly larger
    temp = axis;
    axis_range_x = temp(2)-temp(1);
    axis_range_y = temp(4)-temp(3);
    percent_larger = 0.3;
    axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);

end

%% Fill in basic variables for calculations that follow
% How many points are on the path? How many points are we testing?
Npath = length(path(:,1));
Npoints = length(points(:,1));
Ndims = length(points(1,:));
% point_dimension = length(points(1,:));

% Find the differences, then square X, Y, (Z) terms, and take square root
% to find lengths
path_segment_lengths = sum(diff(path,1,1).^2,2).^0.5;
% Add them up to generate stations
path_stations = [0; cumsum(path_segment_lengths)];
% Make sure last point has a "length"
path_segment_lengths = [path_segment_lengths; path_segment_lengths(end,:)];
path_dimension = length(path(1,:));

%% STEP 1. Find the unit orthogonal vectors of each path segment
% Find unit orthogonal and unit tangent vectors for each segment in the
% path. Notes on this: The first path vector is from points 1 to 2, the second is
% points 2 to 3, etc. Tangential vectors are vectors along the directions
% of the segments Orthogonal vectors are ones that are 90 degrees, in
% positive cross-product direction, from the tangent vectors.
%
% These normal and tangential unit vectors are used to perform dot products
% to determine distances along a segment, or orthogonal to the segment.

% Call the function fcn_Path_findPathOrthogonalVectors to find ortho
% vectors. Set the flag_rounding_type to 1, which is the default - which is
% an ortho projection based on points at start and end of a segment
back_orthogonal_flag_rounding_type = 1;
[unit_orthogonal_vectors_at_midpoints, ~] = ...
    fcn_Path_findPathOrthogonalVectors(path,back_orthogonal_flag_rounding_type, -1);

if path_dimension == 2
    unit_tangential_vectors_at_midpoints = unit_orthogonal_vectors_at_midpoints*[0 -1; 1 0];
elseif  path_dimension == 3 
    unit_tangential_vectors_at_midpoints = unit_orthogonal_vectors_at_midpoints*[0 1 0; -1 0 0; 0 0 1];
else
    error('A 2D or 3D vector is expected');
end




% Show the vectors?
if flag_do_debug

    figure(fig_debug)

    % Plot the tangent vectors for ahead and behind
    for ith_point = 1:Npath-1
        midpoint = (path(ith_point,:) + path(ith_point+1,:))/2;
        quiver(...
            midpoint(1,1),...
            midpoint(1,2),...
            unit_tangential_vectors_at_midpoints(ith_point,1),...
            unit_tangential_vectors_at_midpoints(ith_point,2),...
            0,'g-','filled',...
            'LineWidth',3);
        quiver(...
            midpoint(1,1),...
            midpoint(1,2),...
            unit_orthogonal_vectors_at_midpoints(ith_point,1),...
            unit_orthogonal_vectors_at_midpoints(ith_point,2),...
            0,'c-','filled',...
            'LineWidth',3);
    end

end



% For each of the points, find cross products, dot products, and distances

% Initialize variables
closest_path_points         = zeros(Npoints,Ndims);
s_coordinates               = zeros(Npoints,1);
first_path_point_indicies   = zeros(Npoints,1);
second_path_point_indicies  = zeros(Npoints,1);
percent_along_length        = zeros(Npoints,1);
distances_real              = zeros(Npoints,1);
distances_imaginary         = zeros(Npoints,1);
ortho_vector_used           = zeros(Npoints,path_dimension);


% For each of the points, compare to the path to calculate key values
for ith_point = 1:Npoints
    %% STEP  2. Project from the end of each segment to the query point(s),
    %  This creates two vectors per segment, one at each end-point of each
    %  segment, pointing to each query point.
    projection_vector = ones(Npath,1)*points(ith_point,:) - path;
    segment_start_point_projection_to_test_point = projection_vector(1:end-1,:);
    segment_end_point_projection_to_test_point   = projection_vector(2:end,:);

    %% STEP 3. Take the cross and dot product of point-projection vectors
    % at each end of a line segment with the orthogonal projection vector
    % for the segment.
    crossProduct_ortho_vector_to_segment_start_vector = crossProduct(unit_orthogonal_vectors_at_midpoints,segment_start_point_projection_to_test_point);
    crossProduct_ortho_vector_to_segment_end_vector   = crossProduct(unit_orthogonal_vectors_at_midpoints,segment_end_point_projection_to_test_point);
    dotProduct_ortho_vector_to_segment_start_vector   = dotProduct(unit_orthogonal_vectors_at_midpoints,segment_start_point_projection_to_test_point);
    dotProduct_tangent_vector_to_segment_start_vector = dotProduct(unit_tangential_vectors_at_midpoints,segment_start_point_projection_to_test_point);

    %% STEP 4. Check if the cross-product changes
    %  If the cross-product signs change, then the point must be within the
    %  orthogonal projection of the line segment. Also, calculate the distance
    %  to each segment

    sign_changes = crossProduct_ortho_vector_to_segment_start_vector.*crossProduct_ortho_vector_to_segment_end_vector;
    negative_indicies = find(sign_changes<=0);
    if ~isempty(negative_indicies)
        [~,min_find_index] = min(abs(dotProduct_ortho_vector_to_segment_start_vector(negative_indicies)));
        min_segment_distance = dotProduct_ortho_vector_to_segment_start_vector(negative_indicies(min_find_index));
        segment_index = negative_indicies(min_find_index);
    else
        min_segment_distance = inf;
    end


    %% STEP 5. Check the vertex distances.
    %  There are cases where a point is at the vertex of a segment, and not
    %  within the projection of the line segment. To query these cases, find
    %  the distance of the query point to all points in the line segment,
    %  keeping the minimum.
    distances_squared = sum((path-points(ith_point,:)).^2,2);
    [min_distance_squared,min_vertex_distance_index] = min(distances_squared);
    min_vertex_distance = real(min_distance_squared.^0.5);

    %% Fill in results
    if abs(min_segment_distance) < abs(min_vertex_distance)
        % The segment is closer than any vertex
        closest_path_points(ith_point,:)         = path(segment_index,:) + dotProduct_tangent_vector_to_segment_start_vector(segment_index,:) * unit_tangential_vectors_at_midpoints(segment_index,:);
        s_coordinates(ith_point,:)               = path_stations(segment_index,:) + dotProduct_tangent_vector_to_segment_start_vector(segment_index,:);
        first_path_point_indicies(ith_point,:)   = segment_index;
        second_path_point_indicies(ith_point,:)  = segment_index + 1;
        percent_along_length(ith_point,:)        = dotProduct_tangent_vector_to_segment_start_vector(segment_index,:)/path_segment_lengths(segment_index,:);
        distances_real(ith_point,:)              = min_segment_distance;
        ortho_vector_used(ith_point,:)           = unit_orthogonal_vectors_at_midpoints(segment_index,:);
    else
        % A vertex is closer than any segment

        % This only occurs when the point is in the outside of a bend between front and back
        closest_path_points(ith_point,:)         = path(min_vertex_distance_index,:);
        s_coordinates(ith_point,:)               = path_stations(min_vertex_distance_index,:);
        first_path_point_indicies(ith_point,:)   = min_vertex_distance_index;
        second_path_point_indicies(ith_point,:)  = min_vertex_distance_index;
        percent_along_length(ith_point,:)        = 0;

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
        
        if min_vertex_distance_index == 1
            unit_orthogonal_projection_vector = unit_orthogonal_vectors_at_midpoints(1,:);
        elseif min_vertex_distance_index == Npath
            unit_orthogonal_projection_vector = unit_orthogonal_vectors_at_midpoints(end,:);
        else
            switch flag_rounding_type
                case 1 % Use the rear vector
                    unit_orthogonal_projection_vector = unit_orthogonal_vectors_at_midpoints(min_vertex_distance_index-1,:);

                case 2 % Use the front vector
                    unit_orthogonal_projection_vector = unit_orthogonal_vectors_at_midpoints(min_vertex_distance_index,:);

                case 3 % Use the average of the front and rear vectors
                    average_vector = ...
                        (unit_orthogonal_vectors_at_midpoints(min_vertex_distance_index-1,:) + ...
                        unit_orthogonal_vectors_at_midpoints(min_vertex_distance_index,:))/2;
                    mag_average_vector = sum(average_vector.^2,2).^0.5;
                    unit_orthogonal_projection_vector = average_vector./mag_average_vector;

                case 4
                    % This is a very special case, handled separately
                    error('This method is not coded yet! Try again later, maybe?');
                otherwise
                    error('Unknown flag_rounding_type');
            end
        end
        ortho_vector_used(ith_point,:)           = unit_orthogonal_projection_vector;
        % Calculate real and imaginary portions
        [distances_real(ith_point,:), distances_imaginary(ith_point,:)] = ...
            fcn_INTERNAL_convertOffsetCoordinateToImaginaryNumber(unit_orthogonal_projection_vector, closest_path_points(ith_point,:), points(ith_point,:));

        % If a point is before the start or after the end, need to update
        % coordinates
        if min_vertex_distance_index == 1
            percent_along_length(ith_point,:) = dotProduct_tangent_vector_to_segment_start_vector(1,:)/path_segment_lengths(1,:);
            s_coordinates(ith_point,:)        = path_segment_lengths(1,:)*percent_along_length(ith_point,:);
        elseif min_vertex_distance_index == Npath
            percent_along_length(ith_point,:) = dotProduct_tangent_vector_to_segment_start_vector(end,:)/path_segment_lengths(end,:);
            s_coordinates(ith_point,:)        = path_stations(end,:) + path_segment_lengths(end,:)*(percent_along_length(ith_point,:) - 1);
        end

    end

end


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

    % Is this a 2-D query?
    % if length(path(1,:))==2
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


    % Plot the tangent and ortho vectors
    for ith_point = 1:Npath-1
        midpoint = (path(ith_point,:) + path(ith_point+1,:))/2;
        quiver(...
            midpoint(1,1),...
            midpoint(1,2),...
            unit_tangential_vectors_at_midpoints(ith_point,1),...
            unit_tangential_vectors_at_midpoints(ith_point,2),...
            0,'g-','filled',...
            'LineWidth',3);
        quiver(...
            midpoint(1,1),...
            midpoint(1,2),...
            unit_orthogonal_vectors_at_midpoints(ith_point,1),...
            unit_orthogonal_vectors_at_midpoints(ith_point,2),...
            0,'c-','filled',...
            'LineWidth',3);
    end

    % Plot the reference vector
    for ith_point = 1:Npoints
        quiver(...
            closest_path_points(ith_point,1),...
            closest_path_points(ith_point,2),...
            ortho_vector_used(ith_point,1)*distances_real(ith_point), ...
            ortho_vector_used(ith_point,2)*distances_real(ith_point),...
            0,'Color',[0.5 0.5 0.5],'LineWidth',3);
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


    % elseif length(path(1,:))==3
    %
    %     % Plot the path
    %     plot3(path(:,1),path(:,2),path(:,3),'r-','Linewidth',5);
    %     plot3(path(:,1),path(:,2),path(:,3),'ro','Markersize',20);
    %
    %     axis equal;
    %
    %     % Plot the query point
    %     plot3(points(:,1),points(:,2),points(:,3),'ko');
    %     text( points(:,1),points(:,2),points(:,3),'Query point');
    %
    %     % Plot the closest path points;
    %     plot3(...
    %         path(first_path_point_indicies:second_path_point_indicies,1),...
    %         path(first_path_point_indicies:second_path_point_indicies,2),...
    %         path(first_path_point_indicies:second_path_point_indicies,3),'r*');
    %
    %
    %     % Plot the closest point on path
    %     plot3(...
    %         closest_path_vertex_points(:,1),...
    %         closest_path_vertex_points(:,2),...
    %         closest_path_vertex_points(:,3),...
    %         'go','Markersize',20);
    %     text(...
    %         closest_path_vertex_points(:,1),...
    %         closest_path_vertex_points(:,2),...
    %         closest_path_vertex_points(:,3),...
    %         'Snap Point on Path');
    %
    %
    %     % Connect closest point on path to query point
    %     plot3(...
    %         [points(:,1) closest_path_vertex_points(:,1)],...
    %         [points(:,2) closest_path_vertex_points(:,2)],...
    %         [points(:,3) closest_path_vertex_points(:,3)],...
    %         'g-','Linewidth',2);
    %
    %     view(3);
    %
    % end
end % Ends the flag_do_debug if statement


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
%% Calculate cross products
    function result = crossProduct(v,w)
        result = v(:,1).*w(:,2)-v(:,2).*w(:,1);
    end

%% Calculate dot products
    function result = dotProduct(v,w)
        result = sum(v.*w,2);
    end

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

