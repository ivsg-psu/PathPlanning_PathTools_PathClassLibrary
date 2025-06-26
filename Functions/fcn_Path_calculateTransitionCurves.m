function [TransitionCurves, ...
    closest_path_point1, closest_path_point2,...
    angle_start_radians, angle_end_radians] = ...
    fcn_Path_calculateTransitionCurves(path_1, path_2,radius_of_curve,varargin)
%% fcn_Path_calculateTransitionCurves
% Fits an arc to connect two paths with a user-defined radius. The
% direction of the paths are assumed to be such that the points are ordered
% starting at point 1 of path_1 and progressing forward with each point.
% The direction of travel on path_2 is also assumed to start at its point 1
% and progressing forward with each point. The resulting arc thus connects
% a path of positive travel from path_1 to positive travel on path_2.
%
% The function also returns: the center of the transition curve, the
% intersection points of the transition curve with the two given paths,
% and the angles that the intersection points make with the
% x-axis.
%
% FORMAT:
%
%       [ TransitionCurves, ...
%         closest_path_point1, closest_path_point2,...
%         angle_point1_deg,angle_point2_deg] = ...
%         fcn_Path_calculateTransitionCurves(...
%         path_1, ...
%         path_2,...
%         radius_of_curve, ...
%         (number_of_points_on_curve),...
%         (plot_color),...
%         (plot_line_width),...
%         (plot_text),...
%         (fig_num))
%
% INPUTS:
%
%       path_1 : data array in XY coordinates in [N x 2] that will be the
%       entry path into the transition curve
%
%       path_2 : data array in XY coordinates in [N x 2] that will be the
%       exit path out of the transition curve
%
%       radius_of_curve : radius of the transition curve between the two
%       line segmnets. Is positive for a curve to the left side according
%       to the direction of travel and negative to the right side of the
%       direction of travel.
%
%      (OPTIONAL INPUTS)
%
%      number_of_points_on_curve: the number of points the user wants on
%      the transition curve which will be equidistant throughout the curve
%
%      plot_color: the color of the plot of the transition curve as a string (eg. 'yellow' or 'y') or as a
%      RGB matrix (eg. [1 1 0] for yellow)
%
%      plot_line_width: the line width of the plot of the transition curve
%      in the forn of a number
%
%      plot_text: text to be shown on the plot in the form of a string (eg.
%      'Transition curve 1') in the color defined by plot_color
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%      TransitionCurves: array of XY coordinates of the transition curve
%      between the given line segments and with the given radius. The array
%      always starts on Intersection_point1 and ends at
%      Intersection_point2.
%
%      Intersection_point1 : Point on line segment 1 that intersects the
%      Intersection_point1 : Point on path_1 that intersects the
%      transition curve orthogonally, this is the first point on the
%      transition curve
%
%      Intersection_point2 : Point on path_2 that intersects the
%      transition curve orthogonally, this is the last point on the
%      transition curve
%
%      angle_start_radians: Angle of the intersection point between the
%      transition curve and path_1 measured from the center of
%      the transition curve, relative to the x-axis
%
%      angle_end_radians: Angle of the intersection point between the
%      transition curve and path_2 measured from the center of
%      the transition curve, relative to the x-axis
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_Path_calculateTransitionCurves.m for a full
%       test suite.
%
% This function was written on 2023_07_17 by V. Wagh
% Questions or comments? vbw5054@psu.ed


% Revision history:
% 2023_07_17 by V. Wagh
% -- start writing function
% 2023_07_21 by V. Wagh
% -- added inputs of direction for both line segments
% -- added outputs of intersection points and angles
% -- used atan2 to get the transition curve
% 2023_07_27 by X. Cao
% -- fixed the angle calculation bugs in atan2
% 2023_07_30 by V. Wagh
% -- added the plot_color, plot_line_width and plot_text as optional inputs
% -- removed segment1 and 2 directions as inputs
% 2023_07_31 to 2023_08_02 by S. Brennan
% -- minor reformatting of comments
% -- renamed "line_segment" to "path" throughout as the code was written
% for collections of line segments (paths), not for just line segments
% -- switched angle_point1_radians, angle_point2_radians instead of degrees
% as MATLAB default is to use radians for all angle calculations. Also, it
% appeared that the code was actually giving radians, not degrees (BUG!)
% -- fixed bug where the line segments at start/end were not converted into
% unit vectors before multiplying by a distance, which would have caused
% severe errors.
% -- fixed bug where the cross-product was giving the wrong angle if angles
% were larger than 90 degrees. Required use of both cross and dot product,
% and functionalized this
% -- functionalized creation of unit vectors
% -- fixed lots of the plotting to make it cleaner
% -- removed the offset calculation from path library, and instead moved it
% internal as sub-function as the Path library was using an offset option
% (from vertex) that isn't the same as a true offset
% -- found that, for concave path areas, the intersection calculations were
% wrong. Fixed this bug by calculating intersections one segment at a time.
% This required significant for loops and could be improved with
% optimization, but would require a rewrite of core path library functions
% (not hard, just not urgent)
% 2023_08_09 by S. Brennan
% -- Converted this function over to Path library from LoadWZ library
% 2025_06_23 - S. Brennan
% -- Updated debugging and input checks

% TO-DO
% (none)

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
        narginchk(3,MAX_NARGIN);

        % Check the path_1 input
        fcn_DebugTools_checkInputsToFunctions(path_1, 'path');

        % Check the path_2 input
        fcn_DebugTools_checkInputsToFunctions(path_2, 'path');

        % Check the radius_of_curve input
        fcn_DebugTools_checkInputsToFunctions(radius_of_curve, '1column_of_numbers',[1 1]);

    end


end

number_of_points_on_curve = 50; % Default is 50 points on the transition curve
if 4 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        number_of_points_on_curve = temp;

        % make sure n >= 2
        if number_of_points_on_curve < 2
            error('The number of points on transition curve must be greater than or equal to 2, please increase the number');
        end
    end
end

plot_color = 'red'; % Default color to plot the transition curve is red
if 5 <= nargin
    temp = varargin{2};
    if ~isempty(temp)
        plot_color = temp;
    end
end

plot_line_width = 2; % Default line width to plot the transition curve is 2
if 6 <= nargin
    temp = varargin{3};
    if ~isempty(temp)
        plot_line_width = temp;
    end
end

% Does user want to specify plot_text?
plot_text = 'Transition_Curve'; % Default
if 7 <= nargin
    temp = varargin{4};
    if ~isempty(temp)
        plot_text = temp;
    end
end


% Does user want to show the plots?
flag_do_plots = 0; % Default is to NOT show plots
fig_debug = [];
if (0==flag_max_speed) && (MAX_NARGIN == nargin)
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        fig_num = temp;
        figure(fig_num);
        flag_do_plots = 1;
    end
else
    if flag_do_debug
        fig_debug = 9999;
        flag_do_plots = 1;
    end
end


%% Main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 1:
% extend segment 1 end by amount of radius+10% and extend segment 2
% both sides by amount of radius+10%

% Calculate unit vectors
% FORMAT: unit_vector = fcn_INTERNAL_calcUnitVector(point_start,point_end)

% vector of line segment1 start
unit_v_lineSegment1_start = fcn_INTERNAL_calcUnitVector(path_1(1,:),path_1(2,:));

% unit vector of line segment1 end
unit_v_lineSegment1_end = fcn_INTERNAL_calcUnitVector(path_1(end-1,:),path_1(end,:));

% vector of line segment2 start
unit_v_lineSegment2_start = fcn_INTERNAL_calcUnitVector(path_2(1,:),path_2(2,:));

% vector of line segment2 end
unit_v_lineSegment2_end = fcn_INTERNAL_calcUnitVector(path_2(end-1,:),path_2(end,:));


% distance to extend the lines is 110% of the radius
extension_distance = radius_of_curve * 1.1;

% get the x and y coordinates of the points after line extension
path_1_extended_points_end = path_1(end,:) + (unit_v_lineSegment1_end   * extension_distance); % extending from end
path_2_extended_points_start = path_2(1,:) - (unit_v_lineSegment2_start * extension_distance); % extending from start BACKWARDS
path_2_extended_points_end = path_2(end,:) + (unit_v_lineSegment2_end   * extension_distance); % extending from end

% add the extended points to the array for line segments
path_1_extended = [path_1; path_1_extended_points_end];
path_2_extended = [path_2_extended_points_start; path_2; path_2_extended_points_end];

% Show the results?
if flag_do_debug == 1
    figure(fig_debug);
    clf;
    hold on;
    grid on;
    xlabel('X [m]');
    ylabel('Y [m]');
    axis equal;


    % Show the inputs
    plot(path_1(:,1), path_1(:,2),'-','color',[0 1 0],'LineWidth',7);
    plot(path_2(:,1), path_2(:,2),'-','color',[0 0 1],'LineWidth',5);

    % Show the extensions
    plot(path_1_extended(:,1), path_1_extended(:,2),'-','color',[0 1 0]*0.8,'LineWidth',3);
    plot(path_2_extended(:,1), path_2_extended(:,2),'-','color',[0 0 1]*0.8,'LineWidth',1.5);

    % Show the directions
    quiver(path_1(1,1), path_1(1,2),unit_v_lineSegment1_start(1,1),unit_v_lineSegment1_start(1,2),0,'-','LineWidth',5,'ShowArrowHead','on','Color',[0 1 0],'MaxHeadSize',4)
    quiver(path_2(1,1), path_2(1,2),unit_v_lineSegment2_start(1,1),unit_v_lineSegment2_start(1,2),0,'-','LineWidth',3,'ShowArrowHead','on','Color',[0 0 1],'MaxHeadSize',4)
    title('Original Line Segments and their Extentions','fontsize',24)
end

%% convert the extended paths to traversals
% This is so that we can use them in the PathClass Library function to find
% intersectinos and offsets to the path

% Convert the extended versions
path_1_extended_traversal = fcn_Path_convertPathToTraversalStructure(path_1_extended,-1);
path_2_extended_traversal = fcn_Path_convertPathToTraversalStructure(path_2_extended,-1);


%% 1. calculate if and how the extended paths intersect
% 1a. use function from PathClass to find the intersection
% 1b. check the type of intersection
% 1c. find the segments that are intersecting

% Step 1a. find the intersections
[extended_paths_intersections,...
    s_coordinates_in_traversal_1,...
    s_coordinates_in_traversal_2] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    path_1_extended_traversal,...
    path_2_extended_traversal, -1);

% Show the results?
if flag_do_debug == 1
    figure(fig_debug);
    plot(extended_paths_intersections(:,1), extended_paths_intersections(:,2),'color','black','Marker','o','MarkerSize',10);
end

%% 1b. check the type of intersection

% get the intersection point of the extended paths
number_of_extended_intersections = length(extended_paths_intersections(:,1));

if number_of_extended_intersections == 0
    % If paths never intersect, throw an error!
    error('Paths never intersect, transition curve not possible');
elseif isequal(path_1, path_2)
    % If paths are identical, throw an error
    error('Paths are collinear, curve is not possible between them');

elseif number_of_extended_intersections >= 1
    % THIS IS THE TYPICAL CASE:

    % Check for infinite intersections - this will show up if the first
    % point of intersection, then second point of intersection, etc. are
    % all on path_1, and the sequence is either on the same segment of
    % path_1 or follows the end-points of path_1 (or both). If there are
    % infinite, we want the LAST point of correspondence as the "first"
    % intersection of point.



    % 1.b.1 if they intersect multiple times, use first as intersection pt ( flag as
    % intersection exists)
    extended_intersection_point = extended_paths_intersections(1,:);
    s_coordinate_in_traversal_1 = s_coordinates_in_traversal_1(1);
    s_coordinate_in_traversal_2 = s_coordinates_in_traversal_2(1);
else
    error('Unknown case detected in an if statement - debugging required.')
end

%% 1c. find the segments that are intersecting
% Then use these to find the vectors in/out
prior_segment_path_1 = fcn_INTERNAL_findPriorIntersectingSegment(path_1_extended_traversal,s_coordinate_in_traversal_1);
post_segment_path_2  = fcn_INTERNAL_findPostIntersectingSegment(path_2_extended_traversal,s_coordinate_in_traversal_2);

% Show the results?
if flag_do_debug == 1
    figure(fig_debug);
    plot(...
        [path_1_extended(prior_segment_path_1,1) extended_intersection_point(1,1)],...
        [path_1_extended(prior_segment_path_1,2) extended_intersection_point(1,2)],...
        '-','color',[1 0 1],'LineWidth',2);
    plot(...
        [extended_intersection_point(1,1) path_2_extended(post_segment_path_2+1,1)],...
        [extended_intersection_point(1,2) path_2_extended(post_segment_path_2+1,2)],...
        '-','color',[0 1 1],'LineWidth',2);
end

unit_vector_into_intersection  = fcn_INTERNAL_calcUnitVector(path_1_extended(prior_segment_path_1,:),extended_intersection_point);
unit_vector_outof_intersection = fcn_INTERNAL_calcUnitVector(extended_intersection_point,path_2_extended(post_segment_path_2+1,:));

% Show the results?
if flag_do_debug == 1
    figure(fig_debug);
    start_point_1 = extended_intersection_point - unit_vector_into_intersection;
    quiver(start_point_1(1,1), start_point_1(1,2),unit_vector_into_intersection(1,1),unit_vector_into_intersection(1,2),0,'-','LineWidth',5,'ShowArrowHead','on','Color',[1 0 1],'MaxHeadSize',4)
    quiver(extended_intersection_point(1,1), extended_intersection_point(1,2),unit_vector_outof_intersection(1,1),unit_vector_outof_intersection(1,2),0,'-','LineWidth',3,'ShowArrowHead','on','Color',[0 1 1],'MaxHeadSize',4)
end

%% offset the paths by the distance of the radius in the direction where radius will be positive as default

% 2 find what quadrant are allowed for the curve, find where segment 2 is
% wrt 1 and the intersection point
% 2.1 then take first segment of path_2 that sticks out in the direction
% of the radius and take the point of contact with the segment that sticks
% out and the path_1 as the intersection point
% 2.2 if no parts stick out in correct direction, throw an error
% (error! infinite overlap and unable to resolve a correct side to use)

% 3 find the angle between the intersecting segments to determine which
% side of traversal 2 to use

% Calculate the angle between the unit vectors
% FORMAT: angle_between_vectors_radians = fcn_INTERNAL_findAngleBetweenUnitVectors(from_unit_vector,to_unit_vector)
angle_change_radians = fcn_INTERNAL_findAngleBetweenUnitVectors(unit_vector_into_intersection,unit_vector_outof_intersection);

% get direction to offset the path segments
offset_direction_relative_to_1 = sign(angle_change_radians);
offset_amount = offset_direction_relative_to_1*radius_of_curve;

% offset the path segments in the in the direction of their travel (where
% direction of travel is in the order of the points in the matrix of the
% path segments)
[offset_traversal_1_start, offset_traversal_1_end] = fcn_INTERNAL_findOrthogonalOffsetTraversal(path_1_extended_traversal, offset_amount); % offset for traversal 1
[offset_traversal_2_start, offset_traversal_2_end] = fcn_INTERNAL_findOrthogonalOffsetTraversal(path_2_extended_traversal, offset_amount); % offset for traversal 2

% BELOW uses Path Library function, which is not actually a perfect offset
% - See Example 6 - because it is averaging the end points. These examples
% will fail because of this imperfection!
%
% offset_traversal_1 = fcn_Path_fillOffsetTraversalsAboutTraversal(path_1_extended_traversal, offset_amount,-1); % offset for traversal 1
% offset_traversal_2 = fcn_Path_fillOffsetTraversalsAboutTraversal(path_2_extended_traversal, offset_amount,-1); % offset for traversal 2

% Show the results?
if flag_do_debug == 1
    figure(fig_debug);
    % Plot the offsets, blending them by a bit of white to make them
    % lighter
    for ith_segment = 1:length(offset_traversal_1_start(:,1))
        plot(...
            [offset_traversal_1_start(ith_segment,1) offset_traversal_1_end(ith_segment,1)],...
            [offset_traversal_1_start(ith_segment,2) offset_traversal_1_end(ith_segment,2)],...
            '.-','color',[0 1 0]+0.7*[1 0 1],'LineWidth',3,'MarkerSize',30);

        text(...
            (offset_traversal_1_start(ith_segment,1)+offset_traversal_1_end(ith_segment,1))/2,...
            (offset_traversal_1_start(ith_segment,2)+offset_traversal_1_end(ith_segment,2))/2,...
            sprintf('offset1 %.0d',ith_segment),'color','k','Interpreter','none');
    end

    for ith_segment = 1:length(offset_traversal_2_start(:,1))

        plot(...
            [offset_traversal_2_start(ith_segment,1) offset_traversal_2_end(ith_segment,1)],...
            [offset_traversal_2_start(ith_segment,2) offset_traversal_2_end(ith_segment,2)],...
            '.-','color',[0 0 1]+0.7*[1 1 0],'LineWidth',3,'MarkerSize',30);

        text(...
            (offset_traversal_2_start(ith_segment,1)+offset_traversal_2_end(ith_segment,1))/2,...
            (offset_traversal_2_start(ith_segment,2)+offset_traversal_2_end(ith_segment,2))/2,...
            sprintf('offset2 %.0d',ith_segment),'color','k','Interpreter','none');
    end
    text(...
        offset_traversal_1_start(1,1),...
        offset_traversal_1_start(1,2),...
        'offset1','color','k','Interpreter','none');

    text(...
        offset_traversal_2_start(1,1),...
        offset_traversal_2_start(1,2),...
        'offset2','color','k','Interpreter','none');
end

%% find the intersection point of the offsets, this is the center of the transition curves
% The following is a horrible way to do this, as it is very computational.
% However, it works and is "clean" to understand. It would greatly help to
% rewrite the below search process to be faster.
% This works by taking every segment in path 1, segment by segment, and
% converting each to a traversal, then checking for an intersection with every
% segment in path 2 and keeping the closest intersection if one is found.

flag_an_intersection_was_found = 0;
center_transition_curve = nan(1,2);
closest_distance = inf;
path_1_segment_hit = nan;
path_2_segment_hit = nan;

% Convert all segment_2's into traversals
% Initialize with an empty struct
segment_1_traversals{length(offset_traversal_1_start(:,1))} = struct;
segment_2_traversals{length(offset_traversal_2_start(:,1))} = struct;

% Fill in all the traversals possible for path_1
for ith_segment_in_path1 = 1:length(offset_traversal_1_start(:,1))
    segment_1_traversals{ith_segment_in_path1} = ...
        fcn_Path_convertPathToTraversalStructure([offset_traversal_1_start(ith_segment_in_path1,:); offset_traversal_1_end(ith_segment_in_path1,:)],-1);
end

% Fill in all the traversals possible for path_2
for jth_segment_in_path2 = 1:length(offset_traversal_2_start(:,1))
    segment_2_traversals{jth_segment_in_path2} = ...
        fcn_Path_convertPathToTraversalStructure([offset_traversal_2_start(jth_segment_in_path2,:); offset_traversal_2_end(jth_segment_in_path2,:)],-1);
end

% Now loop through all the path 1 segments, checking each segment from path
% 2 for an intersection, keeping the closest intersection to start of path
% 1.
for ith_segment_in_path1 = 1:length(offset_traversal_1_start(:,1))
    if flag_an_intersection_was_found==0
        segment_1_traversal_to_check = segment_1_traversals{ith_segment_in_path1};

        for jth_segment_in_path2 = 1:length(offset_traversal_2_start(:,1))
            segment_2_traversal_to_check = segment_2_traversals{jth_segment_in_path2};

            [intersection_point,...
                s_coordinates_in_traversal_1,...
                ~] = ...
                fcn_Path_findIntersectionsBetweenTraversals(...
                segment_1_traversal_to_check,...
                segment_2_traversal_to_check, -1);

            if ~isempty(intersection_point)
                flag_an_intersection_was_found = 1;
                if s_coordinates_in_traversal_1<closest_distance
                    center_transition_curve = intersection_point;
                    closest_distance = s_coordinates_in_traversal_1;
                    path_1_segment_hit = ith_segment_in_path1;
                    path_2_segment_hit = jth_segment_in_path2;
                end
            end

        end
    end
end

if isempty(center_transition_curve)
    error('No intersections found - not possible to create a transition curve.');
end

% Show the results?
if flag_do_debug == 1
    figure(fig_debug);
    % plot the intersection point of the offset paths, this is the center of
    % the transition curve
    plot(center_transition_curve(:,1),center_transition_curve(:,2),'Color','k','Marker','*','MarkerSize',20);
    text(center_transition_curve(:,1),center_transition_curve(:,2),'Center used for arc','color','k','Interpreter','none');
end



%% find the closest point from the center of the transition curve to both paths
% these are the intersection points of the curve and the paths
% when the curve is tangent to the paths respectively

% Convert all segments in path 1 and 2 into traversals
% Initialize with an empty struct
path_1_traversals{length(offset_traversal_1_start(:,1))} = struct;
path_2_traversals{length(offset_traversal_2_start(:,1))} = struct;

% Fill in all the traversals possible for path_1
N_segments_path_1 = length(path_1_extended_traversal.X(:,1))-1;
N_segments_path_2 = length(path_2_extended_traversal.X(:,1))-1;

% Fill in all the traversals possible for path_1
for ith_segment_in_path1 = 1:N_segments_path_1
    path_1_traversals{ith_segment_in_path1} = fcn_Path_convertPathToTraversalStructure(...
        [path_1_extended_traversal.X(ith_segment_in_path1,1) path_1_extended_traversal.Y(ith_segment_in_path1,1); ...
        path_1_extended_traversal.X(ith_segment_in_path1+1,1) path_1_extended_traversal.Y(ith_segment_in_path1+1,1)],-1);
end

% Fill in all the traversals possible for path_2
for jth_segment_in_path2 = 1:N_segments_path_2
    path_2_traversals{jth_segment_in_path2} = fcn_Path_convertPathToTraversalStructure(...
        [path_2_extended_traversal.X(jth_segment_in_path2,1) path_2_extended_traversal.Y(jth_segment_in_path2,1); ...
        path_2_extended_traversal.X(jth_segment_in_path2+1,1) path_2_extended_traversal.Y(jth_segment_in_path2+1,1)],-1);
end


% snap point onto path_1
[closest_path_point1,~,~,~,...
    ~,~] = ...
    fcn_Path_snapPointOntoNearestTraversal(center_transition_curve, path_1_traversals{path_1_segment_hit},-1);

% snap point onto path_2
[closest_path_point2,~,~,~,...
    ~,~] = ...
    fcn_Path_snapPointOntoNearestTraversal(center_transition_curve, path_2_traversals{path_2_segment_hit},-1);

% Check results
distance_1 = sum((center_transition_curve - closest_path_point1).^2,2).^0.5;
distance_2 = sum((center_transition_curve - closest_path_point2).^2,2).^0.5;
if abs(distance_1 - distance_2)>0.00001*distance_2
    error('Distances do not match! Stop here...')
end

% Show the results?
if flag_do_debug == 1
    figure(fig_debug);

    % plot the closest point fron the center to the original paths,
    % this is the intersection point of the transition curve
    % and paths, where the curve will be orthogonal to the paths
    plot(closest_path_point1(:,1),closest_path_point1(:,2),'Color','g','Marker','.','MarkerSize',40);
    text(closest_path_point1(:,1),closest_path_point1(:,2),'Transition_start','color','k','Interpreter','none','VerticalAlignment','bottom');
    plot(closest_path_point2(:,1),closest_path_point2(:,2),'Color','b','Marker','.','MarkerSize',40);
    text(closest_path_point2(:,1),closest_path_point2(:,2),'Transition_end','color','k','Interpreter','none','VerticalAlignment','bottom');
end

%% have a transition curve from start to end point with given radius and center of the transition curve


unit_vector_to_point1  = fcn_INTERNAL_calcUnitVector(center_transition_curve,closest_path_point1);
unit_vector_to_point2  = fcn_INTERNAL_calcUnitVector(center_transition_curve,closest_path_point2);

% Find angle between unit vectors
% FORMAT: angle_between_vectors_radians = fcn_INTERNAL_findAngleBetweenUnitVectors(from_unit_vector,to_unit_vector)
angle_change_radians = fcn_INTERNAL_findAngleBetweenUnitVectors(unit_vector_to_point1,unit_vector_to_point2);

% calculation of the angle for the intersection points wrt the x-axis

% make sure n >= 2
if number_of_points_on_curve < 2
    error('This number of points on transition curve is not possible, please increase the number');
else
    % continue with code
end

start_center_diff = closest_path_point1 - center_transition_curve;
angle_start_radians = atan2(start_center_diff(:,2), start_center_diff(:,1));
angle_end_radians = angle_start_radians + angle_change_radians;

% to get n number of equally spaced points on the transition curve
thetas_radians = (linspace(angle_start_radians,angle_end_radians,number_of_points_on_curve))';
TransitionCurves = center_transition_curve + radius_of_curve*[cos(thetas_radians) sin(thetas_radians)]; % data for n points on circle

% Show the results?
if flag_do_debug == 1
    figure(fig_debug);

    % Show the transition curve
    plot(TransitionCurves(:,1),TransitionCurves(:,2),'r.-','LineWidth',2);

    % plot the start and end points of the curve after using atan2
    plot(TransitionCurves(1,1), TransitionCurves(1,2),'color','green','Marker','o','MarkerSize',20);
    plot(TransitionCurves(end,1), TransitionCurves(end,2),'color','blue','Marker','o','MarkerSize',20);
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

if flag_do_plots == 1 % only plot if the user has given a fig_num

    %preparing the figure
    figure(fig_num);
    clf;
    hold on;
    grid on;
    xlabel('X [m]');
    ylabel('Y [m]');
    axis equal;

    % Show the inputs
    plot(path_1(:,1), path_1(:,2),'-','color',[0 1 0],'LineWidth',1.5);
    plot(path_2(:,1), path_2(:,2),'-','color',[0 0 1],'LineWidth',1.5);

    % Show the directions
    quiver(path_1(1,1), path_1(1,2),unit_v_lineSegment1_start(1,1),unit_v_lineSegment1_start(1,2),0,'-','LineWidth',3,'ShowArrowHead','on','Color',[0 1 0],'MaxHeadSize',4)
    quiver(path_2(1,1), path_2(1,2),unit_v_lineSegment2_start(1,1),unit_v_lineSegment2_start(1,2),0,'-','LineWidth',3,'ShowArrowHead','on','Color',[0 0 1],'MaxHeadSize',4)

    % Show the center of the arc
    plot(center_transition_curve(:,1),center_transition_curve(:,2),'Color','k','Marker','*','MarkerSize',20);
    text(center_transition_curve(:,1),center_transition_curve(:,2),'Center used for arc','color','k','Interpreter','none');

    % Show the transition curve
    plot(TransitionCurves(:,1),TransitionCurves(:,2),'.-','Color',plot_color, 'LineWidth',plot_line_width);

    % plot the start and end points of the curve after using atan2
    plot(TransitionCurves(1,1), TransitionCurves(1,2),'color','green','Marker','o','MarkerSize',20);
    plot(TransitionCurves(end,1), TransitionCurves(end,2),'color','blue','Marker','o','MarkerSize',20);
    text(TransitionCurves(1,1), TransitionCurves(1,2),plot_text,'Interpreter','none');
end


if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end
end % Ends main function

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

%% fcn_INTERNAL_calcUnitVector
function unit_vector = fcn_INTERNAL_calcUnitVector(point_start,point_end)
vector_to_calculate    = point_end - point_start;
magnitude_vector_to_calculate = sum(vector_to_calculate.^2,2).^0.5;
unit_vector = vector_to_calculate./magnitude_vector_to_calculate;
end % Ends fcn_INTERNAL_calcUnitVector

%% fcn_INTERNAL_findPriorIntersectingSegment
function prior_intersecting_segment_number = fcn_INTERNAL_findPriorIntersectingSegment(traversal_to_check,s_coordinate)
path_1_segment_intersecting_start = find(traversal_to_check.Station<s_coordinate,1,'last');
path_1_segment_intersecting_end   = find(traversal_to_check.Station>s_coordinate,1,'first');
path_1_segment_intersecing_on     = interp1(traversal_to_check.Station,(1:length(traversal_to_check.Station))',s_coordinate);
prior_intersecting_segment_number       = floor(path_1_segment_intersecing_on);

% If the intersection is on a vertex, then the current intersecting segment
% may be greater than the prior one
if path_1_segment_intersecting_start < prior_intersecting_segment_number
    prior_intersecting_segment_number = path_1_segment_intersecting_start;
end

% Check if the intersection is exactly at the end of the traversal - if so,
% then need to use the prior segment number
if isempty(path_1_segment_intersecting_end) || (prior_intersecting_segment_number == length(traversal_to_check.Station))
    prior_intersecting_segment_number = length(traversal_to_check.Station) - 1;
end
end % Ends fcn_INTERNAL_findPriorIntersectingSegment

%% fcn_INTERNAL_findPostIntersectingSegment
function post_intersecting_segment_number = fcn_INTERNAL_findPostIntersectingSegment(traversal_to_check,s_coordinate)
path_1_segment_intersecting_end   = find(traversal_to_check.Station>s_coordinate,1,'first');
path_1_segment_intersecing_on     = interp1(traversal_to_check.Station,(1:length(traversal_to_check.Station))',s_coordinate);
post_intersecting_segment_number       = floor(path_1_segment_intersecing_on);

% If the intersection is on a vertex, then the rounding down will match the
% current position. If so, we want to use the vertex that follows
if path_1_segment_intersecing_on == post_intersecting_segment_number
    post_intersecting_segment_number = post_intersecting_segment_number+1;
end

% Check if the intersection is exactly at the end of the traversal - if so,
% then need to use the prior segment number
if isempty(path_1_segment_intersecting_end) || (post_intersecting_segment_number >= length(traversal_to_check.Station))
    post_intersecting_segment_number = length(traversal_to_check.Station) - 1;
end
end % Ends fcn_INTERNAL_findPostIntersectingSegment

%% fcn_INTERNAL_findOrthogonalOffsetTraversal
function [offset_segments_start, offset_segments_end] = fcn_INTERNAL_findOrthogonalOffsetTraversal(reference_traversal,offset)

N_segments = length(reference_traversal.Station) - 1;
offset_segments_start = zeros(N_segments,2);
offset_segments_end   = zeros(N_segments,2);



% For each segment, project an offset, and find its coordinates
for ith_segment = 1:N_segments
    segment_coordinates = [reference_traversal.X(ith_segment:ith_segment+1) reference_traversal.Y(ith_segment:ith_segment+1)];
    segment_traversal = fcn_Path_convertPathToTraversalStructure(segment_coordinates,-1);
    segment_stations = segment_traversal.Station;

    % Set the projection type to use. 1 indicates using the prior segment
    % to project a vertex
    projection_type = 1;
    [unit_normal_vector_start, unit_normal_vector_end] = ...
        fcn_Path_findOrthogonalTraversalVectorsAtStations(...
        segment_stations(:,1),segment_traversal,projection_type,-1);

    % Use the unit vectors to find the offsets
    unit_vectors = unit_normal_vector_end - unit_normal_vector_start;
    offset_segment = unit_normal_vector_start + unit_vectors.*offset;


    % Save results into an array of offset segments
    offset_segments_start(ith_segment,:) = offset_segment(1,:);
    offset_segments_end(ith_segment,:)   = offset_segment(2,:);
end

% For debugging
if 1==0
    figure(12323);
    clf
    hold on;
    grid on;
    axis equal;
    plot(reference_traversal.X,reference_traversal.Y,'b.-','MarkerSize',10);
    for ith_segment = 1:N_segments
        plot(...
            [offset_segments_start(ith_segment,1) offset_segments_end(ith_segment,1)],...
            [offset_segments_start(ith_segment,2) offset_segments_end(ith_segment,2)],...
            'r.-','MarkerSize',10);
    end
end

end % Ends fcn_INTERNAL_findOrthogonalOffsetTraversal


%% fcn_INTERNAL_findAngleBetweenUnitVectors
function angle_between_vectors_radians = fcn_INTERNAL_findAngleBetweenUnitVectors(from_unit_vector,to_unit_vector)
cross_result = cross([from_unit_vector(1,1:2) 0],[to_unit_vector(1,1:2) 0]);
dot_result = sum(from_unit_vector(1,1:2).*to_unit_vector(1,1:2),2);
angle_magnitude = real(acos(dot_result));
angle_sign = sign(cross_result(3));
angle_between_vectors_radians = angle_sign*angle_magnitude;
end % Ends fcn_INTERNAL_findAngleBetweenUnitVectors