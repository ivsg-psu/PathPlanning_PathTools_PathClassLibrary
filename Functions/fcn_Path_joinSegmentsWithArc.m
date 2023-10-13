function [arc_points, line_segment_1, line_segment_2,...
    radius_of_arc, circle_center,...
    start_angle_radians, end_angle_radians] = ...
    fcn_Path_joinSegmentsWithArc(segment_1, segment_2,varargin)
%% fcn_Path_joinSegmentsWithArc
% Fits an arc to connect two line segments.
%
% The function also returns: the center of the transition curve, the
% intersection points of the transition curve with the two given paths,
% and the angles that the intersection points make with the
% x-axis.
%
% FORMAT:
%
%       [arc_points, line_segment_1, line_segment_2,...
%        radius_of_arc, circle_center,...
%        start_angle_radians, end_angle_radians] = ...
%         fcn_Path_joinSegmentsWithArc(...
%         segment_1, ...
%         segment_2,...
%         (number_of_points_on_curve),...
%         (fig_num))
%
% INPUTS:
%
%       segment_1 : data array in XY coordinates in [2 x 2] that will be the
%       entry path into the transition curve
%
%       segment_1 : data array in XY coordinates in [2 x 2] that will be the
%       exit path out of the transition curve
%
%      (OPTIONAL INPUTS)
%
%      number_of_points_on_curve: the number of points the user wants on
%      the transition curve which will be equidistant throughout the curve.
%      Default is 5;
%      
%      fig_num: a figure number to plot result
%
% OUTPUTS:
%
%      arc_points: array of XY coordinates of the arc connecting segment_1
%      to segment_2, ordered from segment_1 to segment_2
%
%      line_segment_1: array of XY points connecting end of segment_1 to
%      the start of the arc. It will be [nan nan] if arc starts exactly at
%      end of segment_1.
%
%      line_segment_2: array of XY points connecting end of arc to the
%      start of the segment_2. It will be [nan nan] if arc ends exactly at
%      start of segment_2.
%
%      radius_of_arc : radius of the curve between the two
%       line segments.
%   
%      circle_center: the coordinates for the center of the circle defining
%      the arc.
%
%      start_angle_radians : the angle, in radians, where the arc starts
%      relative to the center of the circle. Angles are measured relative
%      to the x-axis.
%
%      end_angle_radians : the angle, in radians, where the arc ends
%      relative to the center of the circle. Positive values are
%      counter-clockwise, negative values are counter-clockwise.
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_Path_joinSegmentsWithArc.m for a full
%       test suite.
%
% This function was written on 2023_10_09 by S. Brennan
% Questions or comments? sbrennan@psu.edu


% Revision history:
% 2023_10_09 by S. Brennan
% -- start writing function


flag_do_debug = 0; % Flag to show the results for debugging
% flag_do_plots = 0; % Flag to plot the final results (set below, so
% commented out)
flag_check_inputs = 1; % Flag to perform input checking

% Tell user where we are
if flag_do_debug
    st = dbstack; 
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
end
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

end

% Does user want to specify number_of_points_on_curve?
number_of_points_on_curve = 5; % Default is 5 points on the transition curve
if 3 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        number_of_points_on_curve = temp;

        % make sure n >= 2
        if number_of_points_on_curve < 2
            error('The number of points on transition curve must be greater than or equal to 2, please increase the number');
        end
    end
end

% Does user want to specify fig_num?
flag_do_plots = 0; % Default is not to plot the data
if 4 <= nargin
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end

% Setup figures if there is debugging
if flag_do_debug
    fig_debug = 9999; 
else
    fig_debug = []; %#ok<UNRCH> 
end

%% Write main code for plotting
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
% find unit vectors for both segments, determine if connection is to the
% left (positive) or right (negative) in going from segment 1 to 2

% Calculate unit vectors
% FORMAT: unit_vector = fcn_INTERNAL_calcUnitVector(point_start,point_end)

% vector of line segment1 start
unit_v_lineSegment1 = fcn_INTERNAL_calcUnitVector(segment_1(1,:),segment_1(2,:));
unit_v_lineSegment1_orthogonal = unit_v_lineSegment1*[0 1; -1 0];

% vector of line segment2 start
unit_v_lineSegment2 = fcn_INTERNAL_calcUnitVector(segment_2(1,:),segment_2(2,:));
unit_v_lineSegment2_orthogonal = unit_v_lineSegment2*[0 1; -1 0];

cross_product_unit_vectors = crossProduct(unit_v_lineSegment1,unit_v_lineSegment2);
flag_point_to_left = cross_product_unit_vectors>0;

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
    plot(segment_1(:,1), segment_1(:,2),'.-','color',[0 1 0],'LineWidth',7,'MarkerSize',30);
    plot(segment_2(:,1), segment_2(:,2),'.-','color',[0 0 1],'LineWidth',5,'MarkerSize',30);


    % Show the directions
    quiver(segment_1(2,1), segment_1(2,2),unit_v_lineSegment1(1,1),unit_v_lineSegment1(1,2),0,'-','LineWidth',5,'ShowArrowHead','on','Color',[0 1 0],'MaxHeadSize',4)
    quiver(segment_2(1,1), segment_2(1,2),unit_v_lineSegment2(1,1),unit_v_lineSegment2(1,2),0,'-','LineWidth',3,'ShowArrowHead','on','Color',[0 0 1],'MaxHeadSize',4)

    quiver(segment_1(2,1), segment_1(2,2),unit_v_lineSegment1_orthogonal(1,1),unit_v_lineSegment1_orthogonal(1,2),0,'-','LineWidth',5,'ShowArrowHead','on','Color',[0 1 0],'MaxHeadSize',4)
    quiver(segment_2(1,1), segment_2(1,2),unit_v_lineSegment2_orthogonal(1,1),unit_v_lineSegment2_orthogonal(1,2),0,'-','LineWidth',3,'ShowArrowHead','on','Color',[0 0 1],'MaxHeadSize',4)

    plot(segment_1(:,1), segment_1(:,2),'.','color',[0 0 0],'LineWidth',7,'MarkerSize',15);
    plot(segment_2(:,1), segment_2(:,2),'.','color',[0 0 0],'LineWidth',5,'MarkerSize',15);
    title('Original Line Segments and their unit vectors','fontsize',12)
end

%% Find the intersection
% FORMAT: 
%      [distance,location,path_segment, t, u] = ...
%         fcn_Path_findProjectionHitOntoPath(path,...
%         sensor_vector_start,sensor_vector_end,...
%         (flag_search_type),(fig_num))  

[~,~,~, percentage_distance_along_segment_2, percentage_distance_along_segment_1] = ...
    fcn_Path_findProjectionHitOntoPath(segment_2,segment_1(1,:),segment_1(2,:),1);

length_segment_1 = sum((segment_1(2,:)-segment_1(1,:)).^2,2).^0.5;
length_segment_2 = sum((segment_2(2,:)-segment_2(1,:)).^2,2).^0.5;

hit_point_1 = segment_1(1,:) + length_segment_1*percentage_distance_along_segment_1*unit_v_lineSegment1;
hit_point_2 = segment_2(1,:) + length_segment_2*percentage_distance_along_segment_2*unit_v_lineSegment2;

% Check that we get the same answer from both segments for where they cross
if max(abs(hit_point_2 - hit_point_1))>1E-8
    error('Different intersection points calculated from each of the segments. This indicates an error.');
end

% Check that distance is not equal to zero
if percentage_distance_along_segment_1>0 && percentage_distance_along_segment_1<1
    warning('on','backtrace');
    warning('Intersection was found between segments indicating that the segments cannot be joined with an arc, resulting in zero radius.');
    arc_points = [nan nan];
    line_segment_1 = [nan nan];
    line_segment_2 = [nan nan];
    radius_of_arc = 0;
    circle_center = hit_point_1;
    start_angle_radians = nan;
    end_angle_radians = nan;
    return;
end

% Find which one is closer
distance_squared_to_1 = sum((hit_point_1 - segment_1(2,:)).^2,2);
distance_squared_to_2 = sum((hit_point_1 - segment_2(1,:)).^2,2);

if distance_squared_to_2 < distance_squared_to_1
    flag_point_1_is_closest = 0;
    distance = abs(distance_squared_to_2)^0.5;
else
    flag_point_1_is_closest = 1;    
    distance = abs(distance_squared_to_1)^0.5;
end

% Show the results?
if flag_do_debug == 1
    figure(fig_debug);

    % Show the hit_point_1
    plot(hit_point_1(:,1), hit_point_1(:,2),'+','color',[1 1 0],'LineWidth',7,'MarkerSize',30);
end

%% Find the apex angle
angle_between_vectors_radians  = fcn_INTERNAL_findAngleBetweenUnitVectors(unit_v_lineSegment1, unit_v_lineSegment2);
apex_angle_radians = pi - abs(angle_between_vectors_radians);
apex_angle_degrees = apex_angle_radians*180/pi;

%% Find the side distance and center of the circle
radius_of_arc = distance*tan(apex_angle_radians/2);

if flag_point_1_is_closest
    end_point_on_1 = segment_1(2,:);
    if flag_point_to_left
        circle_center = segment_1(2,:) + unit_v_lineSegment1_orthogonal*radius_of_arc;
        start_point_on_2 = circle_center - unit_v_lineSegment2_orthogonal*radius_of_arc;
    else
        circle_center = segment_1(2,:) - unit_v_lineSegment1_orthogonal*radius_of_arc;
        start_point_on_2 = circle_center + unit_v_lineSegment2_orthogonal*radius_of_arc;
    end
else
    start_point_on_2 = segment_2(1,:);
    if flag_point_to_left
        circle_center = segment_2(1,:) + unit_v_lineSegment2_orthogonal*radius_of_arc;
        end_point_on_1 = circle_center - unit_v_lineSegment1_orthogonal*radius_of_arc;
    else
        circle_center = segment_2(1,:) - unit_v_lineSegment2_orthogonal*radius_of_arc;
        end_point_on_1 = circle_center + unit_v_lineSegment1_orthogonal*radius_of_arc;
    end
end

% Show the results?
if flag_do_debug == 1
    figure(fig_debug);


    % Show the circle center
    plot(circle_center(:,1), circle_center(:,2),'+','color',[1 0 0],'LineWidth',7,'MarkerSize',30);

    % Show the start and end points of segments
    plot(end_point_on_1(:,1), end_point_on_1(:,2),'.','color',[1 0 0],'LineWidth',7,'MarkerSize',30);
    plot(start_point_on_2(:,1), start_point_on_2(:,2),'.','color',[1 0 0],'LineWidth',7,'MarkerSize',30);
end

%% Find the start and end angles
% NOTE: this is vectorized. Comes from the geometry library.
if flag_point_to_left
    cross_product_direction = 1;
else
    cross_product_direction = -1;
end
[...
    angles] ...
    = ...
    fcn_geometry_findAngleUsing2PointsOnCircle(...
    circle_center,...
    radius_of_arc,...
    end_point_on_1,...
    start_point_on_2,...
    cross_product_direction);

vector_from_center_to_1 = end_point_on_1 - circle_center;
start_angle_radians = atan2(vector_from_center_to_1(1,2),vector_from_center_to_1(1,1));
end_angle_radians   = start_angle_radians + angles;
angles_to_calculate = linspace(start_angle_radians,end_angle_radians,number_of_points_on_curve)';

%% Find the points and line segments
arc_points = circle_center + radius_of_arc*[cos(angles_to_calculate) sin(angles_to_calculate)];

if flag_point_1_is_closest
    line_segment_1 = [nan nan];
    line_segment_2 = [start_point_on_2; segment_2(1,:)];

else
    line_segment_2 = [nan nan];
    line_segment_1 = [segment_1(2,:); end_point_on_1];
end


% Show the results?
if flag_do_debug == 1
    figure(fig_debug);


    % Show the arc_points
    plot(arc_points(:,1), arc_points(:,2),'.-','color',[1 0 1],'LineWidth',3,'MarkerSize',30);
    plot(arc_points(:,1), arc_points(:,2),'.-','color',[0 0 0],'LineWidth',3,'MarkerSize',15);

    % Show the start and end line segment
    plot(line_segment_1(:,1), line_segment_1(:,2),'.-','color',[0 1 1],'LineWidth',3,'MarkerSize',30);
    plot(line_segment_2(:,1), line_segment_2(:,2),'.-','color',[0 1 1],'LineWidth',3,'MarkerSize',30);
    plot(line_segment_1(:,1), line_segment_1(:,2),'.', 'color',[0 0 0],'LineWidth',3,'MarkerSize',15);
    plot(line_segment_2(:,1), line_segment_2(:,2),'.', 'color',[0 0 0],'LineWidth',3,'MarkerSize',15);

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
    figure(fig_num);
    clf;
    hold on;
    grid on;
    xlabel('X [m]');
    ylabel('Y [m]');
    axis equal;


    % Show the inputs
    plot(segment_1(:,1), segment_1(:,2),'.-','color',[0 1 0],'LineWidth',7,'MarkerSize',30);
    plot(segment_2(:,1), segment_2(:,2),'.-','color',[0 0 1],'LineWidth',5,'MarkerSize',30);


    % Show the directions
    quiver(segment_1(2,1), segment_1(2,2),unit_v_lineSegment1(1,1),unit_v_lineSegment1(1,2),0,'-','LineWidth',5,'ShowArrowHead','on','Color',[0 1 0],'MaxHeadSize',4)
    quiver(segment_2(1,1), segment_2(1,2),unit_v_lineSegment2(1,1),unit_v_lineSegment2(1,2),0,'-','LineWidth',3,'ShowArrowHead','on','Color',[0 0 1],'MaxHeadSize',4)

    quiver(segment_1(2,1), segment_1(2,2),unit_v_lineSegment1_orthogonal(1,1),unit_v_lineSegment1_orthogonal(1,2),0,'-','LineWidth',5,'ShowArrowHead','on','Color',[0 1 0],'MaxHeadSize',4)
    quiver(segment_2(1,1), segment_2(1,2),unit_v_lineSegment2_orthogonal(1,1),unit_v_lineSegment2_orthogonal(1,2),0,'-','LineWidth',3,'ShowArrowHead','on','Color',[0 0 1],'MaxHeadSize',4)

    plot(segment_1(:,1), segment_1(:,2),'.','color',[0 0 0],'LineWidth',7,'MarkerSize',15);
    plot(segment_2(:,1), segment_2(:,2),'.','color',[0 0 0],'LineWidth',5,'MarkerSize',15);

    % Show the circle center
    plot(circle_center(:,1), circle_center(:,2),'+','color',[1 0 0],'LineWidth',2,'MarkerSize',10);

    % Show the start and end points of segments
    plot(end_point_on_1(:,1), end_point_on_1(:,2),'.','color',[1 0 0],'LineWidth',7,'MarkerSize',30);
    plot(start_point_on_2(:,1), start_point_on_2(:,2),'.','color',[1 0 0],'LineWidth',7,'MarkerSize',30);

    % Show the arc_points
    plot(arc_points(:,1), arc_points(:,2),'.-','color',[1 0 1],'LineWidth',3,'MarkerSize',30);
    plot(arc_points(:,1), arc_points(:,2),'.-','color',[0 0 0],'LineWidth',3,'MarkerSize',15);

    % Show the start and end line segment
    plot(line_segment_1(:,1), line_segment_1(:,2),'.-','color',[0 1 1],'LineWidth',3,'MarkerSize',30);
    plot(line_segment_2(:,1), line_segment_2(:,2),'.-','color',[0 1 1],'LineWidth',3,'MarkerSize',30);
    plot(line_segment_1(:,1), line_segment_1(:,2),'.', 'color',[0 0 0],'LineWidth',3,'MarkerSize',15);
    plot(line_segment_2(:,1), line_segment_2(:,2),'.', 'color',[0 0 0],'LineWidth',3,'MarkerSize',15);

    % Show the start and end points of segments
    plot(line_segment_1(1,1), line_segment_1(1,2),'o','color',[0 1 0],'LineWidth',1,'MarkerSize',10);
    plot(line_segment_1(end,1), line_segment_1(end,2),'x','color',[1 0 0],'LineWidth',1,'MarkerSize',10);
    plot(line_segment_2(1,1), line_segment_2(1,2),'o','color',[0 1 0],'LineWidth',1,'MarkerSize',10);
    plot(line_segment_2(end,1), line_segment_2(end,2),'x','color',[1 0 0],'LineWidth',1,'MarkerSize',10);

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

%% Calculate cross products
function result = crossProduct(v,w)
result = v(:,1).*w(:,2)-v(:,2).*w(:,1);
end

%% fcn_INTERNAL_findAngleBetweenUnitVectors
function angle_between_vectors_radians = fcn_INTERNAL_findAngleBetweenUnitVectors(from_unit_vector,to_unit_vector)
cross_result = cross([from_unit_vector(1,1:2) 0],[to_unit_vector(1,1:2) 0]);
dot_result = sum(from_unit_vector(1,1:2).*to_unit_vector(1,1:2),2);
angle_magnitude = real(acos(dot_result));
angle_sign = sign(cross_result(3));
angle_between_vectors_radians = angle_sign*angle_magnitude;
end % Ends fcn_INTERNAL_findAngleBetweenUnitVectors

%% fcn_geometry_findAngleUsing2PointsOnCircle
function [...
    angles] ...
    = ...
    fcn_geometry_findAngleUsing2PointsOnCircle(...
    centers,...
    radii,...
    start_points_on_circle,...
    end_points_on_circle,...
    cross_products)

% fcn_geometry_findAngleUsing2PointsOnCircle -  This function calculates
% the angle from the start_points location to the end_points, in the
% direction of the vector given by is_clockwise.
%
% FORMAT:
%
% [angles] ...
%     = ...
%     fcn_geometry_findAngleUsing2PointsOnCircle(...
%     centers,...
%     radii,...
%     start_points_on_circle,...
%     end_points_on_circle,...
%     cross_products,...
%     varargin)
%
% INPUTS:
%
%      centers: an [N x 2] vector in [x y] of the points of circle centers
%
%      radii: a [N x 1] vector of the radii of the circles (to avoid
%      calculation time)
%
%      start_points_on_circle: an [N x 2] vector in [x y] of the points
%      where sectors start
%
%      end_points_on_circle: an [N x 2] vector in [x y] of the points
%      where sectors end
%
%      cross_products: an [N x 1] vector denoting the cross product
%      direction to follow from input point to output point
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      angles: an [N x 1] vector of the angles, in radians, between input
%      points and output points
%
% DEPENDENCIES:
%
%      fcn_geometry_checkInputsToFunctions
%      fcn_geometry_plotCircle
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_findAngleUsing2PointsOnCircle
% for a full test suite.
%
% This function was written on 2020_05_22 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Use number of radii to calculate the number of centers
num_circles = length(centers(:,1));

%% Step 1: calculate unit vectors for incoming and outgoing points
% Makes things easy later
unit_radial_to_inpoints = (start_points_on_circle - centers)./radii;
unit_radial_to_outpoints = (end_points_on_circle - centers)./radii;

%% Step 2: calculate the dot product angle from in and out unit vectors
dot_product = sum(unit_radial_to_inpoints.*unit_radial_to_outpoints,2);

% Correct the sign of the angle as well to match the cross product input
angles = acos(dot_product).*cross_products;

%% Step 3: calculate the cross products from in to out
cross_in_to_out = cross(...
    [unit_radial_to_inpoints, zeros(num_circles,1)],...
    [unit_radial_to_outpoints, zeros(num_circles,1)]);

%% Step 4: check if the cross product matches
% If the cross_in_to_out is in opposite direction from
% given cross products, then we need to take the refelex angle.
need_reflex_angles = (cross_in_to_out(:,3)...
    .*cross_products)<0;
angles(need_reflex_angles)=2*pi - angles(need_reflex_angles);

end % Ends fcn_geometry_findAngleUsing2PointsOnCircle


