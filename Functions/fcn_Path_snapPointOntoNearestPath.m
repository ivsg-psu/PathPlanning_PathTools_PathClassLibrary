function [closest_path_point,...
    s_coordinate,...
    first_path_point_index,...
    second_path_point_index,...
    percent_along_length,...
    distance_real,...
    distance_imaginary] = ...
    fcn_Path_snapPointOntoNearestPath(point, path,varargin)
% fcn_Path_snapPointOntoNearestPath
% Finds location on a path that is closest to a given query point, e.g.
% snapping the point onto the path.
%
% The solution method is as follows:
%  1. Find the closest point on the path to the query point
%  2. Find path points behind and ahead of the closest point from step 1
%   -> Check for end cases. For start/end, use the adjacent points.
%  3. Find percentage of travel on both the segments using dot products
%  4. Find projected point on the segment that has percentage of travel
%  between 0 and 1
%  5. If percentage of travel is not between 0 and 1 for both the segments,
%  then choose the snap point as closest point, e.g. the result from step 1

% 
% FORMAT: 
%
%    [closest_path_point,...
%     s_coordinate,...
%     first_path_point_index,...
%     second_path_point_index,...
%     percent_along_length,...
%     distance_real,...
%     distance_imaginary] = ...
%      fcn_Path_snapPointOntoNearestPath(point, path,{fig_num})
%
% INPUTS:
%
%      point: a 1x2 vector containing the [X Y] location of the point or
%             a 1x3 vector containing the [X Y Z] location of the point
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
%          vectors.
%
%      figure_number: figure number where results are plotted
%
% OUTPUTS:
%
%      closest_path_point: a 1x2 or 1x3 vector containing the [X Y] or [X Y
%      Z] location of the nearest point on the path
%
%      s_coordinate: a scalar (1x1) representing the s-coordinate distance
%      along the path
%
%      first_path_point_index: a scalar (1x1) representing the index of the
%      starting "row" of the path where the closest point landed.
%
%      second_path_point_index: a scalar (1x1) representing the index of the
%      ending "row" of the path where the closest point landed.
%
%      percent_along_length: a scalar (1x1) of the percent along the
%      segment of the path where the closest point landed.
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

% TO-DO:
% Allow multiple points, e.g.
%      point: a Nx2 vector where N is the number of points, but at least 1. 

flag_do_debug = 0; % Flag to plot the results for debugging
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
% Are the input vectors the right shape?
Npoints = length(path(:,1));

if flag_check_inputs == 1
    % Are there the right number of inputs?
    narginchk(2,4);
    
    % Check the data input
    fcn_Path_checkInputsToFunctions(path, 'path2or3D');
    
    % Check that the dimension of the point and path match
    if length(point(1,:)) ~= length(path(1,:))
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
        fig_num=temp;
        figure(fig_num);
        flag_do_debug = 1;
    end
else
    if flag_do_debug
        fig = figure;  %#ok<UNRCH>
        fig_num = fig.Number;
    end
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

% The solution method is as follows:
%  1. Find the closest point on the path to the query point
%  2. Find path points behind and ahead of the closest point from step 1
%     2a. Check for end cases. For start/end, use the adjacent points.
%     2b. For other cases, just use ends
%  3. Find percentage of travel on both the segments using dot products
%  4. Find projected point on the segment that has percentage of travel
%  between 0 and 1
%  5. If percentage of travel is not between 0 and 1 for both the segments,
%  then choose the projected point as closest point from step 1


traversal = fcn_Path_convertPathToTraversalStructure(path);
path_station = traversal.Station;

% Step 1: find distance from a point to every point on the path
squared_distances_point_to_path = sum((path - point).^2,2); % Calculate distances
[~,closest_path_point_index] = min(squared_distances_point_to_path);  % Grab index of the closest point


% Step 2a: check if we are in one of the end cases
if 1 == closest_path_point_index || Npoints == closest_path_point_index
    % If in here, closest point snapped to the start or end of the path.
    % These are the end cases
    if 1 == closest_path_point_index % At the start?
        % TESTED2
        first_path_point_index  = 1;
        second_path_point_index = 2;
    else % At the end
        % TESTED2
        first_path_point_index  = Npoints-1;
        second_path_point_index = Npoints;
    end
    
    % Calculate the outputs using vector measures on the segment
    [closest_path_point,s_coordinate,percent_along_length,unit_orthogonal_projection_vector] = fcn_INTERNAL_calculateVectorMeasures(point,path(second_path_point_index,:),path(first_path_point_index,:),path_station(first_path_point_index,1));

else % Step 2b:
    % Calculate the outputs using vector measures on the front segment
    [front_closest_path_point,front_s_coordinate,front_percent_along_length,front_unit_orthogonal_projection_vector] = fcn_INTERNAL_calculateVectorMeasures(point,path(closest_path_point_index+1,:),path(closest_path_point_index,:),path_station(closest_path_point_index,1));  

    % Calculate the outputs using vector measures on the back segment   
    [back_closest_path_point,back_s_coordinate,back_percent_along_length,back_unit_orthogonal_projection_vector] = fcn_INTERNAL_calculateVectorMeasures(point,path(closest_path_point_index,:),path(closest_path_point_index-1,:),path_station(closest_path_point_index-1,1));
    
    if 0 > back_percent_along_length 
        % Point is located BEHIND the vector that is the rear-most vector -
        % not possible, as this should only happen if snap point is start
        % and that is caught in previous if statement
        error('ERROR: Point is lying BEHIND the BACK segment on path');
    elseif 1 < back_percent_along_length  
        % To enter here, point is after the end of the "back" vector's end
        % point.
        
        if 1 < front_percent_along_length 
            % Point is ahead of the vector that is in front - not possible
            % as this can only happen if snap point is the last point on
            % the path, and that is caught in previous if statement.
            error('ERROR: Point is lying AHEAD of the FRONT segment on path');
        
        elseif 0 > front_percent_along_length 
            % TESTED2
            % point is BEFORE start of front and also AFTER end of back 
            % This is the special situation when points are before the
            % "front" segment, e.g. the next path segment AND points are
            % after the previous "back" path segment. In this special case,
            % the distance calculation just uses the snap point.
            
            % Calculate the outputs
            first_path_point_index  = closest_path_point_index;
            second_path_point_index = closest_path_point_index;
            percent_along_length    = 0;
            closest_path_point = path(closest_path_point_index,:);
            s_coordinate       = path_station(closest_path_point_index,1);

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
                    unit_orthogonal_projection_vector = back_unit_orthogonal_projection_vector;
                case 2 % Use the front vector
                    unit_orthogonal_projection_vector = front_unit_orthogonal_projection_vector;
                case 3 % Use the average of the front and rear vectors
                    mixed_vector = back_unit_orthogonal_projection_vector+front_unit_orthogonal_projection_vector;
                    magnitude_mixed_vector = sum(mixed_vector.^2,2).^0.5;
                    unit_mixed_vector = mixed_vector/magnitude_mixed_vector;
                    unit_orthogonal_projection_vector = unit_mixed_vector;
                case 4
                    % This is a very special case, handled separately
                    error('Not coded yet!');
                otherwise
                    error('Unknown flag_rounding_type');
            end


        else
            % TESTED
            % Only way to enter here is if point is after the back segement,
            % and in the front vector area. This occurs when the point is
            % on the front vectors "main" area. We just use the front
            % segment result then.

            % Calculate the outputs
            first_path_point_index  = closest_path_point_index;
            second_path_point_index = closest_path_point_index+1;
            percent_along_length    = front_percent_along_length;            
            closest_path_point = front_closest_path_point;
            s_coordinate       = front_s_coordinate;
            unit_orthogonal_projection_vector = front_unit_orthogonal_projection_vector;
        end
    else
        % Only way to enter here is if point is IN the back vector between
        % 0 and 1. 
        
        if 0 > front_percent_along_length 
            % TESTED
            % Point is in back segment, and before front segment. Just use
            % the back segment results.
            
            % Calculate the outputs
            first_path_point_index  = closest_path_point_index-1;
            second_path_point_index = closest_path_point_index;
            percent_along_length    = back_percent_along_length;            
            closest_path_point = back_closest_path_point;
            s_coordinate       = back_s_coordinate;
            unit_orthogonal_projection_vector = back_unit_orthogonal_projection_vector;
        
        else
            % TESTED
            % Only way to enter here is if point is in the back segement,
            % and in front segment. This occurs when the point is
            % on the both the back and front vectors "main" area. We just
            % have to use the point that is closest to the test point.
            distance_squared_front_to_point = sum((point-front_closest_path_point).^2,2);
            distance_squared_back_to_point  = sum((point-back_closest_path_point).^2,2);
            
            % Calculate the outputs
            if distance_squared_front_to_point<=distance_squared_back_to_point
                % Front segment is closer
                first_path_point_index  = closest_path_point_index;
                second_path_point_index = closest_path_point_index+1;
                percent_along_length    = front_percent_along_length;
                closest_path_point = front_closest_path_point;
                s_coordinate       = front_s_coordinate;
                unit_orthogonal_projection_vector = front_unit_orthogonal_projection_vector;
            else
                % Back segment is closer
                first_path_point_index  = closest_path_point_index-1;
                second_path_point_index = closest_path_point_index;
                percent_along_length    = back_percent_along_length;
                closest_path_point = back_closest_path_point;
                s_coordinate       = back_s_coordinate;
                unit_orthogonal_projection_vector = back_unit_orthogonal_projection_vector;
            end
        end % Ends checks on front percent
    end % Ends checks on back percent
end % Ends check for if snap point is at ends



if (percent_along_length>0) && (percent_along_length<1)
    % Use the snap point itself, since it will be on the path.
    [distance_real, distance_imaginary] = ...
        fcn_INTERNAL_convertOffsetCoordinateToImaginaryNumber(unit_orthogonal_projection_vector,closest_path_point, point);
else % Use the point on the path closest to the test point
    [distance_real, distance_imaginary] = ...
    fcn_INTERNAL_convertOffsetCoordinateToImaginaryNumber(unit_orthogonal_projection_vector,path(closest_path_point_index,:), point);
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
if flag_do_debug
    figure(fig_num);
    clf;
    hold on;
    grid on;
    
    % Is this a 2-D query
    if length(path(1,:))==2
        % Plot the path
        plot(path(:,1),path(:,2),'r-','Linewidth',5);
        plot(path(:,1),path(:,2),'r.','Markersize',20);
        
        axis equal;
        
        % Plot the query point
        plot(point(:,1),point(:,2),'k.');
        text(point(:,1),point(:,2),'Query point');
        
        % Plot the closest path points;
        plot(...
            path(first_path_point_index:second_path_point_index,1),...
            path(first_path_point_index:second_path_point_index,2),'m.','MarkerSize',20);
        
        % Plot the reference vector
        quiver(path(closest_path_point_index,1),path(closest_path_point_index,2),...
            unit_orthogonal_projection_vector(1,1),unit_orthogonal_projection_vector(1,2),...
            0,'Color',[0.5 0.5 0.5]);
        midpoint = (path(closest_path_point_index,:) + path(closest_path_point_index,:)+unit_orthogonal_projection_vector)/2;
        text(midpoint(1,1),midpoint(1,2),'Reference vector','LineWidth',2,'VerticalAlignment','bottom');


        % Label the points with distances
        midpoint = (closest_path_point + point)/2;
        if distance_imaginary==0
            distance_string = sprintf('distance = %.2f',distance_real);
        else
            distance_string = sprintf('distance = %.2f + %.2f i',distance_real,distance_imaginary);
        end
        text(midpoint(1,1),midpoint(1,2),distance_string,'VerticalAlignment','top');
        
        % Plot the closest point on path
        plot(closest_path_point(:,1),closest_path_point(:,2),'g.','Markersize',20);
        text(closest_path_point(:,1),closest_path_point(:,2),'Snap Point on Path');
        
        
        % Connect closest point on path to query point
        plot(...
            [point(:,1) closest_path_point(:,1)],...
            [point(:,2) closest_path_point(:,2)],'g-','Linewidth',2);
        
    elseif length(path(1,:))==3
        
        % Plot the path
        plot3(path(:,1),path(:,2),path(:,3),'r-','Linewidth',5);
        plot3(path(:,1),path(:,2),path(:,3),'ro','Markersize',20);
        
        axis equal;
        
        % Plot the query point
        plot3(point(:,1),point(:,2),point(:,3),'ko');
        text( point(:,1),point(:,2),point(:,3),'Query point');
        
        % Plot the closest path points;
        plot3(...
            path(first_path_point_index:second_path_point_index,1),...
            path(first_path_point_index:second_path_point_index,2),...
            path(first_path_point_index:second_path_point_index,3),'r*');
        
        
        % Plot the closest point on path
        plot3(...
            closest_path_point(:,1),...
            closest_path_point(:,2),...
            closest_path_point(:,3),...
            'go','Markersize',20);
        text(...
            closest_path_point(:,1),...
            closest_path_point(:,2),...
            closest_path_point(:,3),...
            'Snap Point on Path');
        
        
        % Connect closest point on path to query point
        plot3(...
            [point(:,1) closest_path_point(:,1)],...
            [point(:,2) closest_path_point(:,2)],...
            [point(:,3) closest_path_point(:,3)],...
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
function     [closest_path_point,s_coordinate,percent_along_length,unit_orthogonal_projection_vector] = fcn_INTERNAL_calculateVectorMeasures(point,segment_end_point,segment_start_point,segment_start_station)
% Do the dot products - define the vectors first
% See: https://mathinsight.org/dot_product for explanation Basically,
% we are seeing what amount the point_vector points in the direction of
% the path_vector.


path_vector  = segment_end_point - segment_start_point;
path_segment_length  = sum(path_vector.^2,2).^0.5;
unit_orthogonal_projection_vector = path_vector/path_segment_length*[0 1; -1 0]; % Rotate by 90 degrees

start_to_point_vector = point-segment_start_point;
projection_distance  = dot(path_vector,start_to_point_vector)/path_segment_length; % Do dot product
% offset_distance  = dot(ortho_path_vector,start_to_point_vector)/path_segment_length; % Do dot product
percent_along_length = projection_distance/path_segment_length;

% Calculate the remaining outputs
closest_path_point = segment_start_point + path_vector*percent_along_length;
s_coordinate       = segment_start_station + path_segment_length*percent_along_length;
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

