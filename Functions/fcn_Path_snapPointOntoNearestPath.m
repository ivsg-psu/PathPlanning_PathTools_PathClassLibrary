function [closest_path_point,s_coordinate,...
    first_path_point_index,...
    second_path_point_index,...
    percent_along_length] = ...
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
%      [closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
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
%      (optional inputs) 
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
    if nargin < 2 || nargin > 3
        error('Incorrect number of input arguments')
    end
    
    % Check the data input
    fcn_Path_checkInputsToFunctions(path, 'path2or3D');
    
    % Check that the dimension of the point and path match
    if length(point(1,:)) ~= length(path(1,:))
        error('The dimension of the query point, in number of columns, must match the dimension of the path, in number of columnts');
    end
end

% Does user want to show the plots?
if 3 == nargin
    fig_num = varargin{1};
    figure(fig_num);
    flag_do_debug = 1;
else
    if flag_do_debug
        fig = figure;  %#ok<UNRCH>
        fig_num = fig.Number;
    end
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

% The solution method is as follows:
%  1. Find the closest point on the path to the query point
%  2. Find path points behind and ahead of the closest point from step 1
%   -> Check for end cases. For start/end, use the adjacent points.
%  3. Find percentage of travel on both the segments using dot products
%  4. Find projected point on the segment that has percentage of travel
%  between 0 and 1
%  5. If percentage of travel is not between 0 and 1 for both the segments,
%  then choose the projected point as closest point from step 1


traversal = fcn_Path_convertPathToTraversalStructure(path);
path_station = traversal.Station;

% Find distance from a point to every point on the path
distances_point_to_path = sum((path - point).^2,2).^0.5; % Calculate distances
[~,closest_path_point_index] = min(distances_point_to_path);  % Grab index of the closest point

% Be sure to check end cases as well
if 1 == closest_path_point_index || Npoints == closest_path_point_index
    if 1 == closest_path_point_index
        first_path_point_index  = 1;
        second_path_point_index = 2;
    else
        first_path_point_index  = Npoints-1;
        second_path_point_index = Npoints;
    end
    
    % Do the dot products - define the vectors first
    % See: https://mathinsight.org/dot_product for explanation
    % Basically, we are seeing what amount the point_vector points in the
    % direction of the path_vector
    path_vector  = path(second_path_point_index,:)-path(first_path_point_index,:);
    path_segment_length  = sum(path_vector.^2,2).^0.5;
    point_vector = point-path(first_path_point_index,:);
    projection_distance  = dot(path_vector,point_vector)/path_segment_length; % Do dot product
    percent_along_length = projection_distance/path_segment_length;
    
    % Calculate the outputs
    closest_path_point = path(first_path_point_index,:) + path_vector*percent_along_length;
    s_coordinate       = path_station(first_path_point_index,1) + path_segment_length*percent_along_length;
else
    % Do the dot products - define the vectors first
    % See: https://mathinsight.org/dot_product for explanation
    % Basically, we are seeing what amount the front_point_vector points in 
    % the direction of the front_path_vector
    front_path_vector  = path(closest_path_point_index+1,:)-path(closest_path_point_index,:);
    front_path_segment_length  = sum(front_path_vector.^2,2).^0.5;
    front_point_vector = point-path(closest_path_point_index,:);
    front_projection_distance  = dot(front_path_vector,front_point_vector)/front_path_segment_length; % Do dot product
    front_percent_along_length = front_projection_distance/front_path_segment_length;
    
    % Do the dot products - define the vectors first
    % See: https://mathinsight.org/dot_product for explanation
    % Basically, we are seeing what amount the back_point_vector points in 
    % the direction of the back_path_vector
    back_path_vector  = path(closest_path_point_index,:)-path(closest_path_point_index-1,:);
    back_path_segment_length  = sum(back_path_vector.^2,2).^0.5;
    back_point_vector = point-path(closest_path_point_index-1,:);
    back_projection_distance  = dot(back_path_vector,back_point_vector)/back_path_segment_length;    % Do dot product
    back_percent_along_length = back_projection_distance/back_path_segment_length;
    
    if 0 > back_percent_along_length % Point is located BEHIND the vector that is the rear-most vector
        error('ERROR: Point is lying BEHIND the BACK segment on path');
    elseif 1 < back_percent_along_length  % Point is in front of vector that is the rear-most vector
        if 1 < front_percent_along_length % Point is ahead of the vector that is in front
            error('ERROR: Point is lying AHEAD of the FRONT segment on path');
        elseif 0 > front_percent_along_length  % point is BEFORE start of front and AHEAD start of back - this is normal
            first_path_point_index  = closest_path_point_index;
            second_path_point_index = closest_path_point_index;
            percent_along_length    = 0;
            
            % Calculate the outputs
            closest_path_point = path(closest_path_point_index,:);
            s_coordinate       = path_station(closest_path_point_index,1);
        else  % Only way to enter here is if point is ahead of back, and ON the front vector
            first_path_point_index  = closest_path_point_index;
            second_path_point_index = closest_path_point_index+1;
            percent_along_length    = front_percent_along_length;
            
            % Calculate the outputs
            closest_path_point = path(first_path_point_index,:) + ...
                front_path_vector*front_percent_along_length;
            s_coordinate       = path_station(first_path_point_index,1) + ...
                front_path_segment_length*front_percent_along_length;
        end
    else % Only way to enter here is if point is on the back vector
        first_path_point_index  = closest_path_point_index-1;
        second_path_point_index = closest_path_point_index;
        percent_along_length    = back_percent_along_length;
        
        % Calculate the outputs
        closest_path_point = path(first_path_point_index,:) + ...
            back_path_vector*back_percent_along_length;
        s_coordinate       = path_station(first_path_point_index,1) + ...
            back_path_segment_length*back_percent_along_length;
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
if flag_do_debug
    figure(fig_num);
    hold on;
    grid on;
    
    % Is this a 2-D query
    if length(path(1,:))==2
        % Plot the path
        plot(path(:,1),path(:,2),'r-','Linewidth',5);
        plot(path(:,1),path(:,2),'ro','Markersize',20);
        
        axis equal;
        
        % Plot the query point
        plot(point(:,1),point(:,2),'ko');
        text(point(:,1),point(:,2),'Query point');
        
        % Plot the closest path points;
        plot(...
            path(first_path_point_index:second_path_point_index,1),...
            path(first_path_point_index:second_path_point_index,2),'r*');
        
        % % Label the points with distances?
        % for i_point = 1:length(path(:,1))
        %     text(path(i_point,2),path(i_point,3),sprintf('%.2f',distances_point_to_path(i_point)));
        % end
        
        % Plot the closest point on path
        plot(closest_path_point(:,1),closest_path_point(:,2),'go','Markersize',20);
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

