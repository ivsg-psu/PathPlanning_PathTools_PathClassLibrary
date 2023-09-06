function center_path = ...
    fcn_Path_findCenterPathBetweenTwoPaths(...
    first_path,second_path,varargin)

%% fcn_Path_findCenterPathBetweenTwoPaths
% Given two paths, finds the path that is geometrically inbetween both.
% 
% The method is to orthogonally project from the first_path to find
% the distance to the second_path at each station in the
% first_path. For situations where the projection does not hit
% anything, the nearest neighbor is used that has a hit. The projection
% distances are then divided in half to find the apparent centerline
% measured via the first_path. The process is then repeated,
% projecting from the second_path to the first_path, again
% finding the centerline measured from the second_path. For each
% centerline point, the orthogonal vector is kept for that centerline
% point. The resulting points are then sorted using the points and the
% orthogonal vectors, taking the first point as the start point in either
% the first_path or second_path that has no point behind it as
% measured via a dot-product with the centerline vector. Every point
% thereafter is sorted by the projections to create a composite centerline.
%
% FORMAT: 
%
%    center_path = ...
%     fcn_Path_findCenterPathBetweenTwoPaths(...
%     first_path,second_path,(flag_rounding_type),(search_radius),(fig_num))
%
% INPUTS:
%
%      first_path: the first traversal defining one edge of the
%      boundary.
%
%      second_path: the second traversal defining the other edge of
%      the boundary.
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
%      search_radius: the distance to project "from" to search for
%      intersections with the "to" traversal (default is 10 meters).
%
%      fig_num: figure number where results are plotted
%
% OUTPUTS:
%
%      center_path: a Mx2 or Mx3 vector containing the [X Y (Z)] points
%      that bisect the two traversals.
%
% EXAMPLES:
%      
% See the script: script_test_fcn_Path_findCenterPathBetweenTwoPaths
% for a full test suite.
%
% DEPENDENCIES:
%
%     fcn_Path_checkInputsToFunctions
%     fcn_Path_findOrthogonalHitFromTraversalToTraversal
%     fcn_Path_findOrthogonalTraversalVectorsAtStations
%     fcn_Path_findCenterlineVoteFromTraversalToTraversal
%
% This function was written on 2023_09_04 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2023_09_04 by S. Brennan
% -- first write of the code

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
    narginchk(2,5);
    
end


% Does user want to specify the rounding type?
flag_rounding_type = 1; % Default
if 3 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        flag_rounding_type=temp;
    end
end

% Does user want to specify the search_radius?
search_radius = 10; % Default
if 4 <= nargin
    temp = varargin{2};
    if ~isempty(temp)
        search_radius=temp;
    end
end

% Does user want to show the plots?
if 5 == nargin
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
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

% Convert the paths into traversal types
first_traversal   = fcn_Path_convertPathToTraversalStructure(first_path);
second_traversal  = fcn_Path_convertPathToTraversalStructure(second_path);

%% Obtain the projections from the traversals toward each other
[centerline_points_first_to_second,centerline_points_first_to_second_unit_vectors_orthogonal] = ...
    fcn_Path_findCenterlineVoteFromTraversalToTraversal(...
    first_traversal,second_traversal,...
    (flag_rounding_type),(search_radius));

[centerline_points_second_to_first,centerline_points_second_to_first_unit_vectors_orthogonal] = ...
    fcn_Path_findCenterlineVoteFromTraversalToTraversal(...
    second_traversal,first_traversal,...
    (flag_rounding_type),(search_radius));



%% Find which point is "first" in the sorting order
startfirst_to_startsecond_vector = second_path(1,1:2) - first_path(1,1:2);
% Check dot product - positive means that right side is first, negative
% means left
first_path_starts = sum(startfirst_to_startsecond_vector.*centerline_points_first_to_second_unit_vectors_orthogonal(1,:),2)>0;

N_firstPath = length(first_path(:,1));
N_secondPath  = length(second_path(:,1));
N_centerPath = N_firstPath+N_secondPath;

% Initialize all the points that need to be sorted
all_points = [centerline_points_first_to_second(:,1:2); centerline_points_second_to_first(:,1:2)];
all_unit_vectors = [centerline_points_first_to_second_unit_vectors_orthogonal; centerline_points_second_to_first_unit_vectors_orthogonal];
all_points_remaining = ones(1,N_centerPath);

% Tag the first point as "not remaining"
if first_path_starts
    first_index = 1;
else
    first_index = N_firstPath+1;
end

%% Sort all the centerline points to create a center path

% Initialize the first point
next_point = all_points(first_index,:);
next_orthoginal_vector = all_unit_vectors(first_index,:);
all_points_remaining(first_index) = 0;


% Initialize the center path
center_path = nan(N_centerPath,2);
current_point_count = 1;
center_path(current_point_count,:) = next_point;
previous_distance = inf;

% Loop through the points
while any(all_points_remaining==1)
    % Set the points to search
    points_to_search = all_points;
    points_to_search(all_points_remaining==0,1) = nan; % Make one value equal to NaN, which makes the calculations result in Nan that follow

    % Create a vector from current point to all remaining points
    vectors_to_search_points =  points_to_search - next_point;

    % Do dot products 
    dot_products = sum((next_orthoginal_vector.*vectors_to_search_points),2);

    % Check that we don't have any negatives (this is an error)
    if any(dot_products<-1*previous_distance)
        error('A dot product was found further negative than the previous distance - this should never happen!');
    end

    % Find the smallest distance, and keep the index for that point
    [~,closest_index] = min(dot_products);

    % Save the result
    all_points_remaining(closest_index) = 0;
    next_point = all_points(closest_index,:);
    current_point_count = current_point_count+1;
    center_path(current_point_count,:) = next_point;
    next_orthoginal_vector = all_unit_vectors(closest_index,:);
    previous_distance = dot_products(closest_index);

    % For debugging:
    % fprintf(1,'Closest index: %.0d at distance of %f\n',closest_index,previous_distance);
end

% Finally, check direcitonality - are all the vectors aligned with each
% other?
centerline_final_vectors = center_path(2:end,:)-center_path(1:end-1,:);
mag_centerline_final_vectors = sum(centerline_final_vectors.^2,2).^0.5;
unit_centerline_final_vectors = centerline_final_vectors./mag_centerline_final_vectors;

successive_dot_products = sum((unit_centerline_final_vectors(1:end-1,:).*unit_centerline_final_vectors(2:end,:)),2);
if any(successive_dot_products<0)
    error('Misaligned vectors found - this indicates an error');
else
    angles = acos(successive_dot_products)*180/pi;
    fprintf('Mean misalignment in angle is: %.3f degrees.\n',mean(angles));
    fprintf('Standard deviation in misalignment in angle is: %.3f degrees.\n',std(angles));
    fprintf('Max misalignment angle is: %.3f degrees.\n',max(angles));
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
    
    % Plot the input traversals
    plot(first_path(:,1),first_path(:,2),'k.-','Linewidth',3,'MarkerSize',20);      
    plot(second_path(:,1),second_path(:,2),'k.-','Linewidth',3,'MarkerSize',20);      

    % Plot the results
    plot_color = [0 1 0];
    line_width = 3;

    sizeOfMarkers = 10;
    plot(center_path(:,1),center_path(:,2), '-','Color',plot_color,'Linewidth',line_width,'Markersize',sizeOfMarkers);
    plot(center_path(:,1),center_path(:,2), '.','Color',plot_color,'Linewidth',line_width,'Markersize',round(sizeOfMarkers*2));
    plot(center_path(:,1),center_path(:,2), '.','Color',[0 0 0],'Linewidth',line_width,'Markersize',round(sizeOfMarkers*0.5));


    % Make axis slightly larger
    temp = axis;
    axis_range_x = temp(2)-temp(1);
    axis_range_y = temp(4)-temp(3);
    percent_larger = 0.3;
    axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);

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
