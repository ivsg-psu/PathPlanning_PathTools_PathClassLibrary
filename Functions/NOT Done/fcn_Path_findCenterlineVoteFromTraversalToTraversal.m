function [centerline_points_projected,unit_vectors_orthogonal] = ...
    fcn_Path_findCenterlineVoteFromTraversalToTraversal(...
    from_traversal,to_traversal,varargin)

%% fcn_Path_findCenterlineVoteFromTraversalToTraversal
% Given a "from" traversal and a "to" traversal, the method is to
% orthogonally project from the "from" traversal to find the distance to
% the "to" traversal at each station in the "from" traversal. For
% situations where the projection does not hit anything, the nearest
% neighbor is used that has a hit. The projection distances are then
% divided in half to find the apparent centerline measured via the "from"
% traersal. As well, the orthogonal unit vectors for each projection are
% returned. The function allows the user to specify the flag_rounding_type
% and search_radius.
% 
% FORMAT: 
%
%    [centerline_points_projected,unit_vectors_orthogonal] = ...
%     fcn_Path_findCenterlineVoteFromTraversalToTraversal(...
%     from_traversal,to_traversal,(flag_rounding_type),(search_radius),(fig_num))
%
% INPUTS:
%
%      from_traversal: a traversal type wherein the projections are
%      calculated "from".
%
%      to_traversal: a traversal type wherein the projections are
%      calculated "to", to determine the distances "from".
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
%      flag_project_full_distance: this is a flag to determine whether the
%      projection is to the halfway distance between 1 and 2. For the
%      following values:
%
%          flag_project_full_distance = 0 (default): In this case, the returned path
%          is the halfway distance from 1 to 2, resulting in
%          the function returning the centerline between 1 and 2 as
%          measured from path 1.
%
%          flag_project_full_distance = 1: In this case, the returned path
%          is the full distance between 1 and 2, resulting in the function
%          returning the projection of path 1 onto 2.
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%      centerline_points_projected: a Mx2 or Mx3 vector containing the [X
%      Y (Z)] locations of the projected centerline, one for each station
%      point in the "from" traversal.
%
%      unit_vectors_orthogonal: vectors at each centerline point that are
%      unit vectors orthogonal to the projection used in the "from"
%      traversal. 
%
% EXAMPLES:
%      
% See the script: script_test_fcn_Path_findCenterlineVoteFromTraversalToTraversal
% for a full test suite.
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%     fcn_Path_findOrthogonalHitFromTraversalToTraversal
%     fcn_Path_findOrthogonalTraversalVectorsAtStations
%
% This function was written on 2023_09_04 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2023_09_04 by S. Brennan
% -- first write of the code
% 2023_09_06 by S. Brennan
% -- minor typo fix and better comments, including corrections in the
% header
% 2023_09_15 by S. Brennan
% -- added flag_project_full_distance to allow full distance projections

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
    narginchk(2,6);
    
    % fcn_DebugTools_checkInputsToFunctions
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

% Does user want to specify flag_project_full_distance?
flag_project_full_distance = 0; % Default
if 5 <= nargin
    temp = varargin{3};
    if ~isempty(temp)
        flag_project_full_distance=temp;
    end
end

% Does user want to show the plots?
if 6 == nargin
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

% Calculate the closest point and distances from right side to left
% FORMAT:
%
%      [closest_path_points,s_coordinate] = ...
%        fcn_Path_findOrthogonalHitFromTraversalToTraversal(...
%        query_stations,central_traversal,nearby_traversal,...
%        (flag_rounding_type),(search_radius),(fig_num));

[~,distances_between] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(...
    from_traversal.Station,...
    from_traversal,to_traversal,...
    flag_rounding_type,search_radius);

full_distances_between = fillmissing(distances_between, 'nearest');

% Find the centerline points:
% FORMAT: 
%
%      [unit_normal_vector_start, unit_normal_vector_end] = ...
%        fcn_Path_findOrthogonalTraversalVectorsAtStations(...
%        station_queries,central_traversal,...
%        (flag_rounding_type),(fig_num));
[unit_normal_vector_start, unit_normal_vector_end] = ...
       fcn_Path_findOrthogonalTraversalVectorsAtStations(...
       from_traversal.Station,...
       from_traversal,...
       flag_rounding_type);
unit_vectors = unit_normal_vector_end - unit_normal_vector_start;

if 0==flag_project_full_distance
    centerline_points_projected = [from_traversal.X from_traversal.Y]  + (unit_vectors/2).*full_distances_between;
else
    centerline_points_projected = [from_traversal.X from_traversal.Y]  + (unit_vectors).*full_distances_between;
end

unit_vectors_orthogonal = unit_vectors*[0 -1; 1 0];


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
    
    % Plot the input traversals
    plot(from_traversal.X,from_traversal.Y,'.-', 'Color',[0 1 0],'Linewidth',5,'MarkerSize',20);      
    plot(to_traversal.X,to_traversal.Y,'.-', 'Color',[0 0 1],'Linewidth',5,'MarkerSize',20);      

    % Plot the centerline_points_right_to_left and centerline_points_left_to_right
    plot(centerline_points_projected(:,1), centerline_points_projected(:,2), '.-','Linewidth',3,'MarkerSize',20,'Color',[1 0 0]*0.7);         

    % Show the unit vectors at the stations
    quiver(centerline_points_projected(:,1),centerline_points_projected(:,2),...
        unit_vectors_orthogonal(:,1),...
        unit_vectors_orthogonal(:,2),...
        0,'Color',[1 0 1],'Linewidth',2);

    legend('From path', 'To path','Calculated result','Unit vectors of calculated result');
    
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
