function [closest_path_points,closest_distances] = ...
    fcn_Path_FindOrthogonalHitFromPathToPath(stations,central_traversal,nearby_traversal, varargin)
% fcn_Path_FindOrthogonalHitFromPathToPath
% Given a central traversal and a set of stations along that traversal,
% finds the location on a nearby traversal that is closest to the central
% traveral at each station point. Closest is defined via an orthogonal
% projection (or modifications of orthogonal projections) from the central
% traversal outward toward nearby traversals.
%
% FORMAT:
%
%      [closest_path_point,s_coordinate] = ...
%        fcn_Path_FindOrthogonalHitFromPathToPath(...
%        stations,central_traversal,nearby_traversal,...
%        (flag_rounding_type),(search_radius),(fig_num));
%
% INPUTS:
%
%      stations: an N x 1 vector containing the station on the central
%      traversal where the projections should take place
%
%      central_traversal: a traversal structure that specifies the path
%      where projections to other paths are taking place.
%
%      nearby_traversal: a traversal structure that specifies the path
%      where intersections from the central_traversal would hit.
%
%      (OPTIONAL INPUTS)
%      flag_rounding_type: a flag to indicate which type of projection is
%      used, especially when stations are located at the end-points of
%      segments within the nearby_traversal. Note that the very first point
%      always uses projections from the following segement, and the very
%      last point always uses the prior. The flag determines behaviors for
%      endpoints of internal segments. The options include:
%
%          flag_rounding_type = 1;  % This is the default, and indicates that
%          the orthogonal projection of an endpoint is created by the PRIOR
%          segment leading up to each station query point.
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
%      search_radius: how far nearby the search points to look (default is
%      3 times the length of the station in the central_trajectory)
%
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      closest_path_point: a 1x2 vector containing the [X Y] location of
%      the nearest point on the path
%      s_coordinate: a scalar (1x1) representing the s-coordinate distance
%      along the path
%
% EXAMPLES:
%
% See the script: script_test_fcn_Path_FindOrthogonalHitFromPathToPath
% for a full test suite.
%
% This function was written on 2020_11_14 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
%     2020_11_14:
%     -- first write of the code
%     2020_12_25:
%     -- changed plot style for situation where query overlaps the central
%     path (so have different line styles)
%     -- fixed a bug where the index was being used, instead of station, to
%     define a search area
%     -- isolated plotting functionality from debugging functionality

% TO DO:
% Define search radius - need to let user define this as an input!

flag_do_debug = 0; % Flag to debug the results
flag_do_plot = 0; % Flag to plot the results
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'Starting function: %s, in file: %s\n',st(1).name,st(1).file);
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
% (station,central_traversal,nearby_traversal, (flag_projection_type?))

Nstations = length(stations(:,1));
if flag_check_inputs == 1
    % Are there the right number of inputs?
    if nargin < 3 || nargin > 6
        error('Incorrect number of input arguments')
    end
    
    % Check the station input
    if size(stations(1,:))~=1
        error('The station input must be an N x 1 vector');
    end
    
    % Check the central_traversal variables
    try
        X_central = central_traversal.X;
        Y_central = central_traversal.Y;
        Station_central = central_traversal.Station;
    catch
        error('The central_traversal input must be a structure containing the fields X, Y, and Station');
    end
    
    if length(X_central(:,1))~=length(Y_central(:,1))
        error('X and Y fields in the central traversal must have the same length, and must be N x 1 vectors');
    end
    
    if any(stations<Station_central(1)) || any(stations>Station_central(end,1))
        error('The station query locations must be within the range of stations within the central_traversal');
    end
    
    if ~issorted(Station_central,'strictascend')
        error('The station field on the central traversal must be increasing and contain no duplicates');
    end
    
    % Check the nearby_traversal variables
    try
        X_nearby = nearby_traversal.X;
        Y_nearby = nearby_traversal.Y;
        Station_nearby = nearby_traversal.Station;
    catch
        error('The nearby_traversal input must be a structure containing the fields X, Y, and Station');
    end
    
    if length(X_nearby(:,1))~=length(Y_nearby(:,1))
        error('X and Y fields in the nearby_traversal must have the same length, and must be N x 1 vectors');
    end
    
    if ~issorted(Station_nearby,'strictascend')
        error('The station field on the nearby_traversal must be increasing and contain no duplicates');
    end
    
end

flag_rounding_type = 1;
% Does the user want to give a rounding type?
if 4 <= nargin
    flag_rounding_type = varargin{1};
end

% Define search radius?
search_radius = max(Station_central)*3;
if 5 <= nargin
    search_radius = varargin{2};
end


% Does user want to show the plots?
if 6 == nargin
    fig_num = varargin{3};
    figure(fig_num);
    flag_do_plot = 1;
else
    if flag_do_debug
        fig = figure;
        fig_num = fig.Number;
        flag_do_plot = 1;
    end
end

%% Start of main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Find the unit normal vectors at each of the station points
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_FindOrthogonalPathVectorsAtStations(...
    stations,central_traversal,flag_rounding_type);


%% Define the path we are searching from the nearby traversal
% Need to constrain this to be only within a small range of the s-distance?
path_to_check = [nearby_traversal.X nearby_traversal.Y];
path_to_check_SXY = fcn_Path_convertXYtoSXY(path_to_check(:,1),path_to_check(:,2));

% Initialize the location of the hits to zeros
locations_of_hits = zeros(Nstations,2);
distances_of_hits = zeros(Nstations,1);

%% Loop through all the stations, finding the hit points at each station
for i_station=1:Nstations
    
    % Defi(ne the query station
    query_station = stations(i_station,1);
    
    % Define the search region by trimming it down?
    s_coord_start = query_station - search_radius;
    s_coord_end = query_station + search_radius;
    
    % Format: [pathSXY_segment,flag_outside_start, flag_outside_end] = ...
    [pathSXY_segment, ~, ~] = ...
        fcn_Path_findPathSXYSegment(path_to_check_SXY, s_coord_start,s_coord_end);
    path_segment_to_check = pathSXY_segment(:,2:3);
    
    % Find results in the search region
    [distance,location] = ...
        fcn_Path_findProjectionHitOntoPath(...
        path_segment_to_check,...
        unit_normal_vector_start(i_station,:),...
        unit_normal_vector_end(i_station,:),1);
    
    % Check that the distances are not NaN values
    if all(isnan(distance))
        location = [nan nan];
        distance = nan;
    end
    
    % Save results
    locations_of_hits(i_station,:) = location;
    distances_of_hits(i_station,1) = distance;
end


closest_path_points = locations_of_hits;
closest_distances = distances_of_hits;



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
if flag_do_plot
    figure(fig_num);
    hold on;
    grid on;
    
    % Plot the central traversal
    plot(central_traversal.X,central_traversal.Y,'k','Linewidth',3);
    
    % Plot the path
    % plot(path_to_check(:,1),path_to_check(:,2),'bo-','Linewidth',2);
    plot(path_to_check(:,1),path_to_check(:,2),'o-','Linewidth',2);
    
    
    axis equal;
    
    % Plot the station points
    plot(unit_normal_vector_start(:,1),unit_normal_vector_start(:,2),'k.','Markersize',35);
       
    % Plot hit locations
    plot(locations_of_hits(:,1),locations_of_hits(:,2),'r.','Markersize',25);
       
    % Show the unit vectors
    normal_unit_vectors_at_stations = ...
        unit_normal_vector_end - unit_normal_vector_start;
    quiver(unit_normal_vector_start(:,1),unit_normal_vector_start(:,2),...
        normal_unit_vectors_at_stations(:,1),normal_unit_vectors_at_stations(:,2),'g','Linewidth',3);
    
    
    legend('Central traversal','Path to check','Station query points','Hit locations','Unit vectors');
    
    % % Label the points with distances?
    % for i_point = 1:length(path(:,1))
    %     text(path(i_point,2),path(i_point,3),sprintf('%.2f',distances_point_to_path(i_point)));
    % end
    
    
    
    
end % Ends the flag_do_plot if statement

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends the function

