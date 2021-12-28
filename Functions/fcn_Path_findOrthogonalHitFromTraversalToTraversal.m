function [closest_path_points,closest_distances] = ...
    fcn_Path_findOrthogonalHitFromTraversalToTraversal(...
    query_stations,...
    central_traversal,...
    nearby_traversal, ...
    varargin)
% fcn_Path_findOrthogonalHitFromTraversalToTraversal
% Given a central traversal and a set of stations along that traversal,
% finds the location on nearby traversals that are closest to the central
% traveral at each station point. Closest is defined via an orthogonal
% projection (or modifications of orthogonal projections) from the central
% traversal outward toward nearby traversals.
%
% FORMAT:
%
%      [closest_path_point,s_coordinate] = ...
%        fcn_Path_findOrthogonalHitFromTraversalToTraversal(...
%        query_stations,central_traversal,nearby_traversal,...
%        (flag_rounding_type),(search_radius),(fig_num));
%
% INPUTS:
%
%      query_stations: an N x 1 vector, with N>=1, containing the station
%      on the central traversal where the projections should take place
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
% DEPENDENCIES:
%
%      fcn_Path_checkInputsToFunctions
%      fcn_Path_findOrthogonalTraversalVectorsAtStations
%      fcn_Path_convertPathToTraversalStructure
%      fcn_Path_findTraversalStationSegment
%      fcn_Path_findProjectionHitOntoPath
%
% EXAMPLES:
%
% See the script: script_test_fcn_Path_findOrthogonalHitFromTraversalToTraversal
% for a full test suite.
%
% This function was written on 2020_11_14 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
%      2020_11_14:
%      -- first write of the code
%      2020_12_25:
%      -- changed plot style for situation where query overlaps the central
%      path (so have different line styles)
%      -- fixed a bug where the index was being used, instead of station, to
%      define a search area
%      -- isolated plotting functionality from debugging functionality
%      2021_01_07 
%      -- updated the name to fix path notation to traversal or traversals
%      2021_01_09
%      -- added input argument checking
%      -- converted all internal SXY variables to traversals
%      -- updated dependencies
%      2021_12_27
%      -- changed name to singular traversal since code only works with one
%      traversal at a time

% TO DO:
% Define search radius - need to know when a loop is going backward!

flag_do_debug = 0; % Flag to debug the results
flag_do_plot = 0; % Flag to plot the results
flag_check_inputs = 1; % Flag to perform input checking
flag_limit_station_range = 0; % Flag that limits the range of stations to search

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
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
% (station,central_traversal,nearby_traversal, (flag_projection_type?))


if flag_check_inputs == 1
    % Are there the right number of inputs?
    if nargin < 3 || nargin > 6
        error('Incorrect number of input arguments')
    end
    
    % Check the query_stations input
    fcn_Path_checkInputsToFunctions(query_stations, 'station');
    
    % Check the central_traversal input
    fcn_Path_checkInputsToFunctions(central_traversal, 'traversal');
    
    % Check the nearby_traversal input
    fcn_Path_checkInputsToFunctions(nearby_traversal, 'traversal');
        
    % Check cases specific to this function
    if any(query_stations<central_traversal.Station(1)) || any(query_stations>central_traversal.Station(end,1))
        error('The station query locations must be within the range of stations within the central_traversal');
    end
       
end

flag_rounding_type = 1;
% Does the user want to give a rounding type?
if 4 <= nargin
    flag_rounding_type = varargin{1};
end

% Define search radius?
search_radius = central_traversal.Station(end)*3;
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


Nstations = length(query_stations(:,1));

%% Find the unit normal vectors at each of the station points
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalTraversalVectorsAtStations(...
    query_stations,central_traversal,flag_rounding_type);
unit_vector_displacement = unit_normal_vector_end - unit_normal_vector_start;

%% Define the sensor vector from the unit vector
sensor_vector_start = unit_normal_vector_start;
postive_sensor_vector_end = unit_normal_vector_start + unit_vector_displacement*search_radius;
negative_sensor_vector_end = unit_normal_vector_start - unit_vector_displacement*search_radius;

%% Define the path we are searching from the nearby traversal
% Need to constrain this to be only within a small range of the s-distance?
path_to_check = [nearby_traversal.X nearby_traversal.Y];
traversal_to_check = fcn_Path_convertPathToTraversalStructure(path_to_check);

% Initialize the location of the hits to zeros
locations_of_hits = zeros(Nstations,2);
distances_of_hits = zeros(Nstations,1);

%% Loop through all the stations, finding the hit points at each station
for i_station=1:Nstations   
    
    % Check flag to see if we need to limit the station search range
    if flag_limit_station_range
        % Define the query station
        query_station = query_stations(i_station,1);

        % Define the search region by trimming it down?
        s_coord_start = query_station - search_radius;
        s_coord_end = query_station + search_radius;
        
        % Format: [pathSXY_segment,flag_outside_start, flag_outside_end] = ...
        [traversal_segment, ~, ~] = ...
            fcn_Path_findTraversalStationSegment(traversal_to_check, s_coord_start,s_coord_end);
        
        path_segment_to_check = [traversal_segment.X traversal_segment.Y];
    else
        path_segment_to_check = [traversal_to_check.X traversal_to_check.Y];
    end
    
    if flag_do_debug
        % Find results in the search region, plotting results
        [positive_distance,positive_location] = ...
            fcn_Path_findProjectionHitOntoPath(...
            path_segment_to_check,...
            sensor_vector_start(i_station,:),...
            postive_sensor_vector_end(i_station,:),0,222222);
        
        [negative_distance,negative_location] = ...
            fcn_Path_findProjectionHitOntoPath(...
            path_segment_to_check,...
            sensor_vector_start(i_station,:),...
            negative_sensor_vector_end(i_station,:),0,222222);
        
    else
        % Find results in the search region, no plotting
        [positive_distance,positive_location] = ...
            fcn_Path_findProjectionHitOntoPath(...
            path_segment_to_check,...
            sensor_vector_start(i_station,:),...
            postive_sensor_vector_end(i_station,:),0);
        [negative_distance,negative_location] = ...
            fcn_Path_findProjectionHitOntoPath(...
            path_segment_to_check,...
            sensor_vector_start(i_station,:),...
            negative_sensor_vector_end(i_station,:),0);
    end
    
    % Check that the distances are not NaN values
    if all(isnan(positive_distance))
        positive_location = [nan nan];
        positive_distance = nan;
    end
    if all(isnan(negative_distance))
        negative_location = [nan nan];
        negative_distance = nan;
    end

    % Fill in distances
    if ~isnan(positive_distance) && ~isnan(negative_distance)
        if positive_distance >= negative_distance
            locations_of_hits(i_station,:) = positive_location;
            distances_of_hits(i_station,1) = positive_distance;
        else
            locations_of_hits(i_station,:) = negative_location;
            distances_of_hits(i_station,1) = negative_distance;
        end
    elseif ~isnan(positive_distance) && isnan(negative_distance)
        locations_of_hits(i_station,:) = positive_location;
        distances_of_hits(i_station,1) = positive_distance;
    elseif isnan(positive_distance) && ~isnan(negative_distance)
        locations_of_hits(i_station,:) = negative_location;
        distances_of_hits(i_station,1) = negative_distance;
    elseif isnan(positive_distance) && isnan(negative_distance)
        locations_of_hits(i_station,:) = nan;
        distances_of_hits(i_station,1) = nan;
    else
        error('should not enter here!');
    end
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
       
    %     % Show the unit vectors
    %     normal_unit_vectors_at_stations = ...
    %         unit_normal_vector_end - unit_normal_vector_start;
    %     quiver(unit_normal_vector_start(:,1),unit_normal_vector_start(:,2),...
    %         normal_unit_vectors_at_stations(:,1),normal_unit_vectors_at_stations(:,2),0,'g','Linewidth',3);
    %
    %
    %     legend('Central traversal','Path to check','Station query points','Hit locations','Unit vectors');
    
    % Show the sensor vectors
    sensor_vectors_at_stations = ...
        search_radius*(unit_normal_vector_end - unit_normal_vector_start);
    quiver(unit_normal_vector_start(:,1),unit_normal_vector_start(:,2),...
        sensor_vectors_at_stations(:,1),sensor_vectors_at_stations(:,2),0,'g','Linewidth',3);
    
    
    legend('Central traversal','Path to check','Station query points','Hit locations','Sensor vectors');
    
    
    % % Label the points with distances?
    % for i_point = 1:length(path(:,1))
    %     text(path(i_point,2),path(i_point,3),sprintf('%.2f',distances_point_to_path(i_point)));
    % end
    
    
    
    
end % Ends the flag_do_plot if statement

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends the function

