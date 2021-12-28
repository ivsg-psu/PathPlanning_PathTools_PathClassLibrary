function [unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalTraversalVectorsAtStations(station_queries,central_traversal, varargin)

% fcn_Path_findOrthogonalTraversalVectorsAtStations
% Given a central traversal and a set of stations along that traversal,
% finds the unit normal vector on the central traveral at each station
% point.
% 
% FORMAT: 
%
%      [unit_normal_vector_start, unit_normal_vector_end] = ...
%        fcn_Path_findOrthogonalTraversalVectorsAtStations(...
%        station_queries,central_traversal,...
%        (flag_rounding_type),(fig_num));
%
% INPUTS:
%
%      station_queries: an N x 1 vector containing the station on the
%      central traversal where the projections should take place
%
%      central_traversal: a traversal structure that specifies the path
%      where projections to other paths are taking place.
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
%      fig_num: a figure number to plot results. Turns debugging on.
%
% OUTPUTS:
%
%      unit_normal_vector start: a Nx2 vector containing [X1 Y1]
%      coordinates as columns, where the [X1 Y1] represents the location of
%      the start point of the vector, on the path.
%
%      unit_normal_vector_end: a Nx2 vector containing the [X2 Y2] location
%      of the end point of the unit vector. On both outputs, there are N
%      rows, one row for each station.
%
% DEPENDENCIES:
%
%      fcn_Path_checkInputsToFunctions
%
% EXAMPLES:
%      
% See the script: script_test_fcn_Path_findOrthogonalTraversalVectorsAtStations
% for a full test suite.
%
% This function was written on 2020_12_31 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%      2020_12_31:
%      -- first write of the code via modification from 
%      fcn_Path_FindOrthogonalHitFromPathToPath
%      2021_01_07
%      -- renamed to transition from path to traversal notation 
%      2021_12_27:
%      -- corrected dependencies in comments

% TO DO:
% Define search radius - need to let user define this as an input!

flag_do_debug = 0; % Flag to debug the results for debugging
flag_do_plots = 0; % Flag to create plots
flag_check_inputs = 1; % Flag to perform input checking

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
    if nargin < 2 || nargin > 4
        error('Incorrect number of input arguments')
    end
    
    % Check the station_queries input
    fcn_Path_checkInputsToFunctions(station_queries, 'station');
    
    % Check the central_traversal input    
    fcn_Path_checkInputsToFunctions(central_traversal, 'traversal');
    
    if any(station_queries<central_traversal.Station(1)) || any(station_queries>central_traversal.Station(end))
        error('The station query locations must be within the range of stations within the central_traversal');
    end
    
    if ~issorted(central_traversal.Station,'strictascend')
        error('The station field on the central traversal must be increasing and contain no duplicates');        
    end
    
end


flag_rounding_type = 1; % Set default value


% Does the user want to give a different rounding type?
if 3 <= nargin
    flag_rounding_type = varargin{1};
end

% Does user want to show the plots?
if 4 == nargin
    fig_num = varargin{2};
    figure(fig_num);
    flag_do_plots = 1;
else
    if flag_do_debug
        flag_do_plots = 1;
        fig = figure; 
        fig_num = fig.Number;
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

Station_central = central_traversal.Station;
X_central       = central_traversal.X;
Y_central       = central_traversal.Y;


%% Find the station coordinates on the central trajectory 
% Find the points on the central trajectory corresponding to the station
% query. Format for interp1: Vq = interp1(X,V,Xq). The result are X and Y
% locations that are ON the central path defined by X_central and
% Y_central, yet at the station distances given by the query.
X_central_at_stations = interp1(Station_central,X_central,station_queries,'linear');
Y_central_at_stations = interp1(Station_central,Y_central,station_queries,'linear');

%% Find the projection vector at the station coordinates
% Now need to create the projection vector for each station. This will
% depend on the input options. For options where the orthogonal projection
% only uses the prior or subsequent values, we can use the original
% Station_center indices. For options where averaging is occuring, we use
% intermediate vectors. To do this, we change the definition of
% Station_central if needed.
% 
% Here's the rounding options:
%      flag_rounding_type = 1;  % This is the default, and indicates that
%      the orthogonal projection of an endpoint is created by the PRIOR
%      segment.
%
%      flag_rounding_type = 2;  % This indicates that the orthogonal
%      projection of an endpoint is created by the FOLLOWING segment.
%
%      flag_rounding_type = 3;  % This indicates that the orthogonal
%      projection, ONLY at an endpoint is created by averaging both the
%      PRIOR segment and FOLLOWING segment.
%
%      flag_rounding_type = 4;  % This indicates that the orthogonal
%      projections along segments should be calculated at the midpoints of
%      each segment, and the endpoints of segments are an average of the
%      PRIOR and SUBSEQUENT midpoints. All projections are interpolated
%      from their prior and subsequent vectors.

% Figure out the number of indices within the central path, and number them
% from 1 to the length. Then find which indices would likely be closest to
% the query stations. (NOTE: this works because station distances increase
% linearly within each segment).
Indices_central = (1:length(Station_central))';

% Format for interp1: Vq = interp1(X,V,Xq)
Indices_Central_at_Query_Stations = interp1(Station_central, Indices_central,station_queries);

% Find the indices of start and end of the vectors
indices_start = floor(Indices_Central_at_Query_Stations);
indices_end   = ceil(Indices_Central_at_Query_Stations);

% If the start and end indices match, then the station coordinate is on top
% of an endpoint of a segment. In these cases, it can be unclear which
% projection vector to use, and this depends on the rounding type.
indices_on_endpoints = find(indices_start==indices_end);
if ~isempty(indices_on_endpoints)    
    switch flag_rounding_type
        case 1  % Default - use orthogonal projection via prior segment
            indices_start(indices_on_endpoints) = indices_end(indices_on_endpoints) - 1;
        case 2  % Use orthogonal projection of subsequent segment
            indices_end(indices_on_endpoints) = indices_start(indices_on_endpoints) + 1;
        case 3  % Average the two segments but only at the endpoints
            indices_end(indices_on_endpoints) = indices_start(indices_on_endpoints) + 0.5;
            indices_start(indices_on_endpoints) = indices_end(indices_on_endpoints) - 1;
        case 4
            % Do nothing - this is a special case which is filled in below
        otherwise
            error('Unrecognized method in flag_rounding_type.');
    end
end
   
if flag_rounding_type == 4 % Always do averaging
    indices_start = Indices_Central_at_Query_Stations - 0.5;
    indices_end   = Indices_Central_at_Query_Stations + 0.5;
end

%% Make sure all the indices are within acceptable range
indices_start = max(1,indices_start);
indices_end   = min(Indices_central(end),indices_end);


%% Check specifically the start and end locations
start_query_index = find(station_queries==Station_central(1));
if ~isempty(start_query_index)
    indices_start(start_query_index) = 1;
    indices_end(start_query_index) = 2;
end

end_query_index = find(station_queries==Station_central(end));
if ~isempty(end_query_index)
    indices_start(end_query_index) = Indices_central(end) - 1;
    indices_end(end_query_index) = Indices_central(end);
end


%% Create the vectors in the direction of the central trajectory 
% at the station points
X_central_indices_start = interp1(Indices_central, X_central,indices_start,'linear');
Y_central_indices_start = interp1(Indices_central, Y_central,indices_start,'linear');
X_central_indices_end   = interp1(Indices_central, X_central,indices_end,  'linear');
Y_central_indices_end   = interp1(Indices_central, Y_central,indices_end,  'linear');


%% Convert the tangent vectors to normal unit vectors
delta_X_central_at_stations = X_central_indices_end - X_central_indices_start;
delta_Y_central_at_stations = Y_central_indices_end - Y_central_indices_start;
tangent_vectors_at_stations = [delta_X_central_at_stations delta_Y_central_at_stations];  
magnitudes = sum(tangent_vectors_at_stations.^2,2).^0.5;
tangent_unit_vectors_at_stations = tangent_vectors_at_stations./magnitudes;
normal_unit_vectors_at_stations = tangent_unit_vectors_at_stations*[0 1; -1 0];

%% Create unit vectors
unit_normal_vector_start = [X_central_at_stations Y_central_at_stations]; 
unit_normal_vector_end   = [X_central_at_stations Y_central_at_stations] + normal_unit_vectors_at_stations;


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
    hold on;
    grid on;
    
    % Plot the central traversal
    plot(central_traversal.X,central_traversal.Y,'k','Linewidth',3);      
    
    % Make the axis normal shaped
    axis equal;
    
    % Plot the station points
    plot(X_central_at_stations,Y_central_at_stations,'r.','Markersize',30);
            
    % Show the unit vectors
    quiver(X_central_at_stations,Y_central_at_stations,normal_unit_vectors_at_stations(:,1),normal_unit_vectors_at_stations(:,2),0,'g','Linewidth',3);

    % Add a legend
    legend('Central traversal','Station query points','Unit vectors');  
    
end % Ends the flag_do_plots if statement

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends the function

