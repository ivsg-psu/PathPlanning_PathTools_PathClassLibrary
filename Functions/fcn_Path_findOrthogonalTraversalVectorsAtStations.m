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
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
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
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_Path_findPathOrthogonalVectors
%
% EXAMPLES:
%      
% See the script: script_test_fcn_Path_findOrthogonalTraversalVectorsAtStations
% for a full test suite.
%
% This function was written on 2020_12_31 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2020_12_31:
% -- first write of the code via modification from
% fcn_Path_FindOrthogonalHitFromPathToPath
% 2021_01_07
% -- renamed to transition from path to traversal notation
% 2021_12_27:
% -- corrected dependencies in comments
% 2023_04_29:
% -- added capability for interpolated results at endpoints that are
% undefined, using imaginary inputs. See flag type 5.
% 2023_08_27:
% -- fixed bug when many points are filled with the first or last
% vector, when previous version of code could only handle one point
% -- Cleaned up input checking a bit, allowing empty inputs and using
% narginchk
% 2025_06_23 - S. Brennan
% -- Updated debugging and input checks

% TO-DO
% Define search radius - need to let user define this as an input!


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

        % Check the station_queries input
        fcn_DebugTools_checkInputsToFunctions(station_queries, 'station');

        % Check the central_traversal input
        fcn_DebugTools_checkInputsToFunctions(central_traversal, 'traversal');

        if any(station_queries<central_traversal.Station(1)) || any(station_queries>central_traversal.Station(end))
            error('The station query locations must be within the range of stations within the central_traversal');
        end

        if ~issorted(central_traversal.Station,'strictascend')
            error('The station field on the central traversal must be increasing and contain no duplicates');
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
        fig_debug = 818181; %#ok<NASGU>
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

Station_central = central_traversal.Station;
X_central       = central_traversal.X;
Y_central       = central_traversal.Y;


%% Find the midpoint and joint tangent vectors for all segments 
% Call the function fcn_Path_findPathOrthogonalVectors

[normal_unit_vectors_at_midpoints, normal_unit_vectors_at_joints] = ...
    fcn_Path_findPathOrthogonalVectors([X_central Y_central],flag_rounding_type);

tangent_unit_vectors_at_midpoints = normal_unit_vectors_at_midpoints*[0 -1; 1 0];


%% Find the XY of query station coordinates on the central trajectory 
% Find the points on the central trajectory corresponding to the station
% query. Format for interp1: Vq = interp1(X,V,Xq). The result are X and Y
% locations that are ON the central path defined by X_central and
% Y_central, yet at the station distances given by the query.

X_central_at_stations = interp1(Station_central,X_central,station_queries,'linear');
Y_central_at_stations = interp1(Station_central,Y_central,station_queries,'linear');

%% Find the projection vector at these station coordinates


% Figure out the number of indices within the central path, and number them
% from 1 to the length. Then find which indices would likely be closest to
% the query stations. (NOTE: this works because station distances increase
% linearly within each segment).
Indices_central = (1:length(Station_central))';

% Format for interp1: Vq = interp1(X,V,Xq)
Indices_Central_at_Query_Stations = interp1(Station_central, Indices_central,station_queries);


%% Find the indices of start and end of the vectors
indices_start = floor(Indices_Central_at_Query_Stations);
indices_end   = ceil(Indices_Central_at_Query_Stations);


% Make sure all the indices are within acceptable range
indices_start = max(1,indices_start);
indices_end   = min(Indices_central(end),indices_end);

%% Initialize normal vectors
normal_unit_vectors_at_stations = nan(length(station_queries), 2); 


%% Check specifically the start and end locations
% Are any queries EXACTLY at the start location?
start_query_index = find(station_queries==Station_central(1));
if ~isempty(start_query_index)
    normal_unit_vectors_at_stations(start_query_index,:) = ones(length(start_query_index),1)*normal_unit_vectors_at_midpoints(1,:);
end

% Are any queries EXACTLY at the end location?
end_query_index = find(station_queries==Station_central(end));
if ~isempty(end_query_index)
    normal_unit_vectors_at_stations(end_query_index,:) = ones(length(end_query_index),1)*normal_unit_vectors_at_midpoints(end,:);
end


%% Check for queries on joints
% If the start and end indices match, then the station coordinate is on top
% of a "joint" of two segment. In these cases, it can be unclear which
% projection vector to use, and this depends on the rounding type.
% NOTE: must add 1 to the indicies as the first joint is the start point,
% and not an actual joint

indices_on_joints = find(indices_start==indices_end);
if ~isempty(indices_on_joints)    
    normal_unit_vectors_at_stations(indices_on_joints,:) = normal_unit_vectors_at_joints(indices_start(indices_on_joints),:);
end

%% Fill any that were not one of the above cases
indicies_not_on_joints = isnan(normal_unit_vectors_at_stations(:,1));
normal_unit_vectors_at_stations(indicies_not_on_joints,:) = normal_unit_vectors_at_midpoints(indices_start(indicies_not_on_joints),:);

%% Check for special case for flag type of 4
if flag_rounding_type == 4 % Always do averaging, but need to calculate ratios
    remainder = Indices_Central_at_Query_Stations - floor(Indices_Central_at_Query_Stations);

    % For this special situation, the joint vectors at the ends need to be
    % modified
    normal_unit_vectors_at_joints(1,:) = -1*normal_unit_vectors_at_midpoints(1,:)*[0 -1; 1 0];
    normal_unit_vectors_at_joints(end,:) =  normal_unit_vectors_at_midpoints(end,:)*[0 -1; 1 0];

   
    % Are any at the very end? If so, the floor operation gives the wrong
    % answers.
    segment_numbers = floor(Indices_Central_at_Query_Stations);
    if any(segment_numbers==Indices_central(end))
        indices_at_ends = segment_numbers==Indices_central(end);
        remainder(indices_at_ends) = 1;
        segment_numbers(indices_at_ends) = Indices_central(end-1);
    end

    % Fraction is 0 if at start of segment, if rounding down, 
    % and 1 if at midpoint
    round_down_indices = find(remainder<=0.5);    
    if ~isempty(round_down_indices)
        fraction_to_midpoint = 2*remainder(round_down_indices);
        normal_unit_vectors_at_stations(round_down_indices,:)   = ...
            (fraction_to_midpoint).*normal_unit_vectors_at_midpoints(segment_numbers(round_down_indices),:) ...
            + (1 - fraction_to_midpoint).*normal_unit_vectors_at_joints(segment_numbers(round_down_indices),:);
        mags = sum(normal_unit_vectors_at_stations(round_down_indices,:).^2,2).^0.5;
        normal_unit_vectors_at_stations(round_down_indices,:) = normal_unit_vectors_at_stations(round_down_indices,:)./mags;
    end

    % Fraction is 0 if at end of segment, if rounding up, 
    % and 1 if at midpoint
    round_up_indicies = find(remainder>0.5);
    if ~isempty(round_up_indicies)
        fraction_to_midpoint = 2*(1-remainder(round_up_indicies));
        normal_unit_vectors_at_stations(round_up_indicies,:)   = ...
            (fraction_to_midpoint).*normal_unit_vectors_at_midpoints(segment_numbers(round_up_indicies),:) ...
            + (1-fraction_to_midpoint).*normal_unit_vectors_at_joints(segment_numbers(round_up_indicies)+1,:);
        mags = sum(normal_unit_vectors_at_stations(round_up_indicies,:).^2,2).^0.5;
        normal_unit_vectors_at_stations(round_up_indicies,:) = normal_unit_vectors_at_stations(round_up_indicies,:)./mags;
    end
end

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
    
    % Make the axis normally shaped
    axis equal;   

    

    % Plot the query station points
    plot(X_central_at_stations,Y_central_at_stations,'r.','Markersize',30);

    % Show the unit vectors at the stations
    quiver(X_central_at_stations,Y_central_at_stations,normal_unit_vectors_at_stations(:,1),normal_unit_vectors_at_stations(:,2),0,'r','Linewidth',2);

    % Add a legend
    legend('Central traversal','Station query points','Unit vectors','Location','southwest');

    if 1==0
        % Plot the modpoint tangent and normal vectors
        midpoints = [(X_central(1:end-1,1)+X_central(2:end,1))/2, (Y_central(1:end-1,1)+Y_central(2:end,1))/2 ];
        plot(midpoints(:,1),midpoints(:,2),'b.','MarkerSize',20);
        quiver(midpoints(:,1),midpoints(:,2),tangent_unit_vectors_at_midpoints(:,1),tangent_unit_vectors_at_midpoints(:,2),0,'b','Linewidth',3);
        quiver(midpoints(:,1),midpoints(:,2),normal_unit_vectors_at_midpoints(:,1),normal_unit_vectors_at_midpoints(:,2),0,'b','Linewidth',3);

        % Plot the joint tangent vectors
        plot(X_central(:,1),Y_central(:,1),'k.','MarkerSize',20);
        quiver(X_central(:,1),Y_central(:,1),normal_unit_vectors_at_joints(:,1),normal_unit_vectors_at_joints(:,2),0,'g','Linewidth',3);
    end
end % Ends the flag_do_plots if statement

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

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

