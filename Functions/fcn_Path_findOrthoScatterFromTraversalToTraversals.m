function [closestXs, closestYs, closestDistances] = ...
    fcn_Path_findOrthoScatterFromTraversalToTraversals(...
    reference_station_points, reference_traversal, ...
    all_traversals, varargin)

% fcn_Path_findOrthoScatterFromTraversalToTraversals
% Given a central traversal and a set of stations along that traversal,
% finds the locations on the nearby traversals that are closest to the central
% traveral at each station point. Closest is defined via an orthogonal
% projection (or modifications of orthogonal projections) from the central
% traversal outward toward nearby traversals. Note that this function is
% essentially the multi-traversal version of:
% fcn_Path_findOrthogonalHitFromTraversalToTraversal, which is called
% internally over/over for each traversal.
%
% FORMAT:
%
%      [closestXs, closestYs, closestDistances] = ...
%        fcn_Path_findOrthoScatterFromTraversalToTraversals(...
%        reference_station_points, reference_traversal, all_traversals,...
%        (flag_rounding_type), (search_radius), (fig_num));
%
% INPUTS:
%
%      reference_station_points: an N x 1 vector containing the stations on
%      the reference_traversal where the orthogonal projections should take
%      place
%
%      reference_traversal: a traversal structure that specifies the
%      traversal where projections to other traversals are taking place.
%
%      all_traversals: a traversals type data structure containing a cell
%      array of traversal structures that specifies the traversals where
%      intersections from the reference_traversal would hit.
%
%      (OPTIONAL INPUTS)
%      flag_rounding_type: a flag to indicate which type of projection is
%      used, especially when stations are located at the end-points of
%      segments within the nearby_traversal. Note that the very first point
%      always uses projections from the following segement, and the very
%      last point always uses the prior. The flag determines behaviors for
%      endpoints of internal segments. The options include:
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
%          vectors. Note that this can generate non-orthogonal results.
%
%      search_radius: how far in station distance for the search points to
%      look (default is 3 times the differences in station lengths within
%      all_traversals )
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%      closestXs:  a NxM vector containing the [X] location of
%      the nearest points at the N stations projected orthogonally to the M
%      trajectories with all_traversals
%
%      closestYs:  a NxM vector containing the [Y] location of
%      the nearest points at the N stations projected orthogonally to the M
%      trajectories with all_traversals
%
%      closestDistancess:  a NxM vector containing the distance of
%      the nearest points at the N stations projected orthogonally to the M
%      trajectories with all_traversals. Note that positive distances are
%      those whose cross product from the reference_trajectory to the
%      intersection is positive, negative distances are in the opposite
%      direction
%
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_Path_findOrthogonalHitFromTraversalToTraversal
%      fcn_Path_findOrthogonalTraversalVectorsAtStations
%      fcn_Path_plotTraversalsXY
%
% EXAMPLES:
%
% See the script: script_test_fcn_Path_findOrthoScatterFromTraversalToTraversals
% for a full test suite.
%
% This function was written on 2021_01_02 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2021_01_02:
% -- first write of the code moving this functionality out of
% fcn_Path_findAveragePathViaOrthogonalProjection.m
% 2021_01_07:
% -- renamed function to clarify paths versus traversals
% 2021_01_09:
% -- corrected terminology in comments
% 2021_12_27:
% -- corrected dependencies in comments
% -- fixed name fcn_Path_findOrthogonalHitFromTraversalToTraversal
% 2022_01_03
% -- found a bug in the constrainted search functionality,
% fcn_Path_findOrthogonalHitFromTraversalToTraversal, fixed it!
% -- in the case of single traversal queries, added printing of the
% distance to the hit
% 2025_06_23 - S. Brennan
% -- Updated debugging and input checks

% TO-DO
% (none)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==6 && isequal(varargin{end},-1))
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
        narginchk(3,6);

        % Check the station input
        fcn_DebugTools_checkInputsToFunctions(reference_station_points, 'station');

        % Check the reference_traversal input
        fcn_DebugTools_checkInputsToFunctions(reference_traversal, 'traversal');

        % Check the all_traversals input
        fcn_DebugTools_checkInputsToFunctions(all_traversals, 'traversals');
    end
end

% Does the user want to give a rounding type?
flag_rounding_type = 1;
if 4 <= nargin
    flag_rounding_type = varargin{1};
end

% Define search radius?
search_radius = max(reference_traversal.Station)*3;
if 5 <= nargin
    search_radius = varargin{2};
end

% Does user want to show the plots?
flag_do_plots = 0; % Default is to NOT show plots
if (0==flag_max_speed) && (6 == nargin) 
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        fig_num = temp;
        figure(fig_num);
        flag_do_plots = 1;
    end
else
    if flag_do_debug
        fig_debug = 3454; 
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

% Fil in default values
Nstations = length(reference_station_points(:,1));
Ntraversals = length(all_traversals.traversal);

%% Initialize matrices
% These are used to save the current reference traversal (so we can
% check changes after this is updated) and initialize arrays for the loop

closestXs = zeros(Nstations,Ntraversals);
closestYs = zeros(Nstations,Ntraversals);
closestDistances = zeros(Nstations,Ntraversals);

if flag_do_debug
    % Plot the paths
    fcn_Path_plotTraversalsXY(all_traversals,27272);
end

%% For each traversal, project from reference orthogonally
% Search nearest points from reference path to each of the other
% traversals, saving the results as "closest" values

for ith_traversal = 1:Ntraversals
    nearby_traversal = all_traversals.traversal{ith_traversal};
    
    if flag_do_debug
        [closest_path_points,closest_distances] = ...
            fcn_Path_findOrthogonalHitFromTraversalToTraversal(...
            reference_station_points,...
            reference_traversal,...
            nearby_traversal,...
            flag_rounding_type,...
            search_radius,fig_debug);
    else
        [closest_path_points,closest_distances] = ...
            fcn_Path_findOrthogonalHitFromTraversalToTraversal(...
            reference_station_points,...
            reference_traversal,...
            nearby_traversal,...
            flag_rounding_type,...
            search_radius, -1);
    end
    
    % Save final results as closest points
    closestXs(:,ith_traversal) = closest_path_points(:,1);
    closestYs(:,ith_traversal)  = closest_path_points(:,2);
    closestDistances(:,ith_traversal) = closest_distances;
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
temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end
    
    % Is this 2D or 3D?
    dimension_of_points = 2;
    
    % Find size of plotting domain
    allPointsBeingPlotted = [reference_traversal.X reference_traversal.Y; closestXs(:,1),closestYs(:,1)];
    
    max_plotValues = max(allPointsBeingPlotted);
    min_plotValues = min(allPointsBeingPlotted);
    sizePlot = max(max_plotValues) - min(min_plotValues);
    nudge = sizePlot*0.006; %#ok<NASGU>

    % Find size of plotting domain
    if flag_rescale_axis
        percent_larger = 0.3;
        axis_range = max_plotValues - min_plotValues;
        if (0==axis_range(1,1))
            axis_range(1,1) = 2/percent_larger;
        end
        if (0==axis_range(1,2))
            axis_range(1,2) = 2/percent_larger;
        end
        if dimension_of_points==3 && (0==axis_range(1,3))
            axis_range(1,3) = 2/percent_larger;
        end

        % Force the axis to be equal?
        if 1==0
            min_valuesInPlot = min(min_plotValues);
            max_valuesInPlot = max(max_plotValues);
        else
            min_valuesInPlot = min_plotValues;
            max_valuesInPlot = max_plotValues;
        end

        % Stretch the axes
        stretched_min_vertexValues = min_valuesInPlot - percent_larger.*axis_range;
        stretched_max_vertexValues = max_valuesInPlot + percent_larger.*axis_range;
        axesTogether = [stretched_min_vertexValues; stretched_max_vertexValues];
        newAxis = reshape(axesTogether, 1, []);
        axis(newAxis);

    end
    % goodAxis = axis;

    hold on;
    grid on;

    xlabel('X [m]');
    ylabel('Y [m]');
    axis equal;
    
    % Plot the central traversal
    plot(reference_traversal.X,reference_traversal.Y,'k.-','Linewidth',3,'Markersize',25);
    
    % Plot the paths
    fcn_Path_plotTraversalsXY(all_traversals,fig_num);
        
    % Setup to plot the station points (green dots) and sensor vectors
    % (green and cyan arrows)
    
    % Find the unit normal vectors at each of the station points
    [unit_normal_vector_start, unit_normal_vector_end] = ...
        fcn_Path_findOrthogonalTraversalVectorsAtStations(...
        reference_station_points,reference_traversal,flag_rounding_type, -1);

    unit_vector_displacement = unit_normal_vector_end - unit_normal_vector_start;
    sensor_vector_start = unit_normal_vector_start;
    postive_sensor_vector_end = unit_normal_vector_start + unit_vector_displacement*search_radius;
    negative_sensor_vector_end = unit_normal_vector_start - unit_vector_displacement*search_radius;
    
    % Plot the sensor vector origin as green dots
    plot(sensor_vector_start(:,1),sensor_vector_start(:,2),'g.','Markersize',35);
       
    % Show the sensor vectors as green arrows (+) and cyan (-)
    positive_normal_vectors_at_stations = ...
        postive_sensor_vector_end - sensor_vector_start;
    negative_normal_vectors_at_stations = ...
        negative_sensor_vector_end - sensor_vector_start;
    quiver(sensor_vector_start(:,1),sensor_vector_start(:,2),...
        positive_normal_vectors_at_stations(:,1),positive_normal_vectors_at_stations(:,2),0,'g','Linewidth',3);  % The 0 term is to prevent scaling
    quiver(sensor_vector_start(:,1),sensor_vector_start(:,2),...
        negative_normal_vectors_at_stations(:,1),negative_normal_vectors_at_stations(:,2),0,'c','Linewidth',3);  % The 0 term is to prevent scaling

    % Plot the hits as red circles
    INTERNAL_plot_only_hits(closestXs,closestYs);
        
    % If there is only one query trajectory, print the distance
    if isscalar(all_traversals.traversal)
        
        % Add a legend - note that the line style, etc tries to be consistent
        % with the legend in:
        % fcn_Path_findOrthogonalHitFromTraversalToTraversal
        legend('Central traversal','Path to check','Station query points','Sensor vectors (+)','Sensor vectors (-)','Hit locations');
        
        % Add text to indicate distance result
        text_locations = sensor_vector_start + unit_vector_displacement.*closestDistances/2;
        for ith_distance = 1:length(closestDistances)
            text(text_locations(ith_distance,1),text_locations(ith_distance,2),sprintf('%.2f',closestDistances(ith_distance,1)),'Color',[1 0 0]);
        end
    end

end % Ends the flag_do_plot if statement

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

function INTERNAL_plot_only_hits(all_xdata,all_ydata)

xdata_to_keep = [];
ydata_to_keep = [];
for i_point = 1:length(all_xdata(:,1))
    
    % Grab the ith row of data for x and y, convert into columns
    xdata = all_xdata(i_point,:)';
    ydata = all_ydata(i_point,:)';
    
    % Find indices that are not NaN
    good_x_indices = find(~isnan(xdata));
    good_y_indices = find(~isnan(ydata));
    good_indices = intersect(good_x_indices,good_y_indices);
    
    % Save values
    xdata_to_keep = [xdata_to_keep; NaN; xdata(good_indices)]; %#ok<AGROW>
    ydata_to_keep = [ydata_to_keep; NaN; ydata(good_indices)]; %#ok<AGROW>
end

% Put all the data into a plot
plot(xdata_to_keep,ydata_to_keep,'ro-','Markersize',15,'Linewidth',3);

end
