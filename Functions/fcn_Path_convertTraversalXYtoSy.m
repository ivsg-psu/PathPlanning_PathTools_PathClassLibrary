function [closestXs, closestYs, closestDistances] = ...
    fcn_Path_convertTraversalXYtoSy(...
    reference_station_points, reference_traversal, ...
    all_traversals, varargin)

% fcn_Path_convertTraversalXYtoSy
% Given a reference traversal and a set of stations along that traversal,
% finds the location on each nearby traversal that is closest to the central
% traveral at each station point. Closest is defined via an orthogonal
% projection (or modifications of orthogonal projections) from the central
% traversal outward toward nearby traversals. The results are then
% presented as an array of lateral offsets.
%
% The process to calculate the results is to project orthogonally from the
% reference traversal outward, depending on the flag type, to find the
% intersection from the reference traversal to other traversals. These
% intersections are then used to calculate lateral offsets.
%
% NOTE: this is not a function to convert from XY to Station coordinates,
% as it does not work for all possible XY positions. It will only work for
% very specific orthogonal projections.
%
% NOTE: this is just a reformulation of the function:
% fcn_Path_FindOrthogonalHitFromTraversalToTraversal
% The main difference is the final plotting result, but otherwise it uses
% the same functionality.
%
% FORMAT:
%
%      [closestXs, closestYs, closestDistances] = ...
%        fcn_Path_convertTraversalXYtoSy(...
%        reference_station_points, reference_traversal, all_traversals,...
%        (flag_rounding_type),(search_radius),(fig_num));
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
% EXAMPLES:
%
% See the script: script_test_fcn_Path_convertTraversalXYtoSy
% for a full test suite.
%
% DEPENDENCIES:
%
%      fcn_Path_FindOrthogonalHitFromTraversalToTraversals
%      fcn_Path_findOrthogonalTraversalVectorsAtStations
%      fcn_DebugTools_checkInputsToFunctions
%
% This function was written on 2021_03_21 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2021_03_21:
% -- first write of the code moving this functionality out of
% fcn_Path_findOrthoScatterFromTraversalToTraversals.m
% 2022_01_07:
% -- edited header to rename:
% fcn_Path_findOrthoScatterFromTraversalToTraversals to
% fcn_Path_FindOrthogonalHitFromTraversalToTraversal
% 2023_04_30:
% -- improved plotting by showing XY and Sy side-by-side
% -- edited header to make it more clear that this is NOT XY to Sy
% converter to use.
% 2025_06_23 - S. Brennan
% -- Updated debugging and input checks

% TO-DO
% (none)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 6; % The largest Number of argument inputs to the function
flag_max_speed = 0;
if (nargin==MAX_NARGIN && isequal(varargin{end},-1))
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
        narginchk(3,MAX_NARGIN);

        % Check the station input
        fcn_DebugTools_checkInputsToFunctions(reference_station_points, 'station');

        % Check the reference_traversal input
        fcn_DebugTools_checkInputsToFunctions(reference_traversal, 'traversal');

        % Check the all_traversals input
        fcn_DebugTools_checkInputsToFunctions(all_traversals, 'traversals');


    end
end

% Fil in default values
Nstations = length(reference_station_points(:,1));
Ntraversals = length(all_traversals.traversal);
Station_central = reference_traversal.Station;

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
flag_do_plots = 0; % Default is to NOT show plots
if (0==flag_max_speed) && (MAX_NARGIN == nargin)
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        fig_num = temp;
        figure(fig_num);
        flag_do_plots = 1;
    end
end

if flag_do_debug
    fig_debug = 35455; %#ok<NASGU>
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

%% Initialize matrices
% These are used to save the current reference traversal (so we can
% check changes after this is updated) and initialize arrays for the loop

closestXs = zeros(Nstations,Ntraversals);
closestYs = zeros(Nstations,Ntraversals);
closestDistances = zeros(Nstations,Ntraversals);


%% For each traversal, project from reference orthogonally
% Search nearest points from reference path to each of the other
% traversals, saving the results as "closest" values

for ith_traversal = 1:Ntraversals
    nearby_traversal = all_traversals.traversal{ith_traversal};

    [closest_path_points,closest_distances] = ...
        fcn_Path_findOrthogonalHitFromTraversalToTraversal(...
        reference_station_points,...
        reference_traversal,...
        nearby_traversal,...
        flag_rounding_type,...
        search_radius, -1);

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

    figure(fig_num);
    clf;

    subplot(1,2,1);
    hold on;
    grid on;

    xlabel('X [m]');
    ylabel('Y [m]');
    axis equal;

    % Plot the central traversal
    plot(reference_traversal.X,reference_traversal.Y,'k.-','Linewidth',3,'Markersize',25);

    % Plot the paths
    h_plots = fcn_Path_plotTraversalsXY(all_traversals,fig_num);
    colors = zeros(length(h_plots),3);
    for ith_plot = 1:length(h_plots)
        set(h_plots(ith_plot),'LineWidth',3)
        colors(ith_plot,:) = get(h_plots(ith_plot),'Color');
    end

    % Plot the station points
    % Find the unit normal vectors at each of the station points
    [unit_normal_vector_start, unit_normal_vector_end] = ...
        fcn_Path_findOrthogonalTraversalVectorsAtStations(...
        reference_station_points,reference_traversal,flag_rounding_type,-1);

    plot(unit_normal_vector_start(:,1),unit_normal_vector_start(:,2),'r.','Markersize',25);

    % Show the unit vectors
    normal_unit_vectors_at_stations = ...
        unit_normal_vector_end - unit_normal_vector_start;
    quiver(unit_normal_vector_start(:,1),unit_normal_vector_start(:,2),...
        normal_unit_vectors_at_stations(:,1),normal_unit_vectors_at_stations(:,2),0,'g','Linewidth',3);  % The 0 term is to prevent scaling

    % Plot the hit results
    for ith_traversal = 1:length(closestXs(1,:))
        plot(closestXs(:,ith_traversal),closestYs(:,ith_traversal),'o','Markersize',5,'Linewidth',3,'Color',colors(ith_traversal,:));
        for ith_point = 1:length(closestXs(:,1))
            text(closestXs(ith_point,ith_traversal),closestYs(ith_point,ith_traversal),sprintf('%.0d',ith_point),'FontSize',12,'VerticalAlignment','bottom','Color',colors(ith_traversal,:));
        end
    end


    subplot(1,2,2);
    hold on; grid on;
    xlabel('Station [m]');
    ylabel('y-offset [m]');
    % axis equal;

    % Plot the results
    for ith_traversal = 1:length(closestDistances(1,:))
        plot(reference_station_points, closestDistances(:, ith_traversal),'.-','Markersize',15,'Linewidth',3,'Color',colors(ith_traversal,:));
        for ith_point = 1:length(closestDistances)
            text(reference_station_points(ith_point,1),closestDistances(ith_point,ith_traversal),sprintf('%.0d',ith_point),...
                'FontSize',12,'VerticalAlignment','bottom','Color',colors(ith_traversal,:));
        end
    end

    % For debugging
    if 1== flag_do_debug  % The following plots the XY results
        figure(fig_num+1);
        clf;
        hold on;
        grid on;

        xlabel('X [m]');
        ylabel('Y [m]');
        axis equal;

        % Plot the central traversal
        plot(reference_traversal.X,reference_traversal.Y,'k.-','Linewidth',3,'Markersize',25);

        % Plot the paths
        fcn_Path_plotTraversalsXY(all_traversals,fig_num);

        % Plot the station points
        % Find the unit normal vectors at each of the station points
        [unit_normal_vector_start, unit_normal_vector_end] = ...
            fcn_Path_findOrthogonalTraversalVectorsAtStations(...
            reference_station_points,reference_traversal,flag_rounding_type,-1);

        plot(unit_normal_vector_start(:,1),unit_normal_vector_start(:,2),'r.','Markersize',35);


        % Plot the results
        for i_point = 1:length(closestXs(:,1))
            plot(closestXs(i_point,:),closestYs(i_point,:),'bo-','Markersize',15,'Linewidth',3);
        end

        % Show the unit vectors
        normal_unit_vectors_at_stations = ...
            unit_normal_vector_end - unit_normal_vector_start;
        quiver(unit_normal_vector_start(:,1),unit_normal_vector_start(:,2),...
            normal_unit_vectors_at_stations(:,1),normal_unit_vectors_at_stations(:,2),0,'g','Linewidth',3);  % The 0 term is to prevent scaling
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

