function [closestXs,closestYs,closestZs,closestYaws] = ...
    fcn_Path_findClosestPointsToTraversal(...
    reference_traversal,data,varargin)
% This function finds the projection point on each traversal from a
% reference path by snapping each vertex of the reference traversal onto
% the nearby traversals. Each "snap" calculates the nearest point within
% each other traversal by measuring the orthogonal distance from the nearby
% traversal, to the point on the reference_traversal. Thus, the "snap
% points" are always orthogonal to each nearby traversal, and represent the
% closest that traversal gets to the reference_traversal when measured at
% the reference_traversal points.
%
% FORMAT:
%
%      [closestXs,closestYs,closestZs,closestYaws] = ...
%      fcn_Path_findClosestPointsToTraversal(path, data,...
%        (flag_yaw), (flag_3D), (fig_num))
%
% INPUTS:
%
%      reference_traversal: a traversal structure that specifies the
%      traversal where projections to other paths are taking place.
%
%      data: a traversals type data structure, namely a structure
%      containing a cell array of traversals, each with subfields of X, Y,
%      etc. in the following form
%           data.traversal{i_path}.X
%      Note that i_path denotes an index into a different traversal.
%      This array of traversals specifies the traversals where
%      intersections from the reference_traversal would hit. Each traversal
%      will be compared separately.
%
%      (OPTIONAL INPUTS)
%
%      (flag_yaw)     = set to 1 to interplate the yaw
%
%      (flag_3D)      = set to 1 to process 3D data (e.g. X, Y, and Z)
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%      [closestXs,closestYs,closestZs,closestYaws]  = ...
%          the lane type of each lane, format:table
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_Path_snapPointToPathViaVectors
%
% EXAMPLES:
%
% See the script: script_test_fcn_Path_findClosestPointsToTraversal
% for a full test suite.
%
% Author:             Liming Gao and S. Brennan (sbrennan@psu.edu)
% Created Date:       2020-07-21

% Revision History
% 2020_02_22 - S. Brennan
% -- take yaw into account
% 2020_11_12 - S. Brennan
% -- (SNB) added more comments
% 2021_01_07 - S. Brennan
% -- Added more comments
% -- Updated name to reflect change from path to traverals
% 2021_01_09 - S. Brennan
% -- corrected terminology in comments
% -- updated dependencies
% -- fixed snap function name (it was wrong)
% -- added input checking
% 2025_06_23 - S. Brennan
% -- Updated debugging and input checks

% TO-DO
% (none)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 5; % The largest Number of argument inputs to the function
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
        narginchk(2,MAX_NARGIN);

        % Check the reference_traversal input
        fcn_DebugTools_checkInputsToFunctions(reference_traversal, 'traversal');

        % Check the data input
        fcn_DebugTools_checkInputsToFunctions(data, 'traversals');
    end
end

% Does user want to specify flag_yaw?
flag_yaw = 0;
if nargin >= 3
    flag_yaw = varargin{1};
end

% Does user want to specify flag_3D?
flag_3D = 0; % whether we consider the Z when do the projection
if nargin >= 4
    flag_3D = varargin{2};
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
else
    if flag_do_debug
        fig_debug = 4848; %#ok<NASGU>
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate the bound we look for similar station. The bound should be no
% more than the difference between longest path and shortest paths. To know
% this difference, we have to find the last station point for each
% traversal.
Ntraversals = length(data.traversal);

last_station = zeros(Ntraversals,1); % extract the max station of each path
for i_traversal = 1:Ntraversals
    last_station(i_traversal)  = data.traversal{i_traversal}.Station(end);
end

% Calculate the maximum and minimum lengths of stations
longest_station = max(max(last_station), reference_traversal.Station(end));
shortest_station = min(min(last_station), reference_traversal.Station(end));

% Add a safety factor to the differences in s-coordinate distances. This
% difference in distance is the largest s-distance we would ever expect
% between one trajectory and another.

safety_factor = 2;
if Ntraversals>1
    station_bound = safety_factor*(longest_station - shortest_station);
else
    station_bound = safety_factor*longest_station;
end


% initializing variables that will record the closest values of each path
% to the reference path
closestXs = zeros(length(reference_traversal.X),Ntraversals);
closestYs = zeros(length(reference_traversal.X),Ntraversals);
closestZs = zeros(length(reference_traversal.X),Ntraversals);
closestYaws = zeros(length(reference_traversal.X),Ntraversals);

% Loop through each point of the reference traversal
for index_reference_path = 1:length(reference_traversal.X) % loop through all the point on the reference path

    % calculating x, y, and z position of one specific point on
    % each traversal given by index_to_match (a common point). This
    % reference point is what will be used to "match" all the other
    % trajectories, e.g. to find the closest value.
    X_ref = reference_traversal.X(index_reference_path);
    Y_ref = reference_traversal.Y(index_reference_path);
    % Z_ref = reference_traversal.Z(index_reference_path);
    Station_ref = reference_traversal.Station(index_reference_path);

    % Define the minimum and maximum stations that we can look for in each
    % of the trajectories. It will be our reference station that we're at,
    % plus and minus the station bound limits. Keep the limits
    station_lower = max(Station_ref - station_bound, 0);

    station_upper = min(Station_ref + station_bound, longest_station);
    % station_upper = min(Station_ref + station_bound, max(last_station));
    % % 2020_11_12 - above line seems to be wrong so corrected it!

    % Loop through all traversals and for each
    % traversal, find the closest point on the first traversal, and
    % use this to calculate station offset and index offset of the
    % ith traversal relative to the first traversal
    for i_traversal = 1:Ntraversals % loop through all traverasals (already excluded path)

        % allowable us to find the indexes within the station range
        index_in_station_range = find(data.traversal{i_traversal}.Station >= station_lower & data.traversal{i_traversal}.Station <= station_upper);

        % pull the X,Y,Z data that are within this station range
        X_tra = data.traversal{i_traversal}.X(index_in_station_range);
        Y_tra = data.traversal{i_traversal}.Y(index_in_station_range);
        Z_tra = data.traversal{i_traversal}.Z(index_in_station_range);

        % Yaw index cannot be larger than N-1
        length_Yaw = length(data.traversal{i_traversal}.Yaw);
        yaw_index_in_station_range = min(index_in_station_range,length_Yaw);
        Yaw_tra = [data.traversal{i_traversal}.Yaw(yaw_index_in_station_range)];

        path = [X_tra, Y_tra];

        [closest_path_point,~,first_path_point_index,second_path_point_index,percent_along_length]...
            = fcn_Path_snapPointToPathViaVectors([X_ref,Y_ref], path, [], -1);

        closestXs(index_reference_path,i_traversal) = closest_path_point(1);
        closestYs(index_reference_path,i_traversal) = closest_path_point(2);

        % interplate Zup
        nearest_Z1 = Z_tra(first_path_point_index);
        nearest_Z2 = Z_tra(second_path_point_index);
        interp_Z = nearest_Z1 + (nearest_Z2-nearest_Z1)*percent_along_length;
        closestZs(index_reference_path,i_traversal) = interp_Z;

        if flag_yaw == 1
            % find the closet yaw
            nearest_yaw1 = Yaw_tra(first_path_point_index);
            nearest_yaw2 = Yaw_tra(second_path_point_index);
            Proj_yaw_ref_to_tra = nearest_yaw1 + (nearest_yaw2-nearest_yaw1)*percent_along_length;
            closestYaws(index_reference_path,i_traversal) = Proj_yaw_ref_to_tra;
        end


        if flag_3D
            error('3D methods are not yet implemented');

            %             % finding the idx of the that is closet to the current reference point
            %             id = knnsearch([X_tra Y_tra Z_tra],[X_ref Y_ref Z_ref],'k',2);  %Find k-nearest neighbors using object
            %             id = sort(id); % sort in Ascending Order
            %
            %             % find the projection point on traversal path
            %             nearest_point1 = [X_tra(id(1)); Y_tra(id(1));Z_tra(id(1))];
            %             nearest_point2 = [X_tra(id(2)); Y_tra(id(2));Z_tra(id(2))];
            %             query_point = [X_ref; Y_ref; Z_ref];
            %         else
            %             % finding the idx of the other traversal that is closet to the current reference point
            %             id = knnsearch([X_tra Y_tra],[X_ref Y_ref],'k',2);  %Find k-nearest neighbors using object
            %             id = sort(id); % sort in Ascending Order
            %
            %             % find the projection point on traversal path
            %             nearest_point1 = [X_tra(id(1)); Y_tra(id(1))];
            %             nearest_point2 = [X_tra(id(2)); Y_tra(id(2))];
            %             query_point = [X_ref; Y_ref];
            %         end
            %
            %         % s = S(id(1));  % The station of the nearest point on the map
            %         % The vector from nearest point1 to nearest point2 on the traversal path
            %         np1_to_np2 = nearest_point2 - nearest_point1;
            %
            %         % The vector from nearest point1 to query point
            %         np1_to_qp = query_point - nearest_point1;
            %
            %         % the length vector
            %         np1_to_np2_Length = sum(np1_to_np2.^2).^0.5;
            %
            %         %pp means proejction point,and np1_to_pp_length is the projection
            %         %length of np1_to_qp to np1_to_np2, i.e. the distance from nearest
            %         %point1 to projection point (the value can be negtive)
            %         np1_to_pp_length = dot(np1_to_qp,np1_to_np2)./ np1_to_np2_Length;
            %
            %         Proj_point_ref_to_tra = nearest_point1 + np1_to_np2* diag(np1_to_pp_length/np1_to_np2_Length); %projection point  of second station on the reference path (first travel)
            %
            %         if flag_3D == 1
            %             closestXs(index_reference_path,i_traversal) = Proj_point_ref_to_tra(1);
            %             closestYs(index_reference_path,i_traversal) = Proj_point_ref_to_tra(2);
            %             closestZs(index_reference_path,i_traversal) = Proj_point_ref_to_tra(3);
            %         else
            %             closestXs(index_reference_path,i_traversal) = Proj_point_ref_to_tra(1);
            %             closestYs(index_reference_path,i_traversal) = Proj_point_ref_to_tra(2);
            %
            %             % interplate Zup
            %             nearest_Z1 = Z_tra(id(1));
            %             nearest_Z2 = Z_tra(id(2));
            %             interp_Z = nearest_Z1 + (nearest_Z2-nearest_Z1)*(np1_to_pp_length/np1_to_np2_Length);
            %             closestZs(index_reference_path,i_traversal) = interp_Z;
            %         end
            %
            %
            %         if flag_yaw == 1
            %             % find the closet yaw
            %             nearest_yaw1 = Yaw_tra(id(1));
            %             nearest_yaw2 = Yaw_tra(id(2));
            %             Proj_yaw_ref_to_tra = nearest_yaw1 + (nearest_yaw2-nearest_yaw1)*(np1_to_pp_length/np1_to_np2_Length);
            %             closestYaws(index_reference_path,i_traversal) = Proj_yaw_ref_to_tra;
        end % Ends flad 3D check
    end % Ends loop through ith traversal
end % Ends loop through the reference path

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
    % Prep the figure for plotting
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end

    % Is this 2D or 3D?
    dimension_of_points = 2;

    % Find size of plotting domain
    allPointsBeingPlotted = [reference_traversal.X reference_traversal.Y; closestXs(:,1),closestYs(:,1)];
    for i_traversal = 1:Ntraversals
        allPointsBeingPlotted = [allPointsBeingPlotted; data.traversal{i_traversal}.X data.traversal{i_traversal}.Y]; %#ok<AGROW>
    end
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

    % INPUTS
    % Plot the reference_traversal
    plot(reference_traversal.X,reference_traversal.Y,'r','Linewidth',2);

    % Plot the path
    for i_traversal = 1:Ntraversals
        plot(data.traversal{i_traversal}.X,data.traversal{i_traversal}.Y,'bo-','Linewidth',2);
    end

    axis equal;

    % Plot the station points
    plot(reference_traversal.X,reference_traversal.Y,'k.','Markersize',15);

    % OUTPUTS
    % Plot hit locations
    plot(closestXs(:,1),closestYs(:,1),'r.','Markersize',30);

    % Show the unit vectors
    quiver(reference_traversal.X,reference_traversal.Y,closestXs(:,1)-reference_traversal.X,closestYs(:,1)-reference_traversal.Y,0,'Linewidth',3);

    legend('Central traversal','Path to check','Station query points','Hit locations');

    % % Label the points with distances?
    % for i_point = 1:length(path(:,1))
    %     text(path(i_point,2),path(i_point,3),sprintf('%.2f',distances_point_to_path(i_point)));
    % end

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end


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
