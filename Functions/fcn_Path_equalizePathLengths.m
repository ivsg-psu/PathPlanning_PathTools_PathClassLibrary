function cellArrayOfEqualizedPaths = ...
    fcn_Path_equalizePathLengths(cellArrayOfUnequalPaths,varargin)
% fcn_Path_equalizePathLengths
% given a cell array of paths that are nominally following the same average
% path, extends the starts and ends of arrays of "short" paths such that
% the projection is fully defined for the "longest" path. The longest path
% is the one whose correction lengths - how much must be added to the
% beginning and end - is smallest before corrections are added.
%
% FORMAT: 
%
%      cellArrayOfEqualizedPaths = ...
%      fcn_Path_equalizePathLengths(...
%            cellArrayOfUnequalPaths,
%            (fig_num));
%
% INPUTS:
%
%     cellArrayOfUnequalPaths: a cell array of paths, with each path
%     containing Nx2 vectors containing [X, Y] positions. Note, the value
%     of N may be different for each entry in the cell array.
%
%     (OPTIONAL INPUTS)
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%     cellArrayOfEqualizedPaths: a cell array of paths, with each path
%     containing Nx2 vectors containing [X, Y] positions. Note, the value
%     of N may be different for each entry in the cell array.
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_Path_findTraversalWithMostData
%      fcn_Path_findOrthogonalTraversalVectorsAtStations
%      fcn_Path_findOrthoScatterFromTraversalToTraversals
%      fcn_Path_cleanPathFromForwardBackwardJogs
%      fcn_Path_plotTraversalsXY
%      fcn_Path_newTraversalByStationResampling
%      fcn_Path_convertPathToTraversalStructure
%
% EXAMPLES:
%      
%     See the script: 
%     script_test_fcn_Path_equalizePathLengths
%     script_demo_fcn_Path_findAverageTraversal
%     for a full test suite and demonstration.
%
% This function was written on 2020_11_15 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%     
% 2025_06_26  - S. Brennan
% -- wrote the code originally 

% TO DO
% (none)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 2; % The largest Number of argument inputs to the function
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

flag_do_debug = 1;

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
        narginchk(1,MAX_NARGIN);

        % Check the cellArrayOfUnequalPaths input
        % fcn_DebugTools_checkInputsToFunctions(cellArrayOfUnequalPaths, 'traversals');

    end
end

% Does user want to show the plots?
flag_do_plots = 0; % Default is to NOT show plots
if (0==flag_max_speed) && (MAX_NARGIN == nargin) 
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        fig_num = temp;
        flag_do_plots = 1;
    end
end

if flag_do_debug
    fig_debug = 78787; %#ok<NASGU>
end



%% Main code starts here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define preliminaries and set up debugging
Npaths = length(cellArrayOfUnequalPaths); % the number of traversals we will be averaging

if flag_do_debug
    figure(fig_debug)
    clf;
    hold on;
    axis equal

    colorsPerData = zeros(Npaths,3);
    for ith_data = 1:Npaths
        h_plot = plot(cellArrayOfUnequalPaths{ith_data}(:,1), cellArrayOfUnequalPaths{ith_data}(:,2),...
            '.-','LineWidth',5,'MarkerSize',20, 'DisplayName',sprintf('Traversal %.0f',ith_data));
        colorsPerData(ith_data,:) = get(h_plot,'Color');
    end
    legend
end  

%%%%%


%% Find all the endPoints and endVectors of all the traversals
allEndPositions   = nan(Npaths,2); % Initialize variables
allEndVectors     = nan(Npaths,2); % Initialize vectors

for ith_path = 1:Npaths
    allEndPositions(ith_path,:) = cellArrayOfUnequalPaths{ith_path}(end,:);
    allEndVectors(ith_path,:)   = cellArrayOfUnequalPaths{ith_path}(end,:) - cellArrayOfUnequalPaths{ith_path}(end-1,:);
end

% Find the longest distance from one path's tail to another
longestDistance = fcn_INTERNAL_findMaxTailDistanceDisparity(allEndPositions);
 
% Using the longest distance, see which traversal extends out the
% furthest from the others. 
index_of_longest = fcn_INTERNAL_findPathThatSticksOutMost(cellArrayOfUnequalPaths, longestDistance);

%% Project paths to be exactly as long as the longest one
% To do this, we calculate the vector projection of each end point that
% would push the end point to align with the longest
longest_pathPoint = allEndPositions(index_of_longest,:);
longest_pathVector = allEndVectors(index_of_longest,:);
orthoToEnd = longest_pathVector*[0 1; -1 0];

% Intersect each path's ending vector with the orthogonal projection of the
% longest point. We do this with the sensor hit method:
% FORMAT:
%      [distance, location, wall_segment, t, u] = ...
%         fcn_Path_findSensorHitOnWall(...
%         wall_start, wall_end,...
%         sensor_vector_start,sensor_vector_end,...
%         (flag_search_return_type), (flag_search_range_type), ...
%         (tolerance), (fig_num))

[~, newEndLocations, wall_segment, ~, ~] = ...
    fcn_Path_findSensorHitOnWall(...
    allEndPositions-allEndVectors, allEndPositions,...
    longest_pathPoint,longest_pathPoint+orthoToEnd,...
    (1), (3), ...
    ([]), (234));

cellArrayOfEqualizedPaths = cell(Npaths,1);
for ith_wall = 1:Npaths
    this_wall = wall_segment(ith_wall);
    if ith_wall~=index_of_longest
        % Make path longer
        cellArrayOfEqualizedPaths{this_wall} = [cellArrayOfUnequalPaths{this_wall}; newEndLocations(this_wall,:)];
    else
        % Keep original path
        cellArrayOfEqualizedPaths{this_wall} = cellArrayOfUnequalPaths{this_wall};
    end
end

if flag_do_debug
    for ith_data = 1:Npaths
        plot(cellArrayOfEqualizedPaths{ith_data}(:,1), cellArrayOfEqualizedPaths{ith_data}(:,2),'.-',...
            'LineWidth',2,'MarkerSize',10,'Color',colorsPerData(ith_data,:)*0.8, 'DisplayName',sprintf('Resampled %.0f',ith_data));
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
if flag_do_plots
    % plot the final XY result

    % Prep the figure for plotting
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end
    
    % Is this 2D or 3D?
    dimension_of_points = 2;

    % Find size of plotting domain
    data
    allPointsBeingPlotted = [(1:length(diff_angles))' diff_angles*180/pi];
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
        if 1==1
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

    fcn_Path_plotTraversalsXY(data,fig_num);
    hold on;
    plot(traversal_average.X,traversal_average.Y,'b.-','Linewidth',4,'Markersize',20,'DisplayName','Average');
    
    legend;
  
    
    if flag_do_debug
        % Plot the path convergence
        figure(2255);
        clf;
        hold on;
        xlabel('Index')
        ylabel('Position change between iterations [m]')
        for ith_iteration = 1:length(iteration_error_X)
            title(sprintf('Iteration %.0d of %.0d',ith_iteration,length(iteration_error_X)));
            plot(average_traversals{ith_iteration}.X,average_traversals{ith_iteration}.Y,'-'); 
            drawnow;
            pause(0.1);
        end
        
        % Plot the error convergence
        figure(2233);
        clf;
        hold on;
        xlabel('Index')
        ylabel('Position change between iterations [m]')
        for ith_iteration = 1:length(iteration_error_X)
            title(sprintf('Iteration %.0d of %.0d',ith_iteration,length(iteration_error_X)));
            plot(total_error{ith_iteration});            
            pause(0.1);
        end
        
        
        % Plot the error convergence
        figure(2244);
        clf;
        
        semilogy(1:length(mean_error),mean_error,'.-','Markersize',30);
        xlabel('Iteration')
        ylabel('Mean change in distance')
        grid on
        
    end
    
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end
end % End of function


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

%% fcn_INTERNAL_findMaxTailDistanceDisparity
function longestDistance = fcn_INTERNAL_findMaxTailDistanceDisparity(allEndPositions)
% Finds longest distance from one path's end to another
% note: do the square-root calculation only once at end, rather than at
% every iteration, since this is slow.

Ntraversals = length(allEndPositions(:,1)); % the number of paths we will be checking

longestDistance = -inf;
for ith_traversal = 1:Ntraversals
    for jth_traversal = 1:Ntraversals
        if jth_traversal~=ith_traversal
            testDistance = sum((allEndPositions(jth_traversal,:)-allEndPositions(ith_traversal,:)).^2,2);
            if testDistance>longestDistance
                longestDistance = testDistance;
            end
        end
    end
end
longestDistance = longestDistance^0.5;

end % Ends fcn_INTERNAL_findMaxTailDistanceDisparity

%% fcn_INTERNAL_findPathThatSticksOutMost
function index_of_longest = fcn_INTERNAL_findPathThatSticksOutMost(cellArrayOfUnequalPaths, longestDistance)

Npaths = length(cellArrayOfUnequalPaths); % the number of paths

% Using the longest distance, see which traversal extends out the
% furthest from the others. This can done by projecting, orthogonally,
% from each traversal's end point and checking if anything was hit,
% using a projection of the longest distance. NOTE: this may not work
% if the trajectories loop back onto each other.

pathIndicesNoHits = zeros(Npaths,1);

for ith_testPath = 1:Npaths
    current_testPath = cellArrayOfUnequalPaths{ith_testPath};    
    current_testTraversal = fcn_Path_convertPathToTraversalStructure(current_testPath, -1);
    endStation = current_testTraversal.Station(end);

    all_orthogonalProjectionHits = nan(Npaths,1);
    for jth_adjacentPath = 1:Npaths
        if jth_adjacentPath~=ith_testPath

            nearby_path = cellArrayOfUnequalPaths{jth_adjacentPath};            
            nearby_traversal = fcn_Path_convertPathToTraversalStructure(nearby_path, -1);

            % Calculate the closest point and distance on the nearby path
            % FORMAT:
            % [closest_path_point,distances] = ...
            %     fcn_Path_findOrthogonalHitFromTraversalToTraversal(stations,...
            %     central_traversal,nearby_traversal,...
            %     flag_rounding_type,search_radius,fig_num);

            [~, distance] = ...
                fcn_Path_findOrthogonalHitFromTraversalToTraversal(endStation,...
                current_testTraversal,nearby_traversal,...
                [], longestDistance, -1);
            all_orthogonalProjectionHits(jth_adjacentPath,1) = distance;
        end
    end
    if all(isnan(all_orthogonalProjectionHits))
        pathIndicesNoHits(ith_testPath) = 1;
    end
end
index_of_longest = find(pathIndicesNoHits);
if isempty(index_of_longest)
    warning('on','backtrace');
    warning('No path found that does not intersect the others. This should not happen');
    error('Unable to continue');
elseif length(index_of_longest)>1
    warning('on','backtrace');
    warning('Multiple paths found that stick out beyond all others.');
    warning('This usually only occurs when paths are severely mis-aligned.');
    warning('This suggests a poor set of paths used for path averaging.');
    warning('The first "stick out" path found will be used, path number: %.0f ',index_of_longest);   
    warning('but caution should be used in accepting the results.');
end
% Force only one index to be kept, if there are more than one
index_of_longest = index_of_longest(1);
end % Ends fcn_INTERNAL_findPathThatSticksOutMost
