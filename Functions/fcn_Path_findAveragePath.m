function [path_average, closestXs, closestYs, closestDistances] = ...
    fcn_Path_findAveragePath(cellArrayOfPaths,varargin)
% fcn_Path_findAveragePath
% finds the average path of several paths. 
% Obtains average by projecting each input path in a cell array to each
% other, and taking these projected distances and projecting each path
% toward these projected distances. The result is that each path "moves"
% toward each other. After each move, any "kinks" in the path are removed,
% and the path is resampled. The user can define the stationInterval for
% resampling, the maximum number of iterations allowed, and the exit
% tolerance where, if the station length difference from the longest to
% shortest "averaged" path is less than exit_tolerance, the function will
% exit.
%
% As additional outputs, for each point in the path, this function also
% finds the intersection point in other paths via orthogonal projection of
% the final path_average, and saves these points into arrays to denote the
% closest X and Y coordinates, and the distances. These X, Y, and distance
% arrays have M columns, one for each path, and N rows, one for each
% station in the final path_average.
%
% FORMAT: 
%     [path_average, closestXs, closestYs, closestDistances] = ...
%         fcn_Path_findAveragePath(cellArrayOfPaths,... 
%                 (stationInterval),...
%                 (max_num_iterations),...
%                 (exit_tolerance),...
%                 (fig_num));
%
% INPUTS:
%
%      cellArrayOfPaths: a cell array of paths to be averaged with each
%      other. Each path is a N x 2 or N x 3 set of coordinates
%      representing the [X Y] or [X Y Z] coordinates, in sequence, of a
%      path. The averaging works best if each path starts and stops in
%      approximately the same area and with similar orientations.
%
%      (OPTIONAL INPUTS)
%
%      stationInterval: the station distance interval of the resulting
%      average path, in meters. Default is: stationInterval = 1; % units are meters 
%
%      max_num_iterations: an integer to specify the maximum number of
%      iterations to perform in the path averaging, where the solution of
%      the first processing step serves as the starting paths for
%      the second iteration, etc. Usually, convergence occurs
%      within 3 to 10 iterations (the default is 40). Note that, within the
%      function, a debug flag can be set to plot and analyze convergence.
%
%      exit_tolerance: the allowable tolerance to make the looping exit
%      before the max_num_iterations is reached. Default is 0.1 meters.
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%      path_average: the resulting path representing the average
%      path. Usually, it is decimated at an even spacing which is currently
%      set to 1 meter.
%
%      closestXs:  a N x M vector containing the [X] location of
%      the nearest points at the N average stations projected orthogonally
%      to the M trajectories within cellArrayOfPaths
%
%      closestYs:  a N x M vector containing the [Y] location of
%      the nearest points at the N average stations projected orthogonally
%      to the M trajectories within cellArrayOfPaths
%
%      closestDistancess:  a N x M vector containing the distance of
%      the nearest points at the N average stations projected orthogonally
%      to the M trajectories within cellArrayOfPaths. Note that positive
%      distances are those whose cross product from the
%      reference_trajectory to the intersection is positive, negative
%      distances are in the opposite direction
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_Path_equalizePathLengths
%      fcn_Path_findOrthoScatterFromTraversalToTraversals
%      fcn_Path_cleanPathFromForwardBackwardJogs
%      fcn_Path_convertPathToTraversalStructure
%
% EXAMPLES:
%      
%     See the script: 
%     script_test_fcn_Path_findAveragePath
%     for a full test suite and demonstration.
%
% This function was written on 2020_11_15 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%     
% 2020_11_15  - S. Brennan
% -- wrote the code originally - lots of bugs
% 2020_12_25  - S. Brennan
% -- added more comments
% 2021_01_01  - S. Brennan
% -- fixed the errors with interpolation
% -- fixed the bug with the end point truncating toward start
% 2021_01_06  - S. Brennan
% -- added functions for input checking
% 2021_01_07  - S. Brennan
% -- deleted unused functions
% -- renamed function to reflect the traversal output, not path
% 2021_01_09:  - S. Brennan
% -- corrected terminology in comments
% -- fixed the input argument notation to be traversals
% 2021_12_27:  - S. Brennan
% -- corrected dependencies in comments
% 2022_01_03:  - S. Brennan
% -- corrected typos in comments
% -- fixed a bug where the Z value is not defined in loop
% 2022_01_06:  - S. Brennan
% -- refactored code, added weighted averaging to prevent iteration
%     bouncing
% 2022_01_10:  - S. Brennan
% -- shut off debugging commentsclose
% 2024_03_14 - S. Brennan
% -- shut off traversal_average.Station calculation as it is giving wrong
% length
% 2025_06_23 - S. Brennan
% -- Updated debugging and input checks
% -- Added use of fcn_Path_equalizePathLengths to fix lengths
% -- Updated the input definition list
% -- Full rewrite of the function
% 2025_07_01 - S. Brennan
% -- Removed traversal input type and replaced with cell array of paths
% -- Renamed function from fcn_Path_findAverageTraversalViaOrthoProjection

% TO DO
% Need to clean up the code - lots of code "lint"

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 5; % The largest Number of argument inputs to the function
flag_max_speed = 0;
if (nargin==MAX_NARGIN && isequal(varargin{end},-1))
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % Flag to plot the results for debugging
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
    debug_fig_num = 999978; 
else
    debug_fig_num = []; 
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

        % Check the cellArrayOfPaths input
        if ~iscell(cellArrayOfPaths)
            error('cellArrayOfPaths input must be a cell type');
        end
        for ith_cell = 1:length(cellArrayOfPaths)
            fcn_DebugTools_checkInputsToFunctions(cellArrayOfPaths{ith_cell}, 'path2or3D');
        end
    end
end

% Check to see if user wants to specify stationInterval?
stationInterval = 1; % units are meters
if nargin >= 2
    temp = varargin{1};
    if ~isempty(temp)
        stationInterval = temp;
        if  1==flag_check_inputs && stationInterval<=0
            error('stationInterval must be greater than zero');
        end
    end
end

% Check to see if the number of iterations was specified?
max_num_iterations = 40;  % the default number of iteration to find the average path 
if nargin >= 3
    temp = varargin{2};
    if ~isempty(temp)
        max_num_iterations = temp+1;
        if 1==flag_check_inputs && max_num_iterations<1
            error('Number of iterations must be at least 1');
        end
    end
end

% Check to see if user wants to specify exit_tolerance?
exit_tolerance = 0.1;
if nargin>=4
    temp = varargin{3};
    if ~isempty(temp)
        exit_tolerance = temp;
        if 1==flag_check_inputs && exit_tolerance <= 0
            error('exit_tolerance must be greater than zero');
        end
    end
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

% Set decimation interval and averaging type
Npaths = length(cellArrayOfPaths); % the number of paths we will be averaging

% Set up debug figure by plotting the inputs
if flag_do_debug
    figure(debug_fig_num)
    clf;
    hold on;
    axis equal

    colorsPerData = zeros(Npaths,3);
    for ith_path = 1:Npaths
        h_plot = plot(cellArrayOfPaths{ith_path}(:,1),cellArrayOfPaths{ith_path}(:,2),'.-','LineWidth',5,'MarkerSize',20, 'DisplayName',sprintf('Path %.0f',ith_path));
        colorsPerData(ith_path,:) = get(h_plot,'Color');
        set(h_plot,'Color',(colorsPerData(ith_path,:)*0.4 + 0.6*[1 1 1]))
    end
    legend
end  

%% Make sure all paths start/stop alongside each other
% As well, grab index of the one that requires least modification

% FORMAT: 
%      [cellArrayOfEqualizedPaths, leastExtensionIndex, bestStartIndex, bestEndIndex] = ...
%      fcn_Path_equalizePathLengths(...
%            cellArrayOfUnequalPaths,...
%            (fig_num));
[cellArrayOfEqualizedPaths, ~, ~, ~] = ...
      fcn_Path_equalizePathLengths(...
            cellArrayOfPaths,...
            (-1));


%% Set up all the paths to have same station interval
[cellArrayOfEquidistantPaths, cellArrayOfEquidistantStations] = fcn_INTERNAL_resamplePaths(cellArrayOfEqualizedPaths, stationInterval); 

if flag_do_debug
    figure(debug_fig_num)
    h_updatedPlots = zeros(Npaths,1);
    for ith_path = 1:Npaths
        h_updatedPlots(ith_path,1) = plot(cellArrayOfEquidistantPaths{ith_path}(:,1), cellArrayOfEquidistantPaths{ith_path}(:,2),'.-',...
            'LineWidth',2,'MarkerSize',10,'Color',colorsPerData(ith_path,:)*0.8, 'DisplayName',sprintf('Resampled %.0f',ith_path));
    end
end    

%% Find locked start/end points
allStarts = nan(Npaths,2);
allEnds   = nan(Npaths,2);
for ith_path = 1:Npaths
    allStarts(ith_path,:) = cellArrayOfEquidistantPaths{ith_path}(1,:);
    allEnds(ith_path,:)   = cellArrayOfEquidistantPaths{ith_path}(end,:);
end

lockedStart = mean(allStarts,1);
lockedEnd = mean(allEnds,1);

if flag_do_debug
    plot(lockedStart(1,1), lockedStart(1,2),'.',...
        'LineWidth',2,'MarkerSize',30,'Color',[0 1 0], 'DisplayName',sprintf('Locked Start'));
    plot(lockedEnd(1,1), lockedEnd(1,2),'.',...
        'LineWidth',2,'MarkerSize',30,'Color',[1 0 0], 'DisplayName',sprintf('Locked End'));
end  

%% Find approximate equivalence in station distances
[cellArrayOfPercentStations, ~] = fcn_INTERNAL_calculateApproximateStationEquivalence(cellArrayOfEquidistantStations);

% Calculate the total error via station variation
finalStations = zeros(Npaths,1);
for ith_path = 1:Npaths
    finalStations(ith_path,1) = cellArrayOfEquidistantStations{ith_path}(end,1);
end


%% Perform iterations to seek an average path
% * Project from the reference path to nearby trajectories to find projections
% * Average projections to find new average path. 
% * Clean up average and store results for next iteration or exit

% Initialize variables to hold results
average_path = cell(max_num_iterations,Npaths);
% total_error{max_num_iterations} = [];

stationError = nan(max_num_iterations, 1);
stationError(1) = max(finalStations) - min(finalStations);
previousMaxStation =  max(finalStations);

% % Set iteration 1 values
for ith_path = 1:Npaths
    average_path{1,ith_path} = cellArrayOfEquidistantPaths{ith_path};
end
% total_error{1}  = allDistancesMatrix;

% Start iterations - at each iteration, find orthogonal projections,
% average results, then recalculate the reference traversal
for ith_iteration =2:max_num_iterations 

    % Show user what we are doing?        
    if flag_do_debug
        figure(debug_fig_num)
        title(sprintf('Iteration: %.0f',ith_iteration))
        fprintf(1,'Averaging paths via iteration: %.0d / %.0d \n',ith_iteration,max_num_iterations);
        fprintf(1,'\tInitial maximum station: %.2f  \n',previousMaxStation);    
        % fprintf(1,'\tOld search radius: %.2f  \n',search_radius);
    end
    
    %% For each path, project from reference orthogonally to the others
    % The allTransverseDistances array is organized as
    % (referenceIndex,toIndex)
    allTransverseDistances = fcn_INTERNAL_findTransverseDistancesAlongStation(cellArrayOfEquidistantPaths, cellArrayOfPercentStations);


    %% Sum the corrections
    pathCorrections = cell(Npaths,1);
    for ith_referencePath = 1:Npaths
        for jth_targetPath = 1:Npaths
            if jth_targetPath == 1
                pathCorrections{ith_referencePath,1} = -1/Npaths .* allTransverseDistances{ith_referencePath}(:,jth_targetPath);
            else
                pathCorrections{ith_referencePath,1} = pathCorrections{ith_referencePath,1} - 1/Npaths .* allTransverseDistances{ith_referencePath}(:,jth_targetPath);
             end
        end
    end

    %% Update each path
    newPaths = cell(Npaths,1);
    for ith_path = 1:Npaths
        thisPath = cellArrayOfEquidistantPaths{ith_path};
        thisCorrection = pathCorrections{ith_path,1};
        orthoProjectionVectors = diff(thisPath,1,1)*[0 1; -1 0];
        mag_orthoProjectionVectors = sum(orthoProjectionVectors.^2,2).^0.5;
        unit_orthoProjectionVectors = orthoProjectionVectors./mag_orthoProjectionVectors;

        updatedPath = [lockedStart; thisPath(2:end-1,:)+thisCorrection(2:end-1,:).*unit_orthoProjectionVectors(1:end-1,:); lockedEnd];
        newPaths{ith_path,1} = updatedPath;

        if flag_do_debug
            set(h_updatedPlots(ith_path,1),'Xdata',updatedPath(:,1),'Ydata',updatedPath(:,2));
            pause(0.02);
        end
    end


    %% Call a special function to remove back-tracking behavior
    % Sometimes averaging can produce results that bounce forward and
    % backward. This function removes these "jumps"
    newPathsNoJogs = cell(Npaths,1);
    for ith_path = 1:Npaths
        rawPath = newPaths{ith_path,1};
        rawPathNoJogs = fcn_Path_cleanPathFromForwardBackwardJogs(rawPath, -1);
        newPathsNoJogs{ith_path,1} = rawPathNoJogs;
        if flag_do_debug
            set(h_updatedPlots(ith_path,1),'Xdata',rawPathNoJogs(:,1),'Ydata',rawPathNoJogs(:,2));
        end
    end
    
    
    %% Set up all the paths to have same station interval
    [cellArrayOfEquidistantPaths, cellArrayOfEquidistantStations] = fcn_INTERNAL_resamplePaths(newPathsNoJogs, stationInterval);

    % Find approximate equivalence in station distances
    [cellArrayOfPercentStations, ~] = fcn_INTERNAL_calculateApproximateStationEquivalence(cellArrayOfEquidistantStations);

    % Save results (for later plotting)
    for ith_path = 1:Npaths
        average_path{ith_iteration,ith_path} = cellArrayOfEquidistantPaths{ith_path};
    end

    if flag_do_debug
        for ith_path = 1:Npaths
            finalPath = cellArrayOfEquidistantPaths{ith_path};
            set(h_updatedPlots(ith_path,1),'Xdata',finalPath(:,1),'Ydata',finalPath(:,2));
        end
    end

    % Calculate the total error via station variation
    finalStations = zeros(Npaths,1);
    for ith_path = 1:Npaths
        finalStations(ith_path,1) = cellArrayOfEquidistantStations{ith_path}(end,1);
    end
    
    % Calculate the change in station between the paths
    stationError(ith_iteration,1) = max(finalStations) - min(finalStations);    
    
    % Show results?
    if flag_do_debug
        fprintf(1,'\t Station differences from iteration %.0d to %.0d: %.3f \n',ith_iteration-1, ith_iteration,stationError(ith_iteration,1));
    end

    if exit_tolerance>stationError(ith_iteration,1)
        % Exit the for loop        
        break;
    end
    
end

lastIteration = find(isnan(stationError),1)-1;
if isempty(lastIteration)
    lastIteration = max_num_iterations;
end

% Use final average path to define "true" s-coordinates of the original trajectories, using projection
path_average = cellArrayOfEquidistantPaths{1};

% traversal_average = fcn_Path_convertPathToTraversalStructure(path_average, -1);  
% Calculate final results
% [closestXs, closestYs, closestDistances] = ...
%     fcn_Path_findOrthoScatterFromTraversalToTraversals(...
%     traversal_average.Station, traversal_average, data,...
%     [], [], -1);
closestXs = zeros(length(path_average(:,1)),Npaths);
closestYs = zeros(length(path_average(:,1)),Npaths);
closestDistances = zeros(length(path_average(:,1)),Npaths);

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
    allPointsBeingPlotted = cellArrayOfEquidistantPaths{1};
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
    goodAxis = axis;

    hold on;
    grid on;
    xlabel('X [m]');
    ylabel('Y [m]');

    % Plot the inputs
    colorsPerData = zeros(Npaths,3);
    for ith_path = 1:Npaths
        h_plot = plot(cellArrayOfPaths{ith_path}(:,1),cellArrayOfPaths{ith_path}(:,2),'.-','LineWidth',5,'MarkerSize',20, 'DisplayName',sprintf('Path %.0f',ith_path));
        colorsPerData(ith_path,:) = get(h_plot,'Color');
        set(h_plot,'Color',(colorsPerData(ith_path,:)*0.4 + 0.6*[1 1 1]))
    end

    % Plot the results
    h_updatedPlots = zeros(Npaths,1);
    for ith_path = 1:Npaths
        h_updatedPlots(ith_path,1) = plot(cellArrayOfEquidistantPaths{ith_path}(:,1), cellArrayOfEquidistantPaths{ith_path}(:,2),'.-',...
            'LineWidth',2,'MarkerSize',10,'Color',colorsPerData(ith_path,:)*0.8, 'DisplayName',sprintf('Resampled %.0f',ith_path));
    end

    legend


    if flag_do_debug
        % Plot the path convergence
        figure(2255);
        clf;
        hold on;
        axis(goodAxis);
        xlabel('X [m]')
        ylabel('Y [m]')
       
        % Plot the inputs
        colorsPerData = zeros(Npaths,3);
        for ith_path = 1:Npaths
            h_plot = plot(cellArrayOfPaths{ith_path}(:,1),cellArrayOfPaths{ith_path}(:,2),'.-','LineWidth',5,'MarkerSize',20, 'DisplayName',sprintf('Path %.0f',ith_path));
            colorsPerData(ith_path,:) = get(h_plot,'Color');
            set(h_plot,'Color',(colorsPerData(ith_path,:)*0.4 + 0.6*[1 1 1]))
        end

        % Plot the results to create handles
        h_updatedPlots = zeros(Npaths,1);
        for ith_path = 1:Npaths
            h_updatedPlots(ith_path,1) = plot(cellArrayOfEquidistantPaths{ith_path}(:,1), cellArrayOfEquidistantPaths{ith_path}(:,2),'.-',...
                'LineWidth',2,'MarkerSize',10,'Color',colorsPerData(ith_path,:)*0.8, 'DisplayName',sprintf('Resampled %.0f',ith_path));
        end

        % Animate the results
        for ith_iteration = 1:lastIteration
            title(sprintf('Iteration %.0d of %.0d',ith_iteration,ith_iteration));
            
            for ith_path = 1:Npaths
                finalPath = average_path{ith_iteration,ith_path};
                set(h_updatedPlots(ith_path,1),'Xdata',finalPath(:,1),'Ydata',finalPath(:,2));
            end
            drawnow;
            pause(0.1);
        end
        
        % % Plot the error convergence
        % figure(2233);
        % clf;
        % hold on;
        % xlabel('Index')
        % ylabel('Position change between iterations [m]')
        % for ith_iteration = 1:length(iteration_error_X)
        %     title(sprintf('Iteration %.0d of %.0d',ith_iteration,length(iteration_error_X)));
        %     plot(total_error{ith_iteration});            
        %     pause(0.1);
        % end
        % 
        
        % Plot the error convergence
        figure(2244);
        clf;
        
        semilogy(1:length(stationError),stationError,'.-','Markersize',30);
        xlabel('Iteration')
        ylabel('Change in total station')
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

%% fcn_INTERNAL_findTransverseDistances
function [distanceMatrix, allDistances] = fcn_INTERNAL_findTransverseDistances(reference_path, cellArrayOfPaths)
Npaths = length(cellArrayOfPaths);
allDistances = cell(Npaths,1);
for jth_path = 1:Npaths
    % FORMAT:
    %    St_points = fcn_Path_convertXY2St(referencePath,XY_points,...
    %    (flag_rounding_type), (fig_num));
    St_points = fcn_Path_convertXY2St(reference_path,cellArrayOfPaths{jth_path},...
        ([]), (-1));
    allDistances{jth_path} = real(St_points(:,2));
end
distanceMatrix = cell2mat(allDistances);
end % Ends fcn_INTERNAL_findTransverseDistances  cellArrayOfPaths

%% fcn_INTERNAL_resamplePaths
function [cellArrayOfEquidistantPaths, cellArrayOfEquidistantStations] = fcn_INTERNAL_resamplePaths(cellArrayOfEqualizedPaths, interval)

% Find the station values for each path
cellArrayOfEqualizedStations = fcn_Path_calcPathStation(cellArrayOfEqualizedPaths,-1);

% cellArrayOfEquidistantStations = cellArrayOfEqualizedStations;
% cellArrayOfEquidistantPaths    = cellArrayOfEqualizedPaths;

Npaths = length(cellArrayOfEqualizedPaths);

% Sets up all the paths to have user-defined station interval
cellArrayOfEquidistantPaths = cell(Npaths,1);
cellArrayOfEquidistantStations = cell(Npaths,1);
for ith_path = 1:Npaths
    % What X, Y, and Station values do we have?
    current_X = cellArrayOfEqualizedPaths{ith_path}(:,1);
    current_Y = cellArrayOfEqualizedPaths{ith_path}(:,2);
    current_Station = cellArrayOfEqualizedStations{ith_path};
    current_StationEnd = current_Station(end);

    % Define the stations we want
    reference_station_points    = (0:interval:current_StationEnd)';
    if ~isequal(reference_station_points(end),current_StationEnd)
        longer_reference_station_points    = [reference_station_points; current_StationEnd];
    else
        longer_reference_station_points = reference_station_points;
    end

    % Find X and Y points at wanted stations
    interp_X       = interp1(current_Station,current_X,longer_reference_station_points,'linear','extrap');
    interp_Y       = interp1(current_Station,current_Y,longer_reference_station_points,'linear','extrap');

    % Save results
    cellArrayOfEquidistantPaths{ith_path} = [interp_X, interp_Y];
    cellArrayOfEquidistantStations{ith_path} = longer_reference_station_points;


end
end % fcn_INTERNAL_resamplePaths

%% fcn_INTERNAL_calculateApproximateStationEquivalence
function [cellArrayOfPercentStations, longestStation] = fcn_INTERNAL_calculateApproximateStationEquivalence(cellArrayOfStations)
% Produces a station cell array such that all the input stations have the
% same total path length (unit length)
flag_doDebug = 0;
if 1==flag_doDebug
    figure(48484);
    clf;
    hold on;
end

Npaths = length(cellArrayOfStations);
cellArrayOfPercentStations = cell(Npaths,1);
longestStation = -inf; % Initialize longest path value
for ith_path = 1:Npaths

    % Normalize the path
    cellArrayOfPercentStations{ith_path,1} = cellArrayOfStations{ith_path}./cellArrayOfStations{ith_path}(end,1);

    % Save longest path
    if cellArrayOfStations{ith_path}(end,1)>longestStation
        longestStation = cellArrayOfStations{ith_path}(end,1);
    end
    if 1==flag_doDebug
        figure(48484);
        subplot(2,1,1);
        hold on;
        plot(cellArrayOfStations{ith_path},'-');
        subplot(2,1,2);
        hold on;
        plot(cellArrayOfPercentStations{ith_path},'-');

    end
end

end % Ends fcn_INTERNAL_calculateApproximateStationEquivalence

%% fcn_INTERNAL_findTransverseDistancesAlongStation
function allTransverseDistances = fcn_INTERNAL_findTransverseDistancesAlongStation(cellArrayOfInputPaths, cellArrayOfPercentStations)
Npaths = length(cellArrayOfInputPaths);

% allTransverseDistances is a cell array organized by the "measured" path
% index. Each cell contains NstationsxNpath vectors, where Nstation is the
% number of stations in the "measured" path, and the columns represent each of
% the "from" reference locations. For example, if path 2 is measured as
% distances relative to each of the paths 1, 2, 3, then the results appear
% in cell{2,1} in columns 1, 2, and 3 respectively.

flag_doDebug = 0;
if 1==flag_doDebug
    figure(23457);
    clf;
    hold on;
    axis equal;

    colorsPerData = zeros(Npaths,3);
    h_plotsLocalSearch = zeros(Npaths,1);
    h_plotsLocalReference = zeros(Npaths,1);

    for ith_path = 1:Npaths
        h_plot = plot(cellArrayOfInputPaths{ith_path}(:,1),cellArrayOfInputPaths{ith_path}(:,2),'.-',...
            'LineWidth',5,'MarkerSize',20, 'DisplayName',sprintf('Path %.0f',ith_path));
        colorsPerData(ith_path,:) = get(h_plot,'Color');
        set(h_plot,'Color',(colorsPerData(ith_path,:)*0.4 + 0.6*[1 1 1]))
        
        h_plot = plot(cellArrayOfInputPaths{ith_path}(:,1),cellArrayOfInputPaths{ith_path}(:,2),'.-',...
            'LineWidth',2,'MarkerSize',20, 'DisplayName',sprintf('Focused %.0f',ith_path),'Color',colorsPerData(ith_path,:));
        h_plotsLocalSearch(ith_path,1) = h_plot;

        h_plot = plot(cellArrayOfInputPaths{ith_path}(:,1),cellArrayOfInputPaths{ith_path}(:,2),'-',...
            'LineWidth',1,'MarkerSize',20, 'DisplayName',sprintf('Reference %.0f',ith_path),'Color',colorsPerData(ith_path,:)*0.8);
        h_plotsLocalReference(ith_path,1) = h_plot;
    end

    legend

end

% Find the variation in station distances
allStations = cell2mat(cellArrayOfPercentStations);
stationDifference = diff(allStations);
stationDifference = stationDifference(stationDifference>0);
maxDifference = max(stationDifference);
percentStationRange = 2*maxDifference;

% Initialize distances
allTransverseDistances = cell(Npaths,1);
for ith_path = 1:Npaths
    allTransverseDistances{ith_path,1} = cellArrayOfPercentStations{ith_path}*nan(1,Npaths);
end

stationSteps = (0:percentStationRange:1)';

for jth_station = 1:length(stationSteps)
    thisStation = stationSteps(jth_station);
    startSearch = max(0,thisStation);
    endSearch   = min(1,thisStation+percentStationRange);

    % Fill in search areas for all paths, relative to this start and end
    % station. NOTE: the local reference paths must be larger than this
    % start/end cycle, otherwise it's possible for the hit to be "missed"
    cellArrayOfSearchableAreas = cell(Npaths,1);
    cellArrayOfSearchableIndices = cell(Npaths,1);
    cellArrayOfLocalReferencePath = cell(Npaths,1);

    extendedStartSearch = max(0,thisStation-5*percentStationRange);
    extendedEndSearch = min(1,thisStation+5*percentStationRange);

    for ith_path = 1:Npaths
        thisPath = cellArrayOfInputPaths{ith_path};

        searchStations = cellArrayOfPercentStations{ith_path,1};
        indicesGood = find((searchStations>=startSearch) .* (searchStations<=endSearch));
        cellArrayOfSearchableIndices{ith_path,1} = indicesGood;
        goodPath = thisPath(indicesGood,:); 
        cellArrayOfSearchableAreas{ith_path,1} = goodPath;

        % Find the extended references
        indicesGood = find((searchStations>=extendedStartSearch) .* (searchStations<=extendedEndSearch));
        try
        goodPath = thisPath(indicesGood,:); %#ok<FNDSB>
        catch
            disp('Stop here');
        end
        cellArrayOfLocalReferencePath{ith_path,1} = goodPath;
    
    end

    if 1==flag_doDebug
        figure(23457);
        for ith_path = 1:Npaths
            set(h_plotsLocalSearch(ith_path,1),'XData',cellArrayOfSearchableAreas{ith_path,1}(:,1),'YData',cellArrayOfSearchableAreas{ith_path,1}(:,2));
            set(h_plotsLocalReference(ith_path,1),'XData',cellArrayOfLocalReferencePath{ith_path,1}(:,1),'YData',cellArrayOfLocalReferencePath{ith_path,1}(:,2));
            pause(0.01);
        end
    end

    % For each path, project from reference orthogonally to the others
    % The allTransverseDistances array is organized as
    % (referenceIndex,toIndex)

    % Loop through each reference path, saving the distances of each of the
    % OTHER paths relative to this reference.
    for ith_pathReference = 1:Npaths
        % Find distance from this path to the others
        [~, tDistance] = fcn_INTERNAL_findTransverseDistances(cellArrayOfLocalReferencePath{ith_pathReference}, cellArrayOfSearchableAreas);

        % Save results of each of the others
        % allTransverseDistances is a cell array organized by the "measured" path
        % index. Each cell contains NstationsxNpath vectors, where Nstation is the
        % number of stations in the "measured" path, and the columns represent each of
        % the "from" reference locations. For example, if path 2 is measured as
        % distances relative to each of the paths 1, 2, 3, then the results appear
        % in cell{2,1} in columns 1, 2, and 3 respectively.
        for jth_targetPath = 1:Npaths
            indicesSearched = cellArrayOfSearchableIndices{jth_targetPath};
            allTransverseDistances{jth_targetPath}(indicesSearched, ith_pathReference) = tDistance{jth_targetPath};
        end
    end
end

end % Ends fcn_INTERNAL_findTransverseDistancesAlongStation