function [stations, differences] = fcn_Path_calcPathStation(path,varargin)
% fcn_Path_calcPathStation
% finds the station for an input path, where station is the distance along
% the path at each point in the path. Also calculates the difference in
% positions along the path from one point to the next. Both station and
% difference are zero for the first point in the path. Both 2D and 3D paths
% are allowed. 
%
% FORMAT: 
%       [stations, differences] = fcn_Path_calcPathStation(path, (fig_num))
%
% INPUTS:
%
%      path: an N x 2 vector of [X Y] positions, with N>=2 OR
%            an N x 3 vector of [X Y Z] positions, with N>=2
%
%      NOTE: if the path input is a cell type, then the code proceeds as if
%      a cell array of path types has been given, and the values output are
%      cell arrays of stations and differences for each path.
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%      stations: the XYZ distance as an N x 1 vector, representing
%           the distance traveled up to the current point (starting with 0
%           at the first point)
% 
%      differences: a [N x 2] or [N x 3] array that is the change in X and Y
%            (front-padded with [0 0 (0)])
%
%      (OPTIONAL INPUTS)
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%      
%       See the script:
%       script_test_fcn_Path_calcPathStation.m for a full
%       test suite. 
%
% This function was written on 2020_11_12 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2025_06_29
% -- first write of the code. Moving code out of traversals types

% TO-DO
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
        narginchk(1,MAX_NARGIN);

        % Check the Path variables
        if ~iscell(path)
            fcn_DebugTools_checkInputsToFunctions(path, 'path2or3D');
        else
            for ith_path = 1:length(path)
                thisPath = path{ith_path};
                fcn_DebugTools_checkInputsToFunctions(thisPath, 'path2or3D');
            end
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
if flag_do_debug
    fig_debug = 8383; %#ok<NASGU>
end


%% Solve for the circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Did user give a cell array?
if iscell(path)
    cellPath = path;
    Npaths = length(path);
    flag_returnCells = 1;
else
    Npaths = 1;
    cellPath = {path};
    flag_returnCells = 0;
end

% Initialize cell arrays
cellDifferences = cell(Npaths,1);
cellStations = cell(Npaths,1);

% Loop through paths, calculating stations and differences for each
for ith_path = 1:Npaths
    thisPath = cellPath{ith_path};

    % Make sure all paths have same length
    if 1==ith_path
        pointDimension = length(thisPath(1,:));
    else
        if(length(thisPath(1,:))~=pointDimension)
            error('Paths in cell array cannot mix 2D and 3D paths.');
        end
    end

    % Calculate the station differences, and force the diff operation to occur
    % along rows
    % diff(X,N,DIM) is the Nth difference function along dimension DIM.
    if length(thisPath(1,:))==3
        thisDifferences = [[0 0 0]; diff(thisPath,1,1)];
    else
        thisDifferences = [[0 0]; diff(thisPath,1,1)];
    end

    % Calculate the stations
    thisStations = cumsum(real(sum(thisDifferences.^2,2).^0.5),'omitnan');

    cellDifferences{ith_path} = thisDifferences;
    cellStations{ith_path} = thisStations;

end

% Save outputs
if 1==flag_returnCells
    stations = cellStations;
    differences = cellDifferences;
else
    stations = cellStations{1};
    differences = cellDifferences{1};
end


%% Any debugging?
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


    %%%%%
    % PLOT STATIONS
    subplot(2,1,1);

    dimension_of_points = 2;

    % Find size of plotting domain
    allPointsBeingPlotted = [];
    cellIndices = cell(Npaths, 1);
    for ith_cell = 1:Npaths
        cellIndices{ith_cell,1} = (1:length(cellStations{ith_cell}(:,1)))';
        allPointsBeingPlotted = [allPointsBeingPlotted; cellIndices{ith_cell} cellStations{ith_cell}]; %#ok<AGROW>
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

    xlabel('Index');
    ylabel('Station [m]');

    % Plot the stations
    plotColors = zeros(Npaths,3);
    for ith_cell = Npaths:-1:1
        h_plot = plot(cellIndices{ith_cell,1}, cellStations{ith_cell,1},'.-','Linewidth',1+ith_cell*2,'Markersize',25,'DisplayName',sprintf('Path station %.0d',ith_cell));
        plotColors(ith_cell,:) = get(h_plot,'Color');
    end
    title('Stations');
    legend('Location','northwest');

    %%%%%
    % PLOT DIFFERENCES
    subplot(2,1,2);

    dimension_of_points = 2;

    % Find size of plotting domain
    allPointsBeingPlotted = [];
    for ith_cell = 1:Npaths
        if length(path(1,:))==3
            allPointsBeingPlotted = [allPointsBeingPlotted; cellIndices{ith_cell,1} cellDifferences{ith_cell,1}(:,1); cellIndices{ith_cell,1} cellDifferences{ith_cell,1}(:,2); cellIndices{ith_cell,1} cellDifferences{ith_cell,1}(:,3)]; %#ok<AGROW>
        else
            allPointsBeingPlotted = [allPointsBeingPlotted; cellIndices{ith_cell,1} cellDifferences{ith_cell,1}(:,1); cellIndices{ith_cell,1} cellDifferences{ith_cell,1}(:,2)]; %#ok<AGROW>
        end
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

    xlabel('Index');
    ylabel('Differences [m]');

    % Plot the differences
    for ith_cell = Npaths:-1:1
        if length(path(1,:))==3
            plot(cellIndices{ith_cell,1}, cellDifferences{ith_cell,1}(:,1),'.-','Linewidth',1+ith_cell*2,'Markersize',25,'DisplayName',sprintf('Path %.0d X',ith_cell), 'Color',plotColors(ith_cell,:));
            plot(cellIndices{ith_cell,1}, cellDifferences{ith_cell,1}(:,2),'.--','Linewidth',1+ith_cell*2,'Markersize',25,'DisplayName',sprintf('Path %.0d Y',ith_cell), 'Color',plotColors(ith_cell,:));
            plot(cellIndices{ith_cell,1}, cellDifferences{ith_cell,1}(:,3),'.:','Linewidth',1+ith_cell*2,'Markersize',25,'DisplayName',sprintf('Path %.0d Z',ith_cell), 'Color',plotColors(ith_cell,:));
        else
            plot(cellIndices{ith_cell,1}, cellDifferences{ith_cell,1}(:,1),'.-','Linewidth',1+ith_cell*2,'Markersize',25,'DisplayName',sprintf('Path %.0d X',ith_cell), 'Color',plotColors(ith_cell,:));
            plot(cellIndices{ith_cell,1}, cellDifferences{ith_cell,1}(:,2),'.--','Linewidth',1+ith_cell*2,'Markersize',25,'DisplayName',sprintf('Path %.0d Y',ith_cell), 'Color',plotColors(ith_cell,:));
        end
    end
    title('Differences');
    legend('Location','northwest');

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