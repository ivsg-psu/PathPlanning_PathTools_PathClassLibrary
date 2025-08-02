function [changeInAngles, edgeLengths] = fcn_Path_calcDiffAnglesBetweenPathSegments(pathVerticesXY,varargin)
% fcn_Path_calcDiffAnglesBetweenSegments
% Calculates the change in angles between path segments. If there are N
% points in the path, there are N-1 segments and thus N-2 angles between
% segments.
%
% Note that this method uses the dot product and the cross product to avoid
% errors caused by angle rollover that occur with arctan calculations. The
% arctan calculation gives the wrong answer with path segments are pointed
% near or at the -180 degree crossover point.
%
% FORMAT: 
%
%       [changeInAngles, edgeLengths] = fcn_Path_calcDiffAnglesBetweenPathSegments(pathVerticesXY,(fig_num))
%
% INPUTS:
%
%      pathVerticesXY: an N x 2 vector with [X Y] data in each row. N must be >= 3.
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
%      changeInAngles: an (N-2) x 1 vector of the change in angles in radians
%
%      edgeLengths: an (N-1) x 1 vector of the lengths of the edges
%
% DEPENDENCIES
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%      
%       See the script:
%       script_test_fcn_Path_calcDiffAnglesBetweenPathSegments.m for a full
%       test suite.
%
% This function was written on 2021_01_03 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2021_01_03
% -- first writing of the code
% 2021_01_06
% -- fixed typos in the comments
% 2021_01_07
% -- fixed typos in the comments, minor header clean-ups
% 2025_06_23 - S. Brennan
% -- Updated debugging and input checks
% 2025_08_02 - S. Brennan
% - In fcn_Path_calcDiffAnglesBetweenSegments
%   % * Minor reordering of code steps to pass out edge lengths
%   % * Allows speed up other codes using same values
%   % * renamed variables for clarity
%   % * updated docstrings
%   % * Removed fill of debug_fig_num, as this is in the new header

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

        % Check the pathVerticesXY variables
        fcn_DebugTools_checkInputsToFunctions(pathVerticesXY, 'paths');
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

% Note that this method uses the dot product and the cross product to avoid
% errors caused by angle rollover. One could use the arctan method of
% calculating angles, which works in general but fails for paths that point
% straight to the left due to the cross-over point for the atan2
% calculation. For example, the following code does NOT always work:
%
%     % angles = atan2(path_average(2:end,2)-path_average(1:end-1,2),path_average(2:end,1)-path_average(1:end-1,1));
%     % angles = [angles; angles(end)];  % Pad the last point twice
%     % changeInAngles2 = abs(diff(angles));

edgeVectors = pathVerticesXY(2:end,:)-pathVerticesXY(1:end-1,:);
edgeLengths = sum(edgeVectors.^2,2).^0.5;

incoming_vector = edgeVectors(1:end-1,:);
outgoing_vector = edgeVectors(2:end,:);
% OLD VERSION:
% incoming_vector = pathVerticesXY(2:end-1,:)-pathVerticesXY(1:end-2,:);
% outgoing_vector = pathVerticesXY(3:end,:)-pathVerticesXY(2:end-1,:);
incoming_dot_outgoing = sum(incoming_vector.*outgoing_vector,2); % Do the dot product
incoming_cross_outgoing = crossProduct(incoming_vector,outgoing_vector);

zeroCross = find(0==incoming_cross_outgoing);
if incoming_dot_outgoing(zeroCross)~=0
    incoming_cross_outgoing(zeroCross) = 1*sign(incoming_dot_outgoing(zeroCross));
end

incoming_mag = edgeLengths(1:end-1,:);
outgoing_mag = edgeLengths(2:end,:);
% OLD VERSION:
% incoming_mag = sum(incoming_vector.^2,2).^0.5;
% outgoing_mag = sum(outgoing_vector.^2,2).^0.5;

% Calculate the change in angle using the dot product
%changeInAngles = [acos(incoming_dot_outgoing./(incoming_mag.*outgoing_mag)); 0];
changeInAngles = sign(incoming_cross_outgoing).*(acos(incoming_dot_outgoing./(incoming_mag.*outgoing_mag)));


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
    temp_h = figure(fig_num); %#ok<NASGU>
    % flag_rescale_axis = 0;
    % if isempty(get(temp_h,'Children'))
    %     flag_rescale_axis = 1;
    % end
    
    % Is this 2D or 3D?
    % dimension_of_points = length(changeInAngles(1,:));

    % Find size of plotting domain
    allPointsBeingPlotted = [(1:length(changeInAngles))' changeInAngles*180/pi];
    max_plotValues = max(allPointsBeingPlotted);
    min_plotValues = min(allPointsBeingPlotted);
    sizePlot = max(max_plotValues) - min(min_plotValues);
    nudge = sizePlot*0.006; %#ok<NASGU>

    % Find size of plotting domain
    % if flag_rescale_axis
    %     percent_larger = 0.3;
    %     axis_range = max_plotValues - min_plotValues;
    %     if (0==axis_range(1,1))
    %         axis_range(1,1) = 2/percent_larger;
    %     end
    %     if (0==axis_range(1,2))
    %         axis_range(1,2) = 2/percent_larger;
    %     end
    %     if dimension_of_points==3 && (0==axis_range(1,3))
    %         axis_range(1,3) = 2/percent_larger;
    %     end
    % 
    %     % Force the axis to be equal?
    %     if 1==0
    %         min_valuesInPlot = min(min_plotValues);
    %         max_valuesInPlot = max(max_plotValues);
    %     else
    %         min_valuesInPlot = min_plotValues;
    %         max_valuesInPlot = max_plotValues;
    %     end
    % 
    %     % Stretch the axes
    %     stretched_min_vertexValues = min_valuesInPlot - percent_larger.*axis_range;
    %     stretched_max_vertexValues = max_valuesInPlot + percent_larger.*axis_range;
    %     axesTogether = [stretched_min_vertexValues; stretched_max_vertexValues];
    %     newAxis = reshape(axesTogether, 1, []);
    %     axis(newAxis);
    % 
    % end
    % goodAxis = axis;

    hold on;
    grid on;

    xlabel('Index');
    ylabel('Angle [deg]');

    % Plot the angle differences
    plot(changeInAngles*180/pi,'k.-','Linewidth',3,'Markersize',25);

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
%% Calculate cross products
function result = crossProduct(v,w)
result = v(:,1).*w(:,2)-v(:,2).*w(:,1);
end

