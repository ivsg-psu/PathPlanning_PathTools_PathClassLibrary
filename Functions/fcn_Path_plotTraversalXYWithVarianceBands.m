function fcn_Path_plotTraversalXYWithVarianceBands(reference_traversal, varargin)
% fcn_Path_plotTraversalXYWithVarianceBands
% Plots a traversal with a variance band around the path
%
% FORMAT:
%
%      fcn_Path_plotTraversalXYWithVarianceBands(...
%            reference_traversal,...
%            (std_deviation),...
%            (fig_num));
%
% INPUTS:
%
%      reference_traversal: the traversal that is being used for randomly
%      generating all the other traversals
%
%      (OPTIONAL INPUTS)
%
%      std_deviation: a positive value representing the standard deviation
%      to use in calculating the deviation distance. Note: the default is
%      to use the variance in angles along the reference_traversal
%      multiplied by the average segment length in the reference traversal.
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%      (none)
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_Path_convertPathToTraversalStructure
%      fcn_Path_calcSingleTraversalStandardDeviation
%      fcn_Path_findOrthogonalTraversalVectorsAtStations
%      fcn_Path_plotTraversalXYWithUpperLowerBands
%
% EXAMPLES:
%
%     See the script: script_test_fcn_Path_plotTraversalXYWithVarianceBands
%     for a full test suite.
%
% This function was written on 2021_01_05 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2021_01_05:
% -- wrote the code originally
% 2021_01_08:
% -- added input checking
% 2022_01_03:
% -- corrected dependency list of functions
% -- broke out fcn_Path_plotTraversalXYWithUpperLowerBands into
% separate functionality (it is very useful for other functions!)
% 2025_06_23 - S. Brennan
% -- Updated debugging and input checks

% TO-DO
% (none)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==3 && isequal(varargin{end},-1))
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
        narginchk(1,3);

        % Check the reference_traversal input
        fcn_DebugTools_checkInputsToFunctions(reference_traversal, 'traversal');

    end
end

% Does the user want to specify standard deviation?
std_deviation = fcn_Path_calcSingleTraversalStandardDeviation(reference_traversal); % Default
if 2 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        std_deviation = temp;
    end
end



% Does user want to show the plots?
flag_do_plots = 0; % Default is to NOT make a plot
fig_num = [];
if (0==flag_max_speed) && (3 == nargin) 
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        fig_num = temp;
        figure(fig_num);
        flag_do_plots = 1;
    end
else
    if flag_do_debug
        fig_debug = 8838; %#ok<NASGU>
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

% Grab key variables
Station_reference = reference_traversal.Station;
Nstations = length(Station_reference(:,1));

% the default number of points to use
num_points = Nstations;

%% Fill in the array of stations.
reference_station_points = Station_reference;

%% Fill in the array of offset distances.
% Each column corresponds to one trajectory.
offsets_from_reference = std_deviation*ones(num_points,1);

%% Find offsets from trajectory
% Set the projection type (see help in function below for details)
flag_rounding_type = 4; % This averages the projection vectors along segments

% Find the unit normal vectors at each of the station points
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalTraversalVectorsAtStations(...
    reference_station_points(:,1),reference_traversal,flag_rounding_type);
unit_vectors = unit_normal_vector_end - unit_normal_vector_start;

%% Calculate random path and traversal and save into final structure
top_path = unit_normal_vector_start + unit_vectors.*offsets_from_reference(:,1);
bottom_path = unit_normal_vector_start - unit_vectors.*offsets_from_reference(:,1);
upper_traversal = fcn_Path_convertPathToTraversalStructure(top_path);
lower_traversal = fcn_Path_convertPathToTraversalStructure(bottom_path);



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
    fcn_Path_plotTraversalXYWithUpperLowerBands( reference_traversal, upper_traversal, lower_traversal, fig_num);
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends main function



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


