function offset_paths = fcn_Path_fillOffsetPathsAboutPath(reference_path, offsets, varargin)
% fcn_Path_fillOffsetPathsAboutPath
% fills in an array of paths about a reference traversal at
% user-defined offset distances.
%
% FORMAT:
%
%      offset_paths = ...
%      fcn_Path_fillOffsetPathsAboutPath(...
%            reference_path,...
%            offsets,...
%            (flag_rounding_type),
%            (fig_num));
%
% INPUTS:
%
%      reference_path: a N x 2 or N x 3 set of coordinates representing the 
%      [X Y] or [X Y Z] coordinates, in sequence, of a path
%
%      offsets: a scalar value or Mx1 array representing the offsets to be
%      used in calculating the resulting paths, producing one path
%      per offset value. Positive values follow the cross-product notation,
%      namely if a path is running from left-to-right, then a positive
%      offset would be above it, negative would be below.
%
%      (OPTIONAL INPUTS)
%
%      flag_rounding_type: determines type of projection, and is passed
%      into fcn_Path_findOrthogonalTraversalVectorsAtStations. See that
%      function for more explanation. Default (empty) is used unless
%      changed via this input.
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
%
% OUTPUTS:
%
%      offset_paths: a structure containing the resulting paths
%      that are generated, one for each given offset
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_Path_calcSingleTraversalStandardDeviation
%      fcn_Path_findOrthogonalTraversalVectorsAtStations
%      fcn_Path_convertPathToTraversalStructure
%      fcn_Path_plotTraversalsXY
%
% EXAMPLES:
%
%     See the script: script_test_fcn_Path_fillOffsetPathsAboutPath
%     for a full test suite.
%
% This function was written on 2021_01_03 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2021_01_24
% -- first write of the code, using
% fcn_Path_fillRandomTraversalsAboutTraversal as a template
% 2022_01_03
% -- minor updates to comments
% 2022_08_20
% -- allow empty figure argument to avoid plotting
% 2023_09_17 by S. Brennan
% -- added flag_rounding_type to inputs
% -- fixed some comments
% 2025_06_23 - S. Brennan
% -- Updated debugging and input checks
% 2025_07_01 - S. Brennan
% -- removed traversal type to convert function to path type, using
% OffsetTraversalsABoutTraversal as template

% TO-DO
% (none)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 4; % The largest Number of argument inputs to the function
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

        % Check the central_path input
        fcn_DebugTools_checkInputsToFunctions(reference_path, 'path2or3D');

        % Check the offsets input (looks like a station type)
        fcn_DebugTools_checkInputsToFunctions(offsets, 'station');

    end
end
% Does user want to specify flag_rounding_type?
flag_rounding_type = []; % Default is to not do any plotting
if 3 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        flag_rounding_type = temp;
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
    fig_debug = 12312; %#ok<NASGU>
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

% Set defaults
Station_reference = fcn_Path_calcPathStation(reference_path,-1);

Nstations = length(Station_reference);
Noffsets = length(offsets);

% the default number of trajectories to use
num_trajectories = length(offsets(:,1));

% the default number of points to use
num_points = Nstations;



%% Fill in the array of reference station points
% Each column corresponds to one trajectory.

% Duplicate the station reference points across all trajectories
reference_station_points = Station_reference*ones(1,num_trajectories);


%% Fill in the array of offset distances.
% Each column corresponds to one trajectory. The matrix multiplication
% below is (Num_points x 1 ) X (1 x offsets) = Numpoints x offsets
offsets_from_reference = ones(num_points,1)*offsets';

%% Find projection from reference orthogonally

% Find the unit normal vectors at each of the station points
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalPathVectorsAtStations(...
    reference_station_points(:,1),reference_path,flag_rounding_type, -1);
unit_vectors = unit_normal_vector_end - unit_normal_vector_start;

%% Perform iterations over trajectories
offset_paths = cell(Noffsets,1);
for ith_trajectory =1:num_trajectories
    % Show user what we are doing?
    if flag_do_debug
        fprintf(1,'Generating offset trajectory: %.0d / %.0d \n',ith_trajectory,num_trajectories);
    end

    % Calculate random path and traversal and save into final structure
    offset_path = unit_normal_vector_start + unit_vectors.*offsets_from_reference(:,ith_trajectory);
    offset_paths{ith_trajectory} = offset_path;
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
    figure(fig_num);
    hold on;

    % Plot the reference trajectory first
    plot(reference_path(:,1),reference_path(:,2),'b.-','Linewidth',4,'Markersize',20,'DisplayName','Reference path');

    % Plot the ofset results
    fcn_Path_plotPathsXY(offset_paths,fig_num);
    title('Reference traversal and offset paths');
    xlabel('X [m]');
    ylabel('Y [m]');

    % Add a legend
    legend('Reference traversal', 'Offset traversals');

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
