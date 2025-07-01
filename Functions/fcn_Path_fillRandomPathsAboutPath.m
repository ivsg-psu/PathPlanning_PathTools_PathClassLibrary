function random_paths = fcn_Path_fillRandomPathsAboutPath(reference_path, varargin)
% fcn_Path_fillRandomPathsAboutPath
% fills in random paths about a reference path. Points are
% generated via orthogonal projection using random normal distribution with
% either a default variance or optional user-defined variance. The station
% points can also be user-specified as randomly distributed uniformly, or
% default to the station points in the reference_path if no optional
% inputs are given. The first and last stations are forced to be the same
% as the reference_path to prevent the route from randomly becoming
% shorter with repeated calls to this function.
%
% FORMAT:
%
%      random_paths = ...
%      fcn_Path_fillRandomPathsAboutPath(...
%            reference_path,...
%            (std_deviation),...
%            (num_paths),...
%            (num_points),...
%            (flag_generate_random_stations),...
%            (spatial_smoothness),...
%            (fig_num));
%
% INPUTS:
%
%      reference_path: the path that is being used for randomly generating
%      all the other paths. A path is a N x 2 or N x 3 set of coordinates
%      representing the [X Y] or [X Y Z] coordinates, in sequence, of a
%      path.
%
%      (OPTIONAL INPUTS)
%
%      std_deviation: a positive value representing the standard deviation
%      to use in calculating the deviation distance (default is to use the
%      variance in angles along the reference_path multiplied by the
%      average segment length in the reference path)
%
%      num_paths: an integer to specify how many paths to
%      generate. (default is 1)
%
%      num_points: an integer to specify the number of points to produce in
%      each random trajectory. (default is same as the reference_trajectory)
%      Note: num_points is not used if flag_gnerate_random_stations is
%      equal to 1, as this flag forces the station points and their number
%      to match the reference.
%
%      flag_generate_random_stations: a flag to specify whether or not to
%      use random stations along the reference trajectory. Options include:
%            1: Use random station distances (default)
%            0: Use the station distances projected at vertices of the
%            reference path
%
%      spatial_smoothness: a distance variable representing the 1/length of
%      expected frequency variation in the data. This represents the cutoff
%      frequency of random swerving. For example: if the spatial smoothness
%      is 20 meters, then swerves tighter than 20 meters would be
%      increasingly rare. This smoothness is implemented via a zero-phase
%      filtering of the random output that maintains the desired standard
%      deviation. Default value is 40 meters. Note: if the reference
%      trajectory is too corsely sampled for a desired spatial smoothness,
%      a warning will be given.
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%      random_paths: a structure containing the resulting paths
%      that are randomly generated about the reference_path
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_Path_calcSinglePathStandardDeviation
%      fcn_Path_calcPathStation
%      fcn_Path_findOrthogonalPathVectorsAtStations
%      fcn_Path_plotTraversalsXY
%
% EXAMPLES:
%
%     See the script: script_test_fcn_Path_fillRandomPathsAboutPath
%     for a full test suite.
%
% This function was written on 2021_01_03 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2021_01_03:
% -- wrote the code originally
% 2021_01_07
% -- added functionalized input checking
% -- fixed typos in comments, plotting at end
% 2021_01_09
% -- fixed function calls that were misnamed due to class edits
% -- updated dependencies
% 2025_06_23 - S. Brennan
% -- Updated debugging and input checks
% 2025_07_01 - S. Brennan
% -- Removed traversal input type and replaced with cell array of paths

% TO-DO
% (none)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 7; % The largest Number of argument inputs to the function
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

        % Check the reference_path variables
        fcn_DebugTools_checkInputsToFunctions(reference_path, 'path2or3D');

    end
end

%      random_paths = ...
%      fcn_Path_fillRandomPathsAboutPath(...
%            reference_path,...
%            (std_deviation),...
%            (num_paths),...
%            (num_points),...
%            (flag_generate_random_stations),...
%            (spatial_smoothness),...
%            (fig_num));

% Does the user want to specify standard deviation?
std_deviation = fcn_Path_calcSinglePathStandardDeviation(reference_path, -1); % the default standard deviation
if 2 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        std_deviation = temp;
    end
end


% Does the user want to specify num_paths?
num_paths = 1; % the default number of paths to use
if 3 <= nargin
    temp = varargin{2};
    if ~isempty(temp)
        num_paths = temp;
    end
end

% Does the user want to specify num_points?
Station_reference = fcn_Path_calcPathStation(reference_path,-1);
Nstations = length(Station_reference);
num_points = Nstations; % the default number of points to use
if 4 <= nargin
    temp = varargin{3};
    if ~isempty(temp)
        num_points = temp;
    end
end

% Does the user want to specify flag_generate_random_stations?
flag_generate_random_stations = 1; % the default on generating random station locations
if 5 <= nargin
    temp = varargin{4};
    if ~isempty(temp)
        flag_generate_random_stations = temp;
    end
end

% Does user want to specify spatial_smoothness?
spatial_smoothness = 40; % the default spatial smoothness Units are meters
if 6 <= nargin
    temp = varargin{5};
    if ~isempty(temp)
        spatial_smoothness = temp;
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
    fig_debug = 34324; %#ok<NASGU>
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

% Define useful variables used in several sections below
maxStation = Station_reference(end);

%% Fill in the array of reference station points
% Each column corresponds to one trajectory.

if flag_generate_random_stations
    % Randomize the station reference points
    reference_station_points = maxStation*rand(num_points,num_paths);
    % Sort the stations
    reference_station_points = sort(reference_station_points,1);
else
    % Duplicate the station reference points across all paths
    reference_station_points = Station_reference*ones(1,num_paths);
end

% Force the first and last points to be the same as the reference path
reference_station_points(1,:) = Station_reference(1);
reference_station_points(end,:) = Station_reference(end);

%% Fill in the array of offset distances.
% Each column corresponds to one trajectory.
offsets_from_reference = std_deviation*randn(num_points,num_paths);

%% Spatially smooth the offsets
spatial_sampling_frequency = Nstations/maxStation; % Units are samples/meter
spatial_smoothness_frequency = 1/spatial_smoothness;  % Units are 1/meter
freq_ratio = spatial_smoothness_frequency/(2*spatial_sampling_frequency);
if freq_ratio>1
    fprintf(1,'Spatial smoothness not possible with this reference trajectory,\n');
    fprintf(1,'which has a Nyquist sampling equivalent to a spatial smoothness \n');
    fprintf(1,'calculated to be approximately: %.2f meters.\n',1/(2*spatial_sampling_frequency));
    fprintf(1,'User-requested smoothness was: %.2f meters \n',spatial_smoothness);
    fprintf(1,'The smoothing of random variations is not possible \n');
    fprintf(1,'if requested smoothness is lower than Nyquist sampling smoothness.\n');
    warning('The resulting random paths may be very irregular.');
else
    % Design a filter
    [B,A] = butter(2,freq_ratio);

    % Loop through each trajectory's offset column
    for ith_path =1:num_paths
        % Grab the raw offsets
        raw_offsets = offsets_from_reference(:,ith_path);
        raw_offsets(1) = 0;
        raw_offsets(end) = 0;

        % Smooth the offsets
        filtered_offsets = filtfilt(B,A,raw_offsets);

        % Fix the amplitude to match desired standard deviation
        filt_std = std(filtered_offsets);
        filtered_offsets = std_deviation/filt_std * filtered_offsets;

        % Save results
        offsets_from_reference(:,ith_path) = filtered_offsets;

        % For debugging:
        if 1==0
            figure(33737);
            clf; hold on;
            plot(raw_offsets,'r');
            plot(filtered_offsets,'b');
            disp('Pause here');
        end
    end
end


%% Perform iterations over paths
random_paths = cell(num_paths,1);
for ith_path =1:num_paths
    % Show user what we are doing?
    if flag_do_debug
        fprintf(1,'Generating random path: %.0d / %.0d \n',ith_path,num_paths);
    end


    %% For each path, project from reference orthogonally
    % Set the projection type (see help in function below for details)
    flag_rounding_type = 4; % This averages the projection vectors along segments

    % Find the unit normal vectors at each of the station points
    [unit_normal_vector_start, unit_normal_vector_end] = ...
        fcn_Path_findOrthogonalPathVectorsAtStations(...
        reference_station_points(:,ith_path),reference_path,flag_rounding_type, -1);
    unit_vectors = unit_normal_vector_end - unit_normal_vector_start;

    %% Calculate random path and path and save into final structure
    random_path = unit_normal_vector_start + unit_vectors.*offsets_from_reference(:,ith_path);
    random_paths{ith_path} = random_path;

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
    plot(reference_path(:,1),reference_path(:,2),'b.-','Linewidth',4,'Markersize',20);

    % Plot the random results
    fcn_Path_plotPathsXY(random_paths,fig_num);
    title('Reference path and random paths');
    xlabel('X [m]');
    ylabel('Y [m]');


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
