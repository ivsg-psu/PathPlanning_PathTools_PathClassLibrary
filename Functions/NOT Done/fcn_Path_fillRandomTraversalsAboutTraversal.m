function random_traversals = fcn_Path_fillRandomTraversalsAboutTraversal(reference_traversal, varargin)
% fcn_Path_fillRandomTraversalsAboutTraversal
% fills in random traversals about a reference traversal. Points are
% generated via orthogonal projection using random normal distribution with
% either a default variance or optional user-defined variance. The station
% points can also be user-specified as randomly distributed uniformly, or
% default to the station points in the reference_traversal if no optional
% inputs are given. The first and last stations are forced to be the same
% as the reference_traversal to prevent the route from randomly becoming
% shorter with repeated calls to this function.
%
% FORMAT:
%
%      random_traversals = ...
%      fcn_Path_fillRandomTraversalsAboutTraversal(...
%            reference_traversal,...
%            (std_deviation),...
%            (num_trajectories),...
%            (num_points),...
%            (flag_generate_random_stations),...
%            (spatial_smoothness),...
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
%      to use in calculating the deviation distance (default is to use the
%      variance in angles along the reference_traversal multiplied by the
%      average segment length in the reference traversal)
%
%      num_trajectories: an integer to specify how many trajectories to
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
%            reference traversal
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
%      random_trajectories: a structure containing the resulting traversals
%      that are randomly generated about the reference_traversal
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
%     See the script: script_test_fcn_Path_fillRandomTraversalsAboutTraversal
%     for a full test suite.
%
% This function was written on 2021_01_03 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
%     2021_01_03:
%     -- wrote the code originally
%     2021_01_07
%     -- added functionalized input checking
%     -- fixed typos in comments, plotting at end
%     2021_01_09
%     -- fixed function calls that were misnamed due to class edits
%     -- updated dependencies

flag_do_debug = 0; % Flag to show the results for debugging
flag_do_plots = 0; % % Flag to plot the final results
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
end


%% check input arguments
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

%      random_trajectories = ...
%      fcn_Path_fillRandomTraversalsAboutTraversal(...
%            reference_traversal,...
%            (std_deviation),...
%            (num_trajectories),...
%            (num_points),...
%            (flag_generate_random_stations),...
%            (spatial_smoothness),...
%            (fig_num));

% Check inputs?
if flag_check_inputs
    % Are there the right number of inputs?
    if nargin < 1 || nargin > 7
        error('Incorrect number of input arguments')
    end
        
    % Check the reference_traversal variables
    fcn_DebugTools_checkInputsToFunctions(reference_traversal, 'traversal');
    
end

Station_reference = reference_traversal.Station;
Nstations = length(Station_reference);


%% Set defaults

% the default standard deviation
std_deviation = fcn_Path_calcSingleTraversalStandardDeviation(reference_traversal);

% the default number of trajectories to use
num_trajectories = 1;

% the default number of points to use
num_points = Nstations;

% the default on generating random station locations
flag_generate_random_stations = 1;

% the default spatial smoothness
spatial_smoothness = 40; % Units are meters

%% Check for variable argument inputs (varargin)
% Does the user want to specify standard deviation?
if 2 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        std_deviation = temp;
    end
end


% Does the user want to specify num_trajectories?
if 3 <= nargin
    temp = varargin{2};
    if ~isempty(temp)
        num_trajectories = temp;
    end
end


% Does the user want to specify num_points?
if 4 <= nargin
    temp = varargin{3};
    if ~isempty(temp)
        num_points = temp;
    end
end


% Does the user want to specify flag_generate_random_stations?
if 5 <= nargin
    temp = varargin{4};
    if ~isempty(temp)
        flag_generate_random_stations = temp;
    end
end

% Does user want to specify spatial_smoothness?
if 6 <= nargin
    temp = varargin{5};
    if ~isempty(temp)
        spatial_smoothness = temp;
    end
end

% Does user want to show the plots?
if 7 == nargin
    fig_num = varargin{6};
    figure(fig_num);
    flag_do_plots = 1;
else
    if flag_do_debug
        fig = figure;
        fig_num = fig.Number;
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

% Define useful variables used in several sections below
maxStation = Station_reference(end);
    
%% Fill in the array of reference station points
% Each column corresponds to one trajectory.

if flag_generate_random_stations
    % Randomize the station reference points
    reference_station_points = maxStation*rand(num_points,num_trajectories);
    % Sort the stations
    reference_station_points = sort(reference_station_points,1);
else
    % Duplicate the station reference points across all trajectories
    reference_station_points = Station_reference*ones(1,num_trajectories);
end

% Force the first and last points to be the same as the reference traversal
reference_station_points(1,:) = Station_reference(1);
reference_station_points(end,:) = Station_reference(end);

%% Fill in the array of offset distances.
% Each column corresponds to one trajectory.
offsets_from_reference = std_deviation*randn(num_points,num_trajectories);

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
    warning('The resulting random trajectories may be very irregular.');
else
    % Design a filter
    [B,A] = butter(2,freq_ratio);
    
    % Loop through each trajectory's offset column
    for ith_trajectory =1:num_trajectories
        % Grab the raw offsets
        raw_offsets = offsets_from_reference(:,ith_trajectory);
        raw_offsets(1) = 0;
        raw_offsets(end) = 0;
        
        % Smooth the offsets
        filtered_offsets = filtfilt(B,A,raw_offsets);
        
        % Fix the amplitude to match desired standard deviation
        filt_std = std(filtered_offsets);
        filtered_offsets = std_deviation/filt_std * filtered_offsets;
        
        % Save results
        offsets_from_reference(:,ith_trajectory) = filtered_offsets;
        
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


%% Perform iterations over trajectories

for ith_trajectory =1:num_trajectories
    % Show user what we are doing?
    if flag_do_debug
        fprintf(1,'Generating random trajectory: %.0d / %.0d \n',ith_trajectory,num_trajectories);
    end
    
    
    %% For each traversal, project from reference orthogonally
    % Set the projection type (see help in function below for details)
    flag_rounding_type = 4; % This averages the projection vectors along segments
    
    % Find the unit normal vectors at each of the station points
    [unit_normal_vector_start, unit_normal_vector_end] = ...
        fcn_Path_findOrthogonalTraversalVectorsAtStations(...
        reference_station_points(:,ith_trajectory),reference_traversal,flag_rounding_type);
    unit_vectors = unit_normal_vector_end - unit_normal_vector_start;
    
    %% Calculate random path and traversal and save into final structure
    random_path = unit_normal_vector_start + unit_vectors.*offsets_from_reference(:,ith_trajectory);
    random_traversal = fcn_Path_convertPathToTraversalStructure(random_path);
    random_traversals.traversal{ith_trajectory} = random_traversal;
    
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
    clf;
    hold on;
    
    % Plot the reference trajectory first
    plot(reference_traversal.X,reference_traversal.Y,'b.-','Linewidth',4,'Markersize',20);
    
    % Plot the random results
    fcn_Path_plotTraversalsXY(random_traversals,fig_num);
    title('Reference traversal and random traversals');
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
