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
%      fig_num: a figure number to plot results.
%
%
% OUTPUTS:
%
%      (none)
%
% DEPENDENCIES:
%
%      fcn_Path_checkInputsToFunctions
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
%     2021_01_05:
%     -- wrote the code originally
%     2021_01_08:
%     -- added input checking
%     2022_01_03:
%     -- corrected dependency list of functions
%     -- broke out fcn_Path_plotTraversalXYWithUpperLowerBands into
%     separate functionality (it is very useful for other functions!)

flag_do_debug = 0; % Flag to show the results for debugging
flag_do_plots = 1; % % Flag to plot the final results
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'Starting function: %s, in file: %s\n',st(1).name,st(1).file);
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
    if nargin < 1 || nargin > 3
        error('Incorrect number of input arguments')
    end
    
    % Check the reference_traversal input
    fcn_Path_checkInputsToFunctions(reference_traversal, 'traversal');
    
    
end

% Grab key variables
Station_reference = reference_traversal.Station;
Nstations = length(Station_reference(:,1));


%% Set defaults
% the default standard deviation
std_deviation = fcn_Path_calcSingleTraversalStandardDeviation(reference_traversal);

% the default number of points to use
num_points = Nstations;

%% Check for variable argument inputs (varargin)

% Does the user want to specify standard deviation?
if 2 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        std_deviation = temp;
    end
end


% Does user want to specify the figure to show the plots?
if 3 == nargin
    fig_num = varargin{2};
    figure(fig_num);
else
    if flag_do_plots
        fig = figure;
        fig_num = fig.Number;
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

%% plot the final XY result
fcn_Path_plotTraversalXYWithUpperLowerBands( reference_traversal, upper_traversal, lower_traversal, fig_num);


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
if flag_do_debug
    % Nothing to put in here!   
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends main function




