function offset_traversals = fcn_Path_fillOffsetTraversalsAboutTraversal(reference_traversal, offsets, varargin)
% fcn_Path_fillOffsetTraversalsAboutTraversal
% fills in an array of traversals about a reference traversal at
% user-defined offset distances.
%
% FORMAT:
%
%      offset_traversals = ...
%      fcn_Path_fillOffsetTraversalsAboutTraversal(...
%            reference_traversal,...
%            offsets,...
%            (flag_rounding_type),
%            (fig_num));
%
% INPUTS:
%
%      reference_traversal: the traversal that is being used for randomly
%      generating all the other traversals
%
%      offsets: a scalar value or Mx1 array representing the offsets to be
%      used in calculating the resulting traversal producing one traversal
%      per offset value. Positive values follow the cross-product notation,
%      namely if a traversal is running from left-to-right, then a positive
%      offset would be above it, negative would be below.
%
%      (OPTIONAL INPUTS)
% 
%      flag_rounding_type: determines type of projection, and is passed
%      into fcn_Path_findOrthogonalTraversalVectorsAtStations. See that
%      function for more explanation. Default (empty) is used unless
%      changed via this input.
%
%      fig_num: a figure number to plot results.
%
%
% OUTPUTS:
%
%      offset_traversals: a structure containing the resulting traversals
%      that are generated, one for each given offset
%
% DEPENDENCIES:
%
%      fcn_Path_checkInputsToFunctions
%      fcn_Path_calcSingleTraversalStandardDeviation
%      fcn_Path_findOrthogonalTraversalVectorsAtStations
%      fcn_Path_convertPathToTraversalStructure
%      fcn_Path_plotTraversalsXY
%
% EXAMPLES:
%
%     See the script: script_test_fcn_Path_fillOffsetTraversalsAboutTraversal
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


% Check inputs?
if flag_check_inputs
    % Are there the right number of inputs?
    narginchk(2,4)
        
    % Check the reference_traversal input
    fcn_Path_checkInputsToFunctions(reference_traversal, 'traversal');
            
    % Check the offsets input (looks like a station type)
    fcn_Path_checkInputsToFunctions(offsets, 'station');
    
end

Station_reference = reference_traversal.Station;
Nstations = length(Station_reference);


%% Set defaults

% the default number of trajectories to use
num_trajectories = length(offsets(:,1));

% the default number of points to use
num_points = Nstations;

% Does user want to specify flag_rounding_type?
flag_rounding_type = []; % Default is to not do any plotting
if 3 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        flag_rounding_type = temp;
    end
end

%% Check for variable argument inputs (varargin)

% Does user want to show the plots?
if 4 == nargin
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        fig_num = temp;
        flag_do_plots = 1;
    end
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
    fcn_Path_findOrthogonalTraversalVectorsAtStations(...
    reference_station_points(:,1),reference_traversal,flag_rounding_type);
unit_vectors = unit_normal_vector_end - unit_normal_vector_start;

%% Perform iterations over trajectories

for ith_trajectory =1:num_trajectories
    % Show user what we are doing?
    if flag_do_debug
        fprintf(1,'Generating offset trajectory: %.0d / %.0d \n',ith_trajectory,num_trajectories);
    end

    % Calculate random path and traversal and save into final structure
    offset_path = unit_normal_vector_start + unit_vectors.*offsets_from_reference(:,ith_trajectory);
    offset_traversal = fcn_Path_convertPathToTraversalStructure(offset_path);
    offset_traversals.traversal{ith_trajectory} = offset_traversal;
    
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
    plot(reference_traversal.X,reference_traversal.Y,'b.-','Linewidth',4,'Markersize',20);
    
    % Plot the ofset results
    fcn_Path_plotTraversalsXY(offset_traversals,fig_num);
    title('Reference traversal and offset traversals');
    xlabel('X [m]');
    ylabel('Y [m]'); 
    
    % Add a legend
    legend('Reference traversal', 'Offset traversals');
    
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends main function

