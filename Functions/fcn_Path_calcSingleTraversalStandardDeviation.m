function std_deviation = fcn_Path_calcSingleTraversalStandardDeviation(reference_traversal, varargin)
% fcn_Path_calcSingleTraversalStandardDeviation
% calculates the standard deviation in the offsets of a single traversal by
% analyzing the variance in angles along a reference_traversal, then
% multiplying these by the average segment length in the reference
% traversal. The resulting standard deviation approximates the
% variance in lateral offset that occurs at the end of each segment, versus
% a line projected from the previous segment.
%
% FORMAT:
%
%      std_deviation = ...
%      fcn_Path_calcSingleTraversalStandardDeviation(...
%            reference_traversal,...
%            (fig_num));
%
% INPUTS:
%
%      reference_traversal: the traversal that is being used for randomly
%      generating all the other traversals
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results.
%
%
% OUTPUTS:
%
%      std_deviation: the standard deviation in the offset of the traversal
%
% DEPENDENCIES:
%
%      fcn_Path_checkInputsToFunctions
%      fcn_Path_calcDiffAnglesBetweenPathSegments
%      fcn_Path_plotTraversalXYWithVarianceBands
%
% EXAMPLES:
%
%     See the script: script_test_fcn_Path_calcSingleTraversalStandardDeviation
%     for a full test suite.
%
% This function was written on 2021_01_05 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
%     2021_01_05:
%     -- wrote the code originally
%     2021_01_06:
%     -- added functions for input checking
%     2021_01_07:
%     -- fixed typos in comments, and in header


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
    if nargin < 1 || nargin > 2
        error('Incorrect number of input arguments')
    end
    
    % Check the reference_traversal input
    fcn_Path_checkInputsToFunctions(reference_traversal, 'traversal');
    
end


%% Check for variable argument inputs (varargin)

% Does user want to show the plots?
if 2 == nargin
    fig_num = varargin{1};
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

% Standard deviation: this is calculated by the standard deviation in path
% angles times the mean segment length. The reason for this is that this
% distance represents the average deviation laterally, right or left, of
% each segment. This is a good first guess to produce paths of similar
% curviness to the original path.

% FIll in useful variables
X_reference = reference_traversal.X;
Y_reference = reference_traversal.Y;
Station_reference = reference_traversal.Station;
reference_path = [X_reference Y_reference];

% Calculate angle changes between points
diff_angles = fcn_Path_calcDiffAnglesBetweenPathSegments(reference_path);
std_angles = std(diff_angles);
mean_segment_length = mean(diff(Station_reference));
std_deviation = std_angles*mean_segment_length;


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
       
    % Plot the results
    fcn_Path_plotTraversalXYWithVarianceBands(reference_traversal,...
    std_deviation,fig_num);
    title(sprintf('Standard deviation found to be: %.2f',std_deviation));
    xlabel('X [m]');
    ylabel('Y [m]'); 
    
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends main function

