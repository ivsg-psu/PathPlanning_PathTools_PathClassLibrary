function fcn_Path_plotTraversalXYWithUpperLowerBands(middle_traversal, upper_traversal, lower_traversal, varargin)
% fcn_Path_plotTraversalXYWithUpperLowerBands
% Plots a traversal with a band defined by an upper and lower traversal.
% All traversals must have the same data length.
%
% FORMAT:
%
%      fcn_Path_plotTraversalXYWithUpperLowerBands(...
%            middle_traversal,...
%            upper_traversal,...
%            lower_traversal,...
%            (fig_num));
%
% INPUTS:
%
%      middle_traversal: the traversal that is being used for the middle
%      plot
%
%      upper_traversal: the traversal that is being used to define the
%      upper band
%
%      lower_traversal: the traversal that is being used to define the
%      lower band
%
%      (OPTIONAL INPUTS)
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
%      fcn_Path_calcSingleTraversalStandardDeviation
%      fcn_Path_findOrthogonalTraversalVectorsAtStations
%
% EXAMPLES:
%
%     See the script: script_test_fcn_Path_plotTraversalXYWithUpperLowerBands
%     for a full test suite.
%
% This function was written on 2022_01_03 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
%     2022_01_03:
%     -- wrote the code originally, using fcn_Path_plotTraversalXYWithVarianceBands

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

%      fcn_Path_plotTraversalXYWithUpperLowerBands(...
%            middle_traversal,...
%            upper_traversal,...
%            lower_traversal,...
%            (fig_num));


% Check inputs?
if flag_check_inputs
    % Are there the right number of inputs?
    if nargin < 3 || nargin > 4
        error('Incorrect number of input arguments')
    end
    
    % Check the middle_traversal input
    fcn_Path_checkInputsToFunctions(middle_traversal, 'traversal');
    
    % Check the middle_traversal input
    fcn_Path_checkInputsToFunctions(upper_traversal, 'traversal');
    
    % Check the middle_traversal input
    fcn_Path_checkInputsToFunctions(lower_traversal, 'traversal');
    
    
end

% Grab key variables
X_middle = middle_traversal.X;
Y_middle = middle_traversal.Y;
Nstations = length(X_middle(:,1));

X_upper = upper_traversal.X;
Y_upper = upper_traversal.Y;
if Nstations~=length(X_upper(:,1))
    error('The number of data points in the upper_traversal must match the middle_traversal');
end

X_lower = lower_traversal.X;
Y_lower = lower_traversal.Y;
if Nstations~=length(X_lower(:,1))
    error('The number of data points in the lower_traversal must match the middle_traversal');
end

% Does user want to specify the figure to show the plots?
if 4 == nargin
    fig_num = varargin{1};
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


%% Generate top and bottom paths
top_path = [X_upper, Y_upper];
bottom_path = [X_lower, Y_lower];


% plot the final XY result
figure(fig_num);

% Check to see if the hold was on?
flag_hold_was_off = 0;
if ~ishold
    flag_hold_was_off = 1;
    hold on;
end

% Plot the reference trajectory first
main_plot_handle = plot(X_middle,Y_middle,'.-','Linewidth',4,'Markersize',20);
plot_color = get(main_plot_handle,'Color');

% % Now make the patch as one object (THIS ONLY WORKS IF NO CROSSINGS)
% x_vector = [top_path(:,1)', fliplr(bottom_path(:,1)')];
% y_vector = [top_path(:,2)', fliplr(bottom_path(:,2)')];
% patch = fill(x_vector, y_vector,[128 193 219]./255);
% set(patch, 'edgecolor', 'none');
% set(patch, 'FaceAlpha', 0.5);

% Now make the patch segment by segment
for i_patch = 2:Nstations
    x_vector = [top_path((i_patch-1):i_patch,1)', fliplr(bottom_path((i_patch-1):i_patch,1)')];
    y_vector = [top_path((i_patch-1):i_patch,2)', fliplr(bottom_path((i_patch-1):i_patch,2)')];
    patch = fill(x_vector, y_vector,plot_color);
    %patch = fill(x_vector, y_vector,(plot_color*0.8 + 0.2*[1 1 1]));
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', 0.2);
end

% Put hold back to the original state
if flag_hold_was_off
    hold off;
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
if flag_do_debug
    % Nothing to put in here!   
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends main function




