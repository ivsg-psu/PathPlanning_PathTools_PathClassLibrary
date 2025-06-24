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
%      (none)
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
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
% 2022_01_03:
% -- wrote the code originally, using fcn_Path_plotTraversalXYWithVarianceBands
% 2025_06_23 - S. Brennan
% -- Updated debugging and input checks

% TO-DO
% (none)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==4 && isequal(varargin{end},-1))
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
        narginchk(3,4);

        % Check the middle_traversal input
        fcn_DebugTools_checkInputsToFunctions(middle_traversal, 'traversal');

        % Check the middle_traversal input
        fcn_DebugTools_checkInputsToFunctions(upper_traversal, 'traversal');

        % Check the middle_traversal input
        fcn_DebugTools_checkInputsToFunctions(lower_traversal, 'traversal');

    end
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


% Does user want to show the plots?
flag_do_plots = 1; % Default is to make a plot
fig_num = [];
if (0==flag_max_speed) && (4 == nargin) 
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        fig_num = temp;
        figure(fig_num);
        flag_do_plots = 1;
    end
else
    if flag_do_debug
        fig = figure;  
        fig_num = fig.Number;
        flag_do_plots = 1;
    end
end
if isempty(fig_num)
    temp = figure;
    fig_num = temp.Number;
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

% Generate top and bottom paths
top_path = [X_upper, Y_upper];
bottom_path = [X_lower, Y_lower];


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
