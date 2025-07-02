function [path_no_pinch_point] = ...
    fcn_Path_removePinchPointInPath(...
    path_with_pinch_point,...
    varargin)
% fcn_Path_removePinchPointInPath
% Given a path with a pinch point - an area where the path
% suddenly bends back on itself before continuing - this function removes
% the pinch point
%
% FORMAT:
%
%     [path_no_pinch_point] = ...
%         fcn_Path_removePinchPointInPath(...
%         path_with_pinch_point
%        (fig_num));
%
% INPUTS:
%
%      path_with_pinch_point: a N x 2 or N x 3 set of coordinates
%      representing the [X Y] or [X Y Z] coordinates, in sequence, that
%      specifies the path with a pinch point
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
%      path_no_pinch_point: a path structure that specifies the
%      path, s-coordinates, etc of a traveral with the pinch points
%      removed.
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_Path_plotTraversalsXY
%
% EXAMPLES:
%
% See the script: script_test_fcn_Path_removePinchPointInPath
% for a full test suite.
%
% This function was written on 2021_01_23 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2021_01_23:
% -- first write of the code
% 2025_06_23 - S. Brennan
% -- Updated debugging and input checks
% 2025_07_01 - S. Brennan
% -- Removed traversal input type and replaced with cell array of paths
% -- Renamed function from fcn_Path_removePinchPointInTraversal

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

        % Check the path_with_pinch_point input
        fcn_DebugTools_checkInputsToFunctions(path_with_pinch_point, 'path2or3D');

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
    fig_debug = 12345; %#ok<NASGU>
end


%% Start of main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nsegments = length(path_with_pinch_point(:,1)) - 1;

if Nsegments<2
    % There needs to be at least 2 segments to self-intersect
    path_no_pinch_point = path_with_pinch_point;
else

    % Loop through all the segments in path 1, checking each for
    % intersections with path 2
    original_path = path_with_pinch_point;

    % Set the flag to indicate that we are ONLY searching if the exact segment
    % crosses or not
    flag_search_type = 0;

    % Initialize values
    no_pinch_path = original_path(1,:);
    remaining_path = original_path(2:end,:);

    % Are there at least 3 points in the remaining path?
    while length(remaining_path(:,1))>=2
        % Define the sensor vector
        sensor_vector_start = no_pinch_path(end,:);
        sensor_vector_end = remaining_path(1,:);
        sensor_vector_length = sum((sensor_vector_end - sensor_vector_start).^2,2).^0.5;

        % Check to see if there are intersections
        [distance,hit_location,path_segments] = ...
            fcn_Path_findProjectionHitOntoPath(...
            remaining_path,...
            sensor_vector_start,sensor_vector_end,...
            flag_search_type, -1);

        % Did we hit anything? If so, save it and set a flag that a pinch
        % point was hit!
        if isnan(distance) || (distance==sensor_vector_length(1,1) && path_segments==1)
            % Nothing hit
            no_pinch_path = [no_pinch_path; sensor_vector_end];  %#ok<AGROW>
            remaining_path = remaining_path(2:end,:);

        else % Hit something

            % Calculate how much s-distance is being "cut" to see if we
            % need to warn the user:

            % First, interpolate the s-coordinate after the 1st hit
            travel = sum((sensor_vector_end - hit_location).^2,2).^0.5;

            % Second: find the s-coordinate for path after this hit
            [~,s_coordinate_after_hit,~,~,~] = ...
                fcn_Path_snapPointToPathViaVectors(...
                hit_location, remaining_path, [], -1);

            % Third: add these values
            s_coordinates_trimmed = (travel + s_coordinate_after_hit);

            % Do we need to warn the user?
            if 1==0
                if s_coordinates_trimmed > 10
                    warning('10 meters or more were trimmed. This is a large amount!');
                end
            end

            % Update the path
            remaining_path = [hit_location; remaining_path(path_segments+1:end,:)];

            % Make sure we didn't just repeat a point in the path (which
            % happens if the loop comes back onto itself)
            if isequal(remaining_path(1,:),remaining_path(2,:))
                remaining_path = remaining_path(2:end,:);
            end

        end % Ends check to see if distance is empty

    end % Ends while loop
    no_pinch_path = [no_pinch_path; remaining_path];

    % Clean up the path by removing repeats - these can occur when the path
    % loops back onto itself
    cleaned_path = no_pinch_path(1,:);
    for ith_row = 2:length(no_pinch_path)
        if ~isequal(no_pinch_path(ith_row,:),no_pinch_path(ith_row-1,:))
            cleaned_path = [cleaned_path; no_pinch_path(ith_row,:)]; %#ok<AGROW>
        end
    end % Ends the for loop over rows, to clean repeats

    path_no_pinch_point = cleaned_path;
end % Ends check to see if there are at least 3 segments




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
    figure(fig_num);
    hold on;
    grid on;

    % Plot the two paths
    clear data
    data{1} = path_with_pinch_point;
    data{2} = path_no_pinch_point;
    fcn_Path_plotPathsXY(data,fig_num);

    legend('Original path with pinch point', 'Traversal with no pinch');
end % Ends the flag_do_plot if statement

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends the function


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
