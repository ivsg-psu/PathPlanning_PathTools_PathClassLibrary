function yaw_angles_in_radians = fcn_Path_calcYawFromPathSegments(Path,varargin)
% fcn_Path_calcYawFromPathSegments
% Calculates the yaw angles of the path segments. If there are N points
% in the path, there are N-1 segments and thus N-1 angles. 
%
% The method used for calculating the path angles is to use the first
% segment of the path to calculate the angle of this segment, using the
% arctangent function (atan2). For any subsequent segments, if any, the
% change in angles is calculated using the dot product and cross product,
% and these angle changes are then cumulatively added to determine yaw
% angles thereafter. This method with the dot product and the cross product
% allows us to avoid errors caused by angle rollover that occurs with the
% arctangent functions. Thus, a path that circles clockwise two times will
% have a yaw angle that starts somewhere between -pi and pi, and decreases
% by abou 4 pi, or two revolutions. Note that the angle convention for yaw
% follows Cartesian definitions, namely yaw is positive when angles
% increase counter-clockwise, and yaw angles start at 0 at the x-axis.
%
% FORMAT: 
%
%     yaw_angles_in_radians = ...
%     fcn_Path_calcYawFromPathSegments(Path, (fig_num))
%
% INPUTS:
%
%     Path: an N x 2 vector with [X Y] data in each row. N must be >= 2.
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
%      yaw_angles_in_radians: the N-1 yaw angles in radians
%
% DEPENDENCIES
%
%      fcn_Path_calcDiffAnglesBetweenPathSegments
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%      
%       See the script:
%       script_test_fcn_Path_calcYawFromPathSegments.m for a full
%       test suite.
%
% This function was written on 2021_01_06 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2021_01_06
% -- wrote the code
% 2021_01_06
% -- minor comment clean-ups
% 2025_06_23 - S. Brennan
% -- Updated debugging and input checks

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

        % Check the Path variables
        fcn_DebugTools_checkInputsToFunctions(Path, 'path');
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
else
    if flag_do_debug
        fig = figure;  
        fig_num = fig.Number;
        flag_do_plots = 1;
    end
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

starting_angle = atan2(Path(2,2)-Path(1,2),Path(2,1)-Path(1,1));

Nangles = length(Path(:,1)) - 1;
if Nangles>1
    diff_angles = fcn_Path_calcDiffAnglesBetweenPathSegments(Path);
    cumulative_angles = [0; cumsum(diff_angles)];
else
    cumulative_angles = 0;
end
yaw_angles_in_radians = cumulative_angles + starting_angle;


%% Any debugging?
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
    xlabel('Index');
    ylabel('Yaw Angle [deg]');

    % Plot the angle differences
    plot(yaw_angles_in_radians*180/pi,'.-','Linewidth',3,'Markersize',25);

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end
end


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

