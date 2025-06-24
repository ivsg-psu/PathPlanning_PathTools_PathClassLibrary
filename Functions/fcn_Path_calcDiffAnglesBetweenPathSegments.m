function diff_angles = fcn_Path_calcDiffAnglesBetweenPathSegments(Path,varargin)
% fcn_Path_calcDiffAnglesBetweenSegments
% Calculates the change in angles between path segments. If there are N
% points in the path, there are N-1 segments and thus N-2 angles between
% segments.
%
% Note that this method uses the dot product and the cross product to avoid
% errors caused by angle rollover that occur with arctan calculations. The
% arctan calculation gives the wrong answer with path segments are pointed
% near or at the -180 degree crossover point.
%
% FORMAT: 
%
%       diff_angles = fcn_Path_calcDiffAnglesBetweenPathSegments(Path,(fig_num))
%
% INPUTS:
%
%      Path: an N x 2 vector with [X Y] data in each row. N must be >= 3.
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
%      diff_angles: tan (N-2) x 1 vector of the change in angles in radians
%
% DEPENDENCIES
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%      
%       See the script:
%       script_test_fcn_Path_calcDiffAnglesBetweenPathSegments.m for a full
%       test suite.
%
% This function was written on 2021_01_03 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2021_01_03
% -- first writing of the code
% 2021_01_06
% -- fixed typos in the comments
% 2021_01_07
% -- fixed typos in the comments, minor header clean-ups
% 2025_06_23 - S. Brennan
% -- Updated debugging and input checks

% TO-DO
% (none)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==2 && isequal(varargin{end},-1))
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
        narginchk(1,2);

        % Check the Path variables
        fcn_DebugTools_checkInputsToFunctions(Path, 'paths');
    end
end

% Does user want to show the plots?
flag_do_plots = 0; % Default is to NOT show plots
if (0==flag_max_speed) && (2 == nargin) 
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

% Note that this method uses the dot product and the cross product to avoid
% errors caused by angle rollover. One could use the arctan method of
% calculating angles, which works in general but fails for paths that point
% straight to the left due to the cross-over point for the atan2
% calculation. For example, the following code does NOT always work:
%
%     % angles = atan2(path_average(2:end,2)-path_average(1:end-1,2),path_average(2:end,1)-path_average(1:end-1,1));
%     % angles = [angles; angles(end)];  % Pad the last point twice
%     % diff_angles2 = abs(diff(angles));

a_vector = Path(2:end-1,:)-Path(1:end-2,:);
b_vector = Path(3:end,:)-Path(2:end-1,:);
a_dot_b = sum(a_vector.*b_vector,2); % Do the dot product
a_cross_b = crossProduct(a_vector,b_vector);

a_mag = sum(a_vector.^2,2).^0.5;
b_mag = sum(b_vector.^2,2).^0.5;

% Calculate the change in angle using the dot product
%diff_angles = [acos(a_dot_b./(a_mag.*b_mag)); 0];
diff_angles = sign(a_cross_b).*(acos(a_dot_b./(a_mag.*b_mag)));


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
    clf;
    hold on;
    grid on;
    xlabel('Index');
    ylabel('Angle [deg]');

    % Plot the angle differences
    plot(diff_angles*180/pi,'k.-','Linewidth',3,'Markersize',25);

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
%% Calculate cross products
function result = crossProduct(v,w)
result = v(:,1).*w(:,2)-v(:,2).*w(:,1);
end

