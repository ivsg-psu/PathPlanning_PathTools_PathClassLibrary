function traversal = fcn_Path_convertPathToTraversalStructure(path,varargin)
% fcn_Path_convertPathToTraversalStructure
% Takes a Path type and creates it to a traversal structure
%
% FORMAT: 
%
%       traversal = fcn_Path_convertPathToTraversalStructure(path,(fig_num))
%
% INPUTS:
%
%      path: an N x 2 vector of [X Y] positions, with N>=2 OR
%            an N x 3 vector of [X Y Z] positions, with N>=2
%
% OUTPUTS:
%
%      traversal: a sttructure containing the following fields
%            - X: an N x 1 vector that is a duplicate of the input X
%            - Y: an N x 1 vector that is a duplicate of the input Y
%            - Z: an N x 1 vector that is a duplicate of the input Z (if 3D) OR
%                 an N x 1 vector that is a zero array the same length as X
%                 (if 2D)
%            - Diff: a [N x 2] or [N x 3] array that is the change in X and Y
%            (front-padded with [0 0 (0)])
%            - Station: the XYZ distance as an N x 1 vector, representing
%            the distance traveled up to the current point (starting with 0
%            at the first point)
%            - Yaw: the calculated yaw angle (radians) of each path
%            segment, where yaw calculations are in the XY plane
%            (note: there are N-1 segments if there are N points, thus
%            there are N-1 Yaw points)
%
%      (OPTIONAL INPUTS)
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_Path_calcDiffAnglesBetweenPathSegments
%      fcn_Path_calcYawFromPathSegments
%      fcn_Path_plotTraversalsXY
%
% EXAMPLES:
%      
%       See the script:
%       script_test_fcn_Path_convertPathToTraversalStructure.m for a full
%       test suite. 
%
% This function was written on 2020_11_12 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2020_11_12
% -- wrote the code
% 2021_01_06
% -- added functions for input checking
% -- added figure variable input option for yaw plotting
% -- cleaned up unnecessary yaw calculation (which was wrong!)
% 2021_01_07
% -- changed input to path type, not X, Y separately
% -- cleaned up comments
% -- fixed error in comment about Yaw being front-padded. It is an
% (N-1) x 1 vector now, not N x 1.
% 2021_03_20
% -- changed input checking to allow 3D path types
% -- fixed possible bug in diff calculation
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
        fcn_DebugTools_checkInputsToFunctions(path, 'path2or3D');
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


%% Solve for the circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fill in the fields typically required on a traversal
traversal.X = path(:,1);
traversal.Y = path(:,2);

% Calculate the station differences, and force the diff operation to occur
% along rows, not columns
if length(path(1,:))==3
    traversal.Z = path(:,3);
    traversal.Diff = [[0 0 0]; diff(path,1,1)];
else
    traversal.Z = 0*path(:,1);
    traversal.Diff = [[0 0]; diff(path,1,1)];
end

traversal.Station = cumsum(sqrt(sum(traversal.Diff.^2,2)));
traversal.Yaw = ...
    fcn_Path_calcYawFromPathSegments(path(:,1:2)); % NOTE: Yaw should be (N-1) x 1


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
    % plot the final XY result
    figure(fig_num); 
    data.traversal{1} = traversal;

    % Plot the XY plot
    subplot(2,1,1);    
    fcn_Path_plotTraversalsXY(data,fig_num);
    hold on;
    xlabel('X [m]');
    ylabel('Y [m]');
    
    
    % Plot the yaw plot
    subplot(2,1,2);
    plot(traversal.Station(1:end-1,1),traversal.Yaw*180/pi);
    hold on;
    xlabel('Station [m]');
    ylabel('Yaw angle [deg]');
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