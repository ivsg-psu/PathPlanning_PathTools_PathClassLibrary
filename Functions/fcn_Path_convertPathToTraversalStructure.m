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
%      fig_num: a figure number to plot results.
%
% DEPENDENCIES:
%
%      fcn_Path_checkInputsToFunctions
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
%     2020_11_12 
%     -- wrote the code
%     2021_01_06
%     -- added functions for input checking
%     -- added figure variable input option for yaw plotting
%     -- cleaned up unnecessary yaw calculation (which was wrong!)
%     2021_01_07
%     -- changed input to path type, not X, Y separately
%     -- cleaned up comments
%     -- fixed error in comment about Yaw being front-padded. It is an
%     (N-1) x 1 vector now, not N x 1. 
%     2021_03_20
%     -- changed input checking to allow 3D path types
%     -- fixed possible bug in diff calculation

flag_do_debug = 0; % Flag to show the results for debugging
flag_do_plots = 0; % Flag to plot the final results
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

if flag_check_inputs == 1
    % Are there the right number of inputs?
    if nargin < 1 || nargin > 2
        error('Incorrect number of input arguments')
    end
    
    % Check the path input
    fcn_Path_checkInputsToFunctions(path, 'path2or3D');

end

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
    fcn_Path_plotTraversalsXY(data,fig_num);
    xlabel('X [m]');
    ylabel('Y [m]');
    
    
    % Plot the yaw plot
    figure(fig_num+100);
    hold on;
    plot(traversal.Station(1:end-1,1),traversal.Yaw*180/pi);
    xlabel('Station [m]');
    ylabel('Yaw angle [deg]');
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end
end
