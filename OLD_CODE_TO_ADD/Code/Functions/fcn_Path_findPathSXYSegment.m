
function [pathSXY_segment,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findPathSXYSegment(...
    pathSXY, s_coord_start,s_coord_end, varargin)
% fcn_Path_findPathSXYSegment
% Finds portion of a SXY-type path that contains the given s_coordinates, starting
% from s_coord_start to s_coord_end
% 
% Format: 
% [closest_path_point,s_coordinate] = fcn_Path_findPathSXYSegment(point, path,varargin)
%
% INPUTS:
%
%      path: a Nx3 vector of [S X Y] path points, where N is the number of
%      points the points on the path, N >= 2. S is the station distance,
%      and X and Y are the XY points of the path coordinates
%
%      s_coord_start: a 1x1 (scalar) indicating the s-coordinate location
%      at which the query starts. The path segment output will start at
%      previous s-value to this station.
%
%      s_coord_end: a 1x1 (scalar) indicating the s-coordinate location
%      at which the query ends. The path segment output will end at
%      subsequent s-value to this station.
%
%      If the query includes s-values outside the path, then the flag_
%      variables are set to 1. Note that this function always returns at
%      least 2 points representing the closest path segment, even if both
%      s-point queries are outside the given path.
%
%      (optional_input) figure_number: plots the results into the given
%      figure
%
% OUTPUTS:
%
%      path_segment: a Mx3 vector of [S X Y] path points, where M is the
%      number of points the points on the path segment, M >= 2. S is the
%      station distance, and X and Y are the XY points of the path
%      coordinates
%
%      flag_outside_start, flag_outside_end: flags that are set equal to 1
%      if the query is outside the s-distance within the given path at
%      either the start, the end, or both
%
% EXAMPLES:
%      
% See the script: 
% SCRIPT_test_fcn_Path_findPathSXYSegment.m
% for a full test suite.
%
% This function was written on 2020_10_14 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%     2020_10_14:
%     - first write of the code
%     2020_11_15:
%     - changed the name to prep for Paths class

flag_do_debug = 0; % Flag to plot the results for debugging
flag_check_inputs = 1; % Flag to perform input checking

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
% Are the input vectors the right shape?
Npoints_in_path = length(pathSXY(:,1));

if flag_check_inputs == 1
    % Are there the right number of inputs?
    if nargin < 3 || nargin > 4
        error('Incorrect number of input arguments')
    end
    
    if Npoints_in_path<2
        error('The path vector must have at least 2 rows, with each row representing a different (x y) point');
    end
    if length(pathSXY(1,:))~=3
        error('The path vector must have 3 columns, with column 1 representing the s-distance, column 2 representing the x portions of the points, column 3 representing the y portions.');
    end
    
    if s_coord_start>s_coord_end
        warning('S coordinates of start and end seem out of order. These will be automatically sorted but the results may be incorrect.');
        s_coord_start_new = s_coord_end;
        s_coord_end = s_coord_start;
        s_coord_start = s_coord_start_new;
    end
end

% Does user have special variable inputs?
if 4 == nargin
    fig_num = varargin{1};
    figure(fig_num);
    flag_do_debug = 1;
else
    if flag_do_debug
        fig = figure;  %#ok<UNRCH>
        fig_num = fig.Number;
    end
end

%% Find the closest point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set default outputs
flag_outside_start = 0;
flag_outside_end   = 0;

% Use find function to grab values
path_segment_start_index = find(pathSXY(:,1)<s_coord_start,1,'last');
path_segment_end_index   = find(pathSXY(:,1)>s_coord_end,1,'first');

% Check if we've gone past the ends of the path
if isempty(path_segment_start_index) % There is no s-coordinate in the path smaller than the start
    path_segment_start_index = 1;
    flag_outside_start = 1;
end
if isempty(path_segment_end_index) % There is no s-coordinate in the path larger than the end
    path_segment_end_index = Npoints_in_path;
    flag_outside_end   = 1;
end

% Check to see if start index is the end of path
if Npoints_in_path == path_segment_start_index  % All path s-coordinates are smaller than the start
    path_segment_start_index = Npoints_in_path - 1;
    path_segment_end_index = Npoints_in_path;
end

% Check to see if end index is the start of path
if 1 == path_segment_end_index  % All path s-coordinates are larger than the end
    path_segment_start_index = 1;
    path_segment_end_index = 2;
end

% Check to see if start and end of segment are the same (degenerate case)
if path_segment_start_index == path_segment_end_index
    path_segment_start_index = path_segment_end_index - 1;
end

% Grab the path segment that is closest to the robot's look-ahead point
pathSXY_segment = pathSXY(path_segment_start_index:path_segment_end_index,:);


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
    figure(fig_num);
    hold on;
    grid on;
    % Plot the path
    plot(pathSXY(:,2),pathSXY(:,3),'r-','Linewidth',5);       
    plot(pathSXY(:,2),pathSXY(:,3),'ro','Linewidth',5);       
    
    axis equal;
    
    % Plot the results
    plot(pathSXY_segment(:,2),pathSXY_segment(:,3),'b-','Linewidth',3);
    plot(pathSXY_segment(1,2),pathSXY_segment(1,3),'b*');   
    plot(pathSXY_segment(end,2),pathSXY_segment(end,3),'b*');   
    text(pathSXY_segment(1,2),pathSXY_segment(1,3),'Start');   
    text(pathSXY_segment(end,2),pathSXY_segment(end,3),'End');     
    
end % Ends the flag_do_debug if statement



end % Ends the function
