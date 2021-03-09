function [traversal_trimmed,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    long_traversal, s_coord_start,s_coord_end, varargin)
% fcn_Path_findTraversalStationSegment
% Finds portion of a long_traversal that contains the given s_coordinates,
% starting from s_coord_start to s_coord_end, and returns that portion as
% another traversal "trimmed" out of the original by finding the path
% segments within the long traversal closest to the queried s-coordinates.
% If the s-coordinates are outside those of the long_traversal, then flags
% are set to 1 for either flag_outside_start, flag_outside_end; otherwise,
% these flags are zero.
%
% Note that this function always returns at least 2 points representing the
% closest path segment, even if both s-point queries are outside the given
% path.
%
% 
% FORMAT: 
%
%      [traversal_trimmed,flag_outside_start, flag_outside_end] = ...
%      fcn_Path_findTraversalStationSegment(...
%      long_traversal, s_coord_start,s_coord_end, 
%      (fig_num))
%
% INPUTS:
%
%      long_traversal:  the traversal that is being used for trimming by
%      the given station coordinates.
%
%      s_coord_start: a 1x1 (scalar) indicating the s-coordinate location
%      at which the query starts. The path segment output will start at
%      previous s-value to this station.
%
%      s_coord_end: a 1x1 (scalar) indicating the s-coordinate location
%      at which the query ends. The path segment output will end at
%      subsequent s-value to this station.
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      traversal_trimmed: a traversal output trimmed out of the original
%      traversal. 
%
%      flag_outside_start, flag_outside_end: flags that are set equal to 1
%      if the query is outside the s-distance within the given path at
%      either the start, the end, or both
%
%
% DEPENDENCIES:
%
%      fcn_Path_checkInputsToFunctions
%      fcn_Path_convertPathToTraversalStructure
%
% EXAMPLES:
%      
% See the script: 
% script_test_fcn_Path_findTraversalStationSegment.m
% for a full test suite.
%
% This function was written on 2020_10_14 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%     2020_10_14:
%     - first write of the code
%     2020_11_15:
%     - changed the name to prep for Paths class
%     2021_01_08
%     -- started updating for new class
%     2021_01_09
%     -- updated name and types to take traversal inputs
%     -- added input checking
%     -- added flag_do_plots    

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
% Are the input vectors the right shape?
Npoints_in_path = length(long_traversal.Station);

if flag_check_inputs == 1
    % Are there the right number of inputs?
    if nargin < 3 || nargin > 4
        error('Incorrect number of input arguments')
    end
    
    % Check the long_traversal input
    fcn_Path_checkInputsToFunctions(long_traversal, 'traversal');
    
    % Check the s_coord_start input
    fcn_Path_checkInputsToFunctions(s_coord_start, 'station');

    % Check the s_coord_end input
    fcn_Path_checkInputsToFunctions(s_coord_end, 'station');
   
    % Check that the start is strictly before the end
    if s_coord_start>s_coord_end
        warning('S coordinates of start and end seem out of order. These will be automatically sorted but the results may be incorrect.');
        s_coord_start_new = s_coord_end;
        s_coord_end = s_coord_start;
        s_coord_start = s_coord_start_new;
    end
end

% Does user want to show the plots?
if 4 == nargin
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
path_segment_start_index = find(long_traversal.Station < s_coord_start,1,'last');
path_segment_end_index   = find(long_traversal.Station >s_coord_end,1,'first');

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
% Are all path s-coordinates are smaller than the start?
if Npoints_in_path == path_segment_start_index  
    path_segment_start_index = Npoints_in_path - 1;
    path_segment_end_index = Npoints_in_path;
end

% Check to see if end index is the start of path?
if 1 == path_segment_end_index  % All path s-coordinates are larger than the end
    path_segment_start_index = 1;
    path_segment_end_index = 2;
end

% Check to see if start and end of segment are the same (degenerate case)
if path_segment_start_index == path_segment_end_index
    path_segment_start_index = path_segment_end_index - 1;
end

% Grab the path segment that is closest to the given s-coordinates
traversal_trimmed_path = ...
    [long_traversal.X(path_segment_start_index:path_segment_end_index,:),...
    long_traversal.Y(path_segment_start_index:path_segment_end_index,:)];

traversal_trimmed = fcn_Path_convertPathToTraversalStructure(...
    traversal_trimmed_path);

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
    % Plot the path
    plot(long_traversal.X,long_traversal.Y,'r-','Linewidth',5);       
    plot(long_traversal.X,long_traversal.Y,'ro','Linewidth',5);       
    
    axis equal;
    
    % Plot the results
    plot(traversal_trimmed.X,traversal_trimmed.Y,'b-','Linewidth',3);
    plot(traversal_trimmed.X(1),traversal_trimmed.Y(1),'b*');   
    plot(traversal_trimmed.X(end),traversal_trimmed.Y(end),'b*');   
    text(traversal_trimmed.X(1),traversal_trimmed.Y(1),'Start');   
    text(traversal_trimmed.X(end),traversal_trimmed.Y(end),'End');     
    
end % Ends the flag_do_debug if statement


if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends the function
