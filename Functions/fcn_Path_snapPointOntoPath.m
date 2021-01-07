function [closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = fcn_Path_snapPointOntoPath(point, path,varargin)
% fcn_Path_snapPointOntoPath
% Finds location on a path that is closest to a given point, e.g. snapping
% the point onto the path
% 
% FORMAT: 
%
%      [closest_path_point,s_coordinate] = fcn_Path_snapPointOntoPath(point, path,varargin)
%
% INPUTS:
%
%      point: a 1x2 vector containing the [X Y] location of the point
%      path: a Nx2 vector of [X Y] path points, where N is the number of points the points on the path, N >= 2. 
%      (optional) figure_number: number where results are plotted
%
% OUTPUTS:
%
%      closest_path_point: a 1x2 vector containing the [X Y] location of
%      the nearest point on the path
%      s_coordinate: a scalar (1x1) representing the s-coordinate distance
%      along the path
%
% EXAMPLES:
%      
%      % BASIC example
%      point = [0.5 0.2];
%      path = [0 0 0; 1 1 0; 2 2 0; 3 2 1];
%      fignum = 222;
%      [closest_path_point,s_coordinate] = ...
%      fcn_Path_snapPointOntoPath(point, path,fignum)
% 
% See the script: script_test_fcn_Path_snapPointOntoPath
% for a full test suite.
%
% This function was written on 2020_10_10 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%     2020_10_10 - first write of the code
%     2020_11_10 - changed function names in prep for DataClean class
%     2020_12_03 - updated some of the plotting/debug details to improve

% TO-DO:
% Allow multiple points, e.g.
%      point: a Nx2 vector where N is the number of points, but at least 1. 

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
Npoints = length(path(:,1));

if flag_check_inputs == 1
    % Are there the right number of inputs?
    if nargin < 2 || nargin > 3
        error('Incorrect number of input arguments')
    end
    
    if Npoints<2
        error('The path vector must have at least 2 rows, with each row representing a different (x y) point');
    end
    if length(point(1,:))~=2
        error('The point vector must have 2 columns, with column 1 representing the x portions of the points, column 2 representing the y portions.');
    end
end

% Does user want to show the plots?
if 3 == nargin
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
 
% The solution method is as follows:
%  1. Find the closest point on the path to the query point
%  2. Find next closest point to that one to the query point
%  3. Sort these to find which is first/second
%  4. Find percentage of travel using dot products


% Find distance from point to path
distances_point_to_path = sum((path(:,2:3)-point).^2,2).^0.5; % Calculate distances
[~,closest_path_point_index] = min(distances_point_to_path); % Grab closest point

% Find next closest ADJACENT point - be sure to check end cases as well
if 1 == closest_path_point_index
    next_closest_path_point_index = 2;
elseif Npoints == closest_path_point_index
    next_closest_path_point_index = Npoints-1;
elseif distances_point_to_path(closest_path_point_index+1)<distances_point_to_path(closest_path_point_index-1)
    next_closest_path_point_index = closest_path_point_index+1;
else
    next_closest_path_point_index = closest_path_point_index-1;
end
first_path_point_index = min(closest_path_point_index,next_closest_path_point_index);
second_path_point_index = max(closest_path_point_index,next_closest_path_point_index);

% Do the dot products - define the vectors first
% See: https://mathinsight.org/dot_product for explanation
% Basically, we are seeing what amount the point_vector points in the
% direction of the path_vector
path_vector = path(second_path_point_index,2:3)-path(first_path_point_index,2:3);
path_segment_length = sum(path_vector.^2,2).^0.5;
point_vector = point-path(first_path_point_index,2:3);
projection_distance = dot(path_vector,point_vector)/path_segment_length; % Do dot product
percent_along_length = projection_distance/path_segment_length;

% Calculate the outputs
closest_path_point = path(first_path_point_index,2:3) + path_vector*percent_along_length;
s_coordinate = path(first_path_point_index,1) + path_segment_length*percent_along_length;


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
    plot(path(:,2),path(:,3),'r-','Linewidth',5);       
    plot(path(:,2),path(:,3),'ro','Markersize',20);       
    
    axis equal;
    
    % Plot the query point
    plot(point(:,1),point(:,2),'ko');
    text(point(:,1),point(:,2),'Query point');
    
    % Plot the closest path points;
    plot(...
        path(first_path_point_index:second_path_point_index,2),...
        path(first_path_point_index:second_path_point_index,3),'r*');       
    
    % % Label the points with distances?
    % for i_point = 1:length(path(:,1))
    %     text(path(i_point,2),path(i_point,3),sprintf('%.2f',distances_point_to_path(i_point)));
    % end
    
    % Plot the closest point on path
    plot(closest_path_point(:,1),closest_path_point(:,2),'go','Markersize',20);
    text(closest_path_point(:,1),closest_path_point(:,2),'Snap Point on Path');
    
    
    % Connect closest point on path to query point
    plot(...
        [point(:,1) closest_path_point(:,1)],...
        [point(:,2) closest_path_point(:,2)],'g-','Linewidth',2);
    
    
end % Ends the flag_do_debug if statement



end % Ends the function

