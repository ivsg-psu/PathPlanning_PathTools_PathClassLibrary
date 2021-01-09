function [distance,location] = fcn_Path_findProjectionHitOntoPath(path,sensor_vector_start,sensor_vector_end,varargin)   
% fcn_Path_findProjectionHitOntoPath calculates hits between a sensor
% projection and a path, returning the distance and location of the hit.
%
% FORMAT: 
%
%      [distance,location] = ...
%         fcn_Path_findProjectionHitOntoPath(path,...
%         sensor_vector_start,sensor_vector_end,...
%         (flag_search_type),(fig_num))  
%
% INPUTS:
%
%      path: an N x 2 vector containing the X,Y points of the path to be
%      checked for intersections
%
%      sensor_vector_start: a 1 x 2 vector containing the X,Y points of the
%      sensor's start location
%
%      sensor_vector_end: a 1 x 2 vector containing the X,Y points of the
%      sensor's end location
%
%      (OPTIONAL INPUTS)
%      flag_search_type: an integer specifying the type of search.
%
%            0: return distance and location only if the given
%            sensor_vector overlaps the path (this is the default)
%
%            1: return distane and location if any projection of the sensor
%            vector, in any direction, hits the path (in other words, if
%            there is any intersection). Note that distance returned will
%            be negative if the nearest intersection is in the opposite
%            direction of the given sensor vector.
%
%      fig_num: a figure number to plot results. Turns debugging on.
%
% OUTPUTS:
%
%      distance: a 1 x 1 scalar representing the distance to the closest
%      intersection of the sensor with the path
%
%      location: a 1 x 2 vector of the X,Y location of intersection point
%
% EXAMPLES:
%      
%       See the script: script_test_fcn_Path_findProjectionHitOntoPath.m
%       for a full test suite. 
%
% Adopted from https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
% This function was written on 2020_11_14 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%      2020_11_14 
%      - wrote the code
%      2020_12_29 
%      - fixed bug with missing sensor vector if starts on a path
%      - added better comments
%      2020_12_31 
%      - fixed the input arguments so that they are more clear (start/end)
%      - fixed the incorrect location of the debug echo at top of the code
%      - fixed flag usage to decouple plotting with debugging
%      2021_01_08 
%      -- Added input check on path type
 

%% Set up for debugging
flag_do_debug = 0; % Flag to plot the results for debugging
flag_do_plot = 0; % Flag to plot the results for debugging
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
% Set default values
flag_search_type = 0;

% check input arguments
if flag_check_inputs == 1
    if nargin < 3 || nargin > 5
        error('Incorrect number of input arguments.')
    end
    
    % Check path input
    fcn_Path_checkInputsToFunctions(path, 'path');    
end

if 4 <= nargin
    flag_search_type = varargin{1};
end

if 5 == nargin
    fig_num = varargin{2};
    figure(fig_num);
    flag_do_plot = 1;
end

if flag_do_debug
    flag_do_plot = 1;
end


%% Calculations begin here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define each path as a set of walls that can be hit
wall_start = [path(1:end-1,1) path(1:end-1,2)];
wall_end   = [path(2:end,1) path(2:end,2)];

% Define p, q, r and s vectors
p = wall_start;
q = sensor_vector_start;
r = wall_end - wall_start;
s = sensor_vector_end - sensor_vector_start;

% Define useful intermediate terms
r_cross_s = crossProduct(r,s);
q_minus_p =  q - p;
q_minus_p_cross_s = crossProduct(q_minus_p,s);
q_minus_p_cross_r = crossProduct(q_minus_p,r);

% Are any of these parallel?
parallel_indices = find(0==r_cross_s);
if any(parallel_indices)
    r_cross_s(parallel_indices) = 1; % They are colinear or parallel, so make dummy length
end

% Calculate t and u, where t is distance along the path, and u is distance
% along the sensor.

t = q_minus_p_cross_s./r_cross_s; % Distance along the path
u = q_minus_p_cross_r./r_cross_s; % Distance along the sensor

% Fix any situations that are parallel, as these will give wrong
% calculations
t(parallel_indices) = inf;
u(parallel_indices) = inf;

% Initialize all intersections to infinity
intersections = NaN*ones(length(p(:,1)),2);

% Note: could speed this up with nested if logical statements, but only if
% are doing checks on single segments at a time. Since doing many segments
% at once, need to use vector form.
if 0 == flag_search_type
    good_vector = ((0<=t).*(1>=t).*(0<=u).*(1>=u));
elseif 1 == flag_search_type
    good_vector = ((0<=t).*(1>=t));
else
    error('Incorrect flag_search_type entered');
end

good_indices = find(good_vector>0);

if ~isempty(good_indices)
    result = p + t.*r; 
    intersections(good_indices,:) = result(good_indices,:);    
end

% Find the distances via Euclidian distance to the sensor's origin
% note: a faster way to do this might be to just
% calculate t*r as a length
distances_squared = sum((intersections - sensor_vector_start).^2,2);

% Keep just the minimum distance
[closest_distance_squared,closest_index] = min(distances_squared);

distance = closest_distance_squared^0.5*sign(u(closest_index));
location = intersections(closest_index,:);


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
if flag_do_plot
    
    % Set up the figure
    figure(fig_num);
    clf;
    hold on;
    axis equal;
    grid on; grid minor;
    
    % Plot the path in black
    plot(path(:,1),path(:,2),'k.-','Linewidth',1);
    handle_text = text(path(1,1),path(1,2),'Path');
    set(handle_text,'Color',[0 0 0]);
       
    % Plot the sensor vector
    quiver(q(:,1),q(:,2),s(:,1),s(:,2),'r','Linewidth',3);
    handle_text = text(q(:,1),q(:,2),'Sensor');
    set(handle_text,'Color',[1 0 0]);
    
    axis_size = axis;
    y_range = axis_size(4)-axis_size(3);
        
    % Plot any hits in blue
    plot(location(:,1),location(:,2),'bo','Markersize',30);
    handle_text = text(location(:,1),location(:,2)-0.05*y_range,sprintf('Hit at distance: %.2f',distance));
    set(handle_text,'Color',[0 0 1]);
    
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end
end



%% Calculate cross products
function result = crossProduct(v,w)
result = v(:,1).*w(:,2)-v(:,2).*w(:,1);
end



