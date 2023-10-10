function [distance,location,path_segment, t, u] = fcn_Path_findProjectionHitOntoPath(path,sensor_vector_start,sensor_vector_end,varargin)   
% fcn_Path_findProjectionHitOntoPath calculates hits between a sensor
% projection and a path, returning the distance and location of the hit.
%
% FORMAT: 
%
%      [distance,location,path_segment, t, u] = ...
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
%            0: return distance and location of first intersection only if
%            the given sensor_vector overlaps the path (this is the
%            default)
%
%            1: return distance and location of first intersection if any
%            projection of the sensor vector, in any direction, hits the
%            path (in other words, if there is any intersection). Note that
%            distance returned will be negative if the nearest intersection
%            is in the opposite direction of the given sensor vector.
%
%            2: returns distances and locations as M x 1 and M x 2 vectors
%            respectively, where the M rows represent ALL the detected
%            intersections. In cases where the sensor vector completely
%            overlaps a path segment, the start and end of overlap are
%            given.
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
%      path_segment: the segment number of the path that was hit (1 is the
%      first segment, 2 is the second, etc)
%
%       t and u: t is distance along the path, and u is distance
%       along the sensor, as fractions of the input vector lengths.    
%
% DEPENDENCIES:
%
%      fcn_Path_checkInputsToFunctions
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
%      2021_01_23
%      -- Added flag_search_type = 2 option, to allow multiple cross points
%      to be returned
%      -- Fixed bug with partially overlapping vectors not returning a
%      result
%      -- Added path segment output so that we can ID which segment was hit
%      2021_01_24
%      -- Fixed bug with overlapping colinear where two path segments
%      identified when there is only one
%      2021_12_27
%      -- Added better comments on flags

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

% Does user wish to specify search type?
if 4 <= nargin
    flag_search_type = varargin{1};
end


% Does user want to show the plots?
if 5 == nargin
    fig_num = varargin{2};
    figure(fig_num);
    flag_do_plot = 1;
else
    if flag_do_debug
        fig = figure;
        fig_num = fig.Number;
        flag_do_plot = 1;
    end
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
path_segments = (1:length(wall_start(:,1)))';

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
parallel_non_intersecting_indices = find((0==r_cross_s).*(0~=q_minus_p_cross_r));
if any(parallel_non_intersecting_indices)
    r_cross_s(parallel_non_intersecting_indices) = 1; % They are colinear or parallel, so make dummy length
end

% Are any of these colinear?
colinear_indices = find((0==r_cross_s).*(0==q_minus_p_cross_r));
if any(colinear_indices)
    r_cross_s(colinear_indices) = 1; % They are colinear or parallel, so make dummy length
    r_dot_r = sum(r.*r,2);
    q_minus_p_dot_r = sum(q_minus_p.*r,2);
    s_dot_r = sum(s.*r,2);
    t0 = q_minus_p_dot_r./r_dot_r;
    t1 = t0 + s_dot_r./r_dot_r;
    
    % Keep only the good indices
    %     % For debugging:
    %     conditions = [-0.5 -0.4; -0.5 0; -0.5 .2; 0 0.2; 0.2 0.4; 0.2 1; 0.2 1.2; 1 1.2; 1.2 1.3; -0.5 1.2]
    %     t0 = conditions(:,1);
    %     t1 = conditions(:,2);
    %     t0_inside = (t0>=0)&(t0<=1);
    %     t1_inside = (t1>=0)&(t1<=1);
    %     t0_t1_surround = (t0<0)&(t1>1) | (t1<0)&(t0>1);
    %     any_within = t0_inside | t1_inside | t0_t1_surround;
    %     fprintf(1,'t0 inside flag:\n');
    %     [conditions t0_inside]
    %     fprintf(1,'t1 inside flag:\n');
    %     [conditions t1_inside]
    %     fprintf(1,'surround flag:\n');
    %     [conditions t0_t1_surround]
    %     fprintf(1,'any_wighin flag:\n');
    %     [conditions any_within]
    
    % Check whether there is overlap by seeing of the t0 and t1 values are
    % within the interval of [0 1], endpoint inclusive
    t0_inside = (t0>=0)&(t0<=1);
    t1_inside = (t1>=0)&(t1<=1);
    t0_t1_surround = (t0<0)&(t1>1) | (t1<0)&(t0>1);
    any_within = t0_inside | t1_inside | t0_t1_surround;
    good_indices = find(any_within);
    good_colinear_indices = intersect(colinear_indices,good_indices);
    
    % Fix the ranges to be within 0 and 1
    t0(good_colinear_indices) = max(0,t0(good_colinear_indices));
    t0(good_colinear_indices) = min(1,t0(good_colinear_indices));

    t1(good_colinear_indices) = max(0,t1(good_colinear_indices));
    t1(good_colinear_indices) = min(1,t1(good_colinear_indices));

end

% Calculate t and u, where t is distance along the path, and u is distance
% along the sensor.

t = q_minus_p_cross_s./r_cross_s; % Distance along the path
u = q_minus_p_cross_r./r_cross_s; % Distance along the sensor

% Fix any situations that are parallel and non-intersecting, as these will
% give wrong calculation results from the t and u calculations above
t(parallel_non_intersecting_indices) = inf;
u(parallel_non_intersecting_indices) = inf;

% Fix any situations that are colinear. For these, we save the start point
% as the point where overlap starts, and the end point where overlap ends.
% Note that this creates NEW intersections beyond the number of segements
if any(colinear_indices)
    % Shut off colinear ones to start
    t(colinear_indices) = inf;
    u(colinear_indices) = inf;
    
    % Correct the t values
    u(good_colinear_indices) = 1;
    t(good_colinear_indices) = t0(good_colinear_indices);

    % Do we need to add more hit points?
    indices_hit_different_point = find(t0~=t1);
    more_indices = intersect(indices_hit_different_point,good_colinear_indices);
         
    % Make p and r, t and u longer so that additional hit points are
    % calculated in special case of overlaps
    p = [p; p(more_indices,:)];
    r = [r; r(more_indices,:)];
    u = [u; u(more_indices)];
    t = [t; t1(more_indices)];   
    path_segments = [path_segments; path_segments(more_indices)];

end

% Initialize all intersections to infinity
intersections = NaN*ones(length(p(:,1)),2);

% Note: could speed this up with nested if logical statements, but only if
% are doing checks on single segments at a time. Since doing many segments
% at once, need to use vector form.
if 0 == flag_search_type
    good_vector = ((0<=t).*(1>=t).*(0<=u).*(1>=u));
elseif 1 == flag_search_type
    good_vector = ((0<=t).*(1>=t));
elseif 2 == flag_search_type
    good_vector = ((0<=t).*(1>=t).*(0<=u).*(1>=u));
else
    error('Incorrect flag_search_type entered');
end

% Keep only the indices that work
good_indices = find(good_vector>0);

if ~isempty(good_indices)
    result = p + t.*r; 
    intersections(good_indices,:) = result(good_indices,:);    
end

% Find the distances via Euclidian distance to the sensor's origin
% note: a faster way to do this might be to just
% calculate t*r as a length
distances_squared = sum((intersections - sensor_vector_start).^2,2);

if flag_search_type ~=2
    % Keep just the minimum distance
    [closest_distance_squared,closest_index] = min(distances_squared);
    
    distance = closest_distance_squared^0.5*sign(u(closest_index));
    location = intersections(closest_index,:);
    path_segment = path_segments(closest_index);
else
    % Return all the results
    good_indices = find(~isnan(distances_squared));
    distance = distances_squared(good_indices).^0.5.*sign(u(good_indices));
    location = intersections(good_indices,:);
    path_segment = path_segments(good_indices);
end

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
    plot(path(:,1),path(:,2),'k.-','Linewidth',5);
    handle_text = text(path(1,1),path(1,2),'Path');
    set(handle_text,'Color',[0 0 0]);
       
    % Plot the sensor vector
    quiver(q(:,1),q(:,2),s(:,1),s(:,2),'r','Linewidth',3);
    plot(sensor_vector_end(:,1),sensor_vector_end(:,2),'r.','Markersize',10);
    
    handle_text = text(q(:,1),q(:,2),'Sensor');
    set(handle_text,'Color',[1 0 0]);
    
    axis_size = axis;
    y_range = axis_size(4)-axis_size(3);
        
    % Plot any hits in blue
    for i_result = 1:length(distance)
        plot(location(i_result,1),location(i_result,2),'bo','Markersize',30);
        handle_text = text(location(i_result,1),location(i_result,2)-0.05*y_range,sprintf('Hit at distance: %.2f',distance(i_result)));
        set(handle_text,'Color',[0 0 1]);
    end
    
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end
end



%% Calculate cross products
function result = crossProduct(v,w)
result = v(:,1).*w(:,2)-v(:,2).*w(:,1);
end



