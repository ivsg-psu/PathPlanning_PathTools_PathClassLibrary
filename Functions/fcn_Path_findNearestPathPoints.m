function closest_path_point_indicies = ...
    fcn_Path_findNearestPathPoints(query_points, path,varargin)
% fcn_Path_findNearestPathPoints
% Finds the closest_path_point_indicies repressenting the indices of the
% path that are closest to a set of query points
%
% The solution method is vectorized by finding the minimum of the squared
% distances of each query point to each path point. Thus, if there are M
% query points and N path points, there are N x M calculations of squared
% distances, and M repetitions of minmum sorting operations across N data
% sets. This can be time consuming, so avoid calling this function
% with a large path or with very large numbers of query points.
% 
% FORMAT: 
%
%     closest_path_point_indicies = ...
%     fcn_Path_findNearestPathPoints(query_points, path, (fig_num))
%
% INPUTS:
%
%     query_points: a Mx2 or Mx3 vector containing the [X Y (Z)] locations
%     of query points where M is the number of rows of different points
%
%     path: a Nx2 or Nx3 vector of [X Y (Z)] path points, where N is the
%     number of points the points on the path, with N >= 2.
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
%      closest_path_point_indicies: a Mx1 vector containing the indices of
%      the path that is nearest to each of the query_points.
%
% EXAMPLES:
%      
% See the script: script_test_fcn_Path_findNearestPathPoints
% for a full test suite.
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_Path_convertPathToTraversalStructure
%
% This function was written on 2023_06_02 by S. Brennan in support of the
% fcn_Path_snapPointOntoNearestPath function expansion.
%
% Questions or comments? sbrennan@psu.edu 

% Revision history:   
% 2023_06_02 by sbrennan@psu.edu
% -- first write of the code
% 2025_06_23 - S. Brennan
% -- Updated debugging and input checks

% TO-DO
% (none)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==3 && isequal(varargin{end},-1))
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
        narginchk(2,3);

        % Check the data input
        fcn_DebugTools_checkInputsToFunctions(path, 'path2or3D');

        % Check that the dimension of the point and path match
        if length(query_points(1,:)) ~= length(path(1,:))
            error('The dimension of the query_points input, in number of columns, must match the column dimension of the path');
        end
    end
end

% Does user want to show the plots?
flag_do_plots = 0; % Default is to NOT show plots
if (0==flag_max_speed) && (3 == nargin) 
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        fig_num = temp;
        figure(fig_num);
        flag_do_plots = 1;
    end
else
    if flag_do_debug
        fig_debug = 4848; %#ok<NASGU>
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The solution method is as follows:
%  1. If there is just one query point, calculate that point alone
%  2. If there are more than one query point, then 
%  a. duplicate the N path points by the M query points to create a NxM
%  column vector
%  b. duplicate the M query points by the N path points to create a NxM
%  column vector
%  c. Take the difference between the two vectors, and square this
%  difference.
%  d. Rearrange the differences into M rows by N columns, 
%  e. Find the minimum along each row dimension, saving the index of
%  minimum.


%  1. If there is just one query point, calculate that point alone
if isscalar(query_points(:,1))
    squared_distances_point_to_path = sum((path - query_points).^2,2); % Calculate distances
    [~,closest_path_point_indicies] = min(squared_distances_point_to_path);  % Grab index of the closest point
else   
    %  2. If there are more than one query point, then
    %  a. duplicate the N path points by the M query points to create a NxM
    %  column vector
    
    % How many duplications do we need?
    Mqueries = length(query_points(:,1));
    path_points_duplicated = repmat(path,Mqueries,1);
    
    %  b. duplicate the M query points by the N path points to create a NxM
    %  column vector
    Npathpoints = length(path(:,1));
    query_points_indicies = (1:Mqueries);    
    query_points_indicies_duplicated_matrix = repmat(query_points_indicies,Npathpoints,1);
    query_points_Indicies_duplicated = reshape(query_points_indicies_duplicated_matrix,Npathpoints*Mqueries,1);
    query_points_duplicated = query_points(query_points_Indicies_duplicated,:);
    
    %  c. Take the difference between the two vectors, and square this
    %  difference.
    vector_differences_squared = sum((path_points_duplicated-query_points_duplicated).^2,2);
    
    %  d. Rearrange the differences into M rows by N columns,
    distance_squared_transposed = reshape(vector_differences_squared,Npathpoints,Mqueries);
    distances_squared_matrix = distance_squared_transposed';
    
    %  e. Find the minimum along each row dimension, saving the index of
    %  minimum.
    [~,closest_path_point_indicies] = min(distances_squared_matrix,[],2);  %#ok<UDIM> % Grab index of the closest point along column direction
end


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
    clf;
    hold on;
    grid on;
    
    % Is this a 2-D query
    if length(path(1,:))==2
        % Plot the path
        plot(path(:,1),path(:,2),'r-','Linewidth',3);
        plot(path(:,1),path(:,2),'r.','Markersize',20);
        
        axis equal;
        
        % Plot the query points
        plot(query_points(:,1),query_points(:,2),'k.','Markersize',20);
        % text(query_points(:,1),query_points(:,2),'Query point');
        
        % Show vectors from the query points to closest path points
        for ith_query = 1:length(query_points(:,1))
            closest_point_vector = path(closest_path_point_indicies(ith_query,1),:) - query_points(ith_query,:);
            quiver(query_points(ith_query,1),query_points(ith_query,2),...
                closest_point_vector(1,1),closest_point_vector(1,2),...
                0,'-','Linewidth',3);
        end

        
    elseif length(path(1,:))==3
        
        % Plot the path
        plot3(path(:,1),path(:,2),path(:,3),'r-','Linewidth',3);
        plot3(path(:,1),path(:,2),path(:,3),'r.','Markersize',20);
        
        axis equal;
        
        % Plot the query points
        plot3(query_points(:,1),query_points(:,2),query_points(:,3),'ko');
        % text( query_points(:,1),query_points(:,2),query_points(:,3),'Query point');
        
        % Show vectors from the query points to closest path points
        for ith_query = 1:length(query_points(:,1))
            closest_point_vector = path(closest_path_point_indicies(ith_query,1),:) - query_points(ith_query,:);
            quiver3(query_points(ith_query,1),query_points(ith_query,2),query_points(ith_query,3),...
                closest_point_vector(1,1),closest_point_vector(1,2),closest_point_vector(1,3),...
                0,'-','Linewidth',3);
        end
        

    end
end % Ends the flag_do_debug if statement



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

% %% fcn_INTERNAL_calculateVectorMeasures    
% function     [closest_path_point,s_coordinate,percent_along_length,unit_orthogonal_projection_vector] = fcn_INTERNAL_calculateVectorMeasures(point,segment_end_point,segment_start_point,segment_start_station)
% % Do the dot products - define the vectors first
% % See: https://mathinsight.org/dot_product for explanation Basically,
% % we are seeing what amount the point_vector points in the direction of
% % the path_vector.
% 
% 
% path_vector  = segment_end_point - segment_start_point;
% path_segment_length  = sum(path_vector.^2,2).^0.5;
% unit_orthogonal_projection_vector = path_vector/path_segment_length*[0 1; -1 0]; % Rotate by 90 degrees
% 
% start_to_point_vector = point-segment_start_point;
% projection_distance  = dot(path_vector,start_to_point_vector)/path_segment_length; % Do dot product
% % offset_distance  = dot(ortho_path_vector,start_to_point_vector)/path_segment_length; % Do dot product
% percent_along_length = projection_distance/path_segment_length;
% 
% % Calculate the remaining outputs
% closest_path_point = segment_start_point + path_vector*percent_along_length;
% s_coordinate       = segment_start_station + path_segment_length*percent_along_length;
% end % Ends fcn_INTERNAL_calculateVectorMeasures
% 
% 
% %% fcn_INTERNAL_convertOffsetCoordinateToImaginaryNumber
% function     [real_distance, imag_distance] = ...
%     fcn_INTERNAL_convertOffsetCoordinateToImaginaryNumber(...
%     unit_projection_vector_of_real_axis,...
%     origin_point, ...
%     point_to_convert)
% 
% % This function is used to convert a point's location into an imaginary
% % distance, where the real portion of the distance represents the distance in
% % the direction of a given unit vector, and the imaginary portion of the
% % distance represents the orthogonal component. The conversion assumes a
% % coordinate orientation of the Re/Im plane such that the real axis is
% % aligned with the given unit projection vector representing the real axis.
% %
% % The reason for this function is that, when snapping a point onto a
% % reference traversal, there are situations - particularly before the start
% % of the traversal, after the end of the traversal, and along sharp corners
% % of the traversal, where the point cannot be correctly represented solely
% % by the orthogonal projection distance. For example, in sharp corners, the
% % orthogonal projection is undefined and must be chosen by the user. The
% % conversion of orthogonal distance thus allows correct and invertible
% % converstions from path offsets in XY to Sy representations, and vice
% % versa.
% 
% unit_projection_vector_of_imag_axis = unit_projection_vector_of_real_axis*[0 1; -1 0]; % Rotate by 90 degrees
% 
% start_to_point_vector = point_to_convert-origin_point;
% real_distance  = dot(unit_projection_vector_of_real_axis,start_to_point_vector); % Do dot product
% imag_distance  = dot(unit_projection_vector_of_imag_axis,start_to_point_vector); % Do dot product
% 
% % % Measure the distance from the closest_path_point station
% % distance_imaginary = distance_imaginary - path_station(closest_path_point_index,1);
% % 
% % % Need to flip the sign so that the projection is correctly positioned
% % % relative to the origin
% % distance_imaginary = -1.0*distance_imaginary;
% end % Ends fcn_INTERNAL_convertOffsetCoordinateToImaginaryNumber
% 

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
