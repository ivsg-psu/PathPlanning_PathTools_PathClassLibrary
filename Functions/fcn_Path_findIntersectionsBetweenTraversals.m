function [intersection_points,...
    s_coordinates_in_traversal_1,...
    s_coordinates_in_traversal_2] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    traversal_1,...
    traversal_2, ...
    varargin)
% fcn_Path_findIntersectionsBetweenTraversals
% Given two traversals, finds the intersection points between the traverals
% and returns the results as points, station coordinates in 1, and station
% coordinates in 2
%
% FORMAT:
%
%      [intersection_points,...
%       s_coordinates_in_traversal_1,...
%       s_coordinates_in_traversal_2] = ...
%        fcn_Path_findIntersectionsBetweenTraversals(...
%        traversal_1,traversal_2,...
%        (fig_num));
%
% INPUTS:
%
%      traversal_1: a traversal structure that specifies the path,
%      s-coordinates, etc of the first traveral
%
%      traversal_2: a traversal structure that specifies the path,
%      s-coordinates, etc of the second traveral
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      intersection_points: a Mx2 vector containing the [X Y] locations of
%      the M intersections between two traversals. Returns an empty matrix
%      if there are no intersections
%
%      s_coordinates_in_traversal_1: a scalar (1x1) representing the
%      s-coordinate distance along the first traversal where the
%      intersection(s) take place
%
%      s_coordinates_in_traversal_2: a scalar (1x1) representing the
%      s-coordinate distance along the second traversal where the
%      intersection(s) take place
%
% DEPENDENCIES:
%
%      fcn_Path_checkInputsToFunctions
%      fcn_Path_plotTraversalsXY
%
% EXAMPLES:
%
% See the script: script_test_fcn_Path_findIntersectionsBetweenTraversals
% for a full test suite.
%
% This function was written on 2021_01_23 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
%      2021_01_23:
%      -- first write of the code


flag_do_debug = 0; % Flag to debug the results
flag_do_plot = 0; % Flag to plot the results
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
% (station,central_traversal,nearby_traversal, (flag_projection_type?))


if flag_check_inputs == 1
    % Are there the right number of inputs?
    if nargin < 2 || nargin > 3
        error('Incorrect number of input arguments')
    end
    
    % Check the traversal_1 input
    fcn_Path_checkInputsToFunctions(traversal_1, 'traversal');
    
    % Check the traversal_2 input
    fcn_Path_checkInputsToFunctions(traversal_2, 'traversal');
              
end


% Does user want to show the plots?
if 3 == nargin
    fig_num = varargin{1};
    figure(fig_num);
    flag_do_plot = 1;
else
    if flag_do_debug
        fig = figure;
        fig_num = fig.Number;
        flag_do_plot = 1;
    end
end

%% Start of main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through all the segments in traversal 1, checking each for
% intersections with traversal 2
path_of_traversal_1 = [traversal_1.X traversal_1.Y];
path_of_traversal_2 = [traversal_2.X traversal_2.Y];

% Set the flag to indicate that we will allow multiple crossings and
% request all this information
flag_search_type = 2;

% Initialize values
raw_intersection_points = [];
raw_s_coordinates_in_traversal_1 = [];
raw_s_coordinates_in_traversal_2 = [];
raw_all_path_segments_1 = [];
raw_all_path_segments_2 = [];

for i_segment_path_1 = 2:length(traversal_1.X)
    sensor_vector_start = [...
        traversal_1.X(i_segment_path_1-1) traversal_1.Y(i_segment_path_1-1)];
    sensor_vector_end = [...
        traversal_1.X(i_segment_path_1) traversal_1.Y(i_segment_path_1)];
    
    % Check to see if there are intersections
    [distance,location,path_segments] = ...
        fcn_Path_findProjectionHitOntoPath(...
        path_of_traversal_2,...
        sensor_vector_start,sensor_vector_end,...
        flag_search_type);
    
    % Did we hit anything? If so, save it
    if ~isempty(distance)
        raw_intersection_points = [raw_intersection_points; location];  %#ok<*AGROW>
        raw_all_path_segments_1 = [raw_all_path_segments_1; (i_segment_path_1-1)*ones(length(path_segments(:,1)),1)];
        raw_all_path_segments_2 = [raw_all_path_segments_2; path_segments];
        
        % For each of the hits, find the s-coordinates of the hits
        for ith_hit = 1:length(location(:,1))
            point = location(ith_hit,:);
            
            % Find the s-coordinate for traversal 1
            [~,s_coordinate_1_this_hit,first_path_point_index,second_path_point_index,percent_along_length] = ...
                fcn_Path_snapPointToPathViaVectors(point, path_of_traversal_1); %#ok<ASGLU>
            
            % Interpolate the s-coordinate in 1
            s_diff_1_this_hit = ...
                traversal_1.Station(i_segment_path_1) - ...
                traversal_1.Station(i_segment_path_1-1);
            s_diff_calc = sum((sensor_vector_end - sensor_vector_start).^2,2).^0.5;
            
            %             if s_diff_1_this_hit ~= s_diff_calc(1,1)
            %                 warning('S-coordiantes in trajectory 1 do not appear to match distances between points');
            %             end
            
            % Percent along length
            travel = sum((point - sensor_vector_start).^2,2).^0.5;
            percent_along_length = travel/s_diff_calc;
            
            s_coordinate_1_this_hit_interp = ...
                traversal_1.Station(i_segment_path_1 - 1) + ...
                s_diff_1_this_hit * percent_along_length;
        
            %             if (s_coordinate_1_this_hit_interp ~= s_coordinate_1_this_hit)
            %                 warning('Interpolated and calculated s-coordinate distances did not match. This usually indicates a poorly formed traversal or a situation where a traversal loops back onto itself.');
            %             end
            
            % Update the s-coordinates for traversal 1
            raw_s_coordinates_in_traversal_1 = [...
                raw_s_coordinates_in_traversal_1;
                s_coordinate_1_this_hit_interp];
            

            % Find the s-coordinate for traversal 2
            [~,s_coordinate_2_this_hit,first_path_point_index,second_path_point_index,percent_along_length] = ...
                fcn_Path_snapPointToPathViaVectors(point, path_of_traversal_2); %#ok<ASGLU>

            % Interpolate the s-coordinate in 2
            s_diff_2_this_hit = ...
                traversal_2.Station(path_segments(ith_hit)+1) - ...
                traversal_2.Station(path_segments(ith_hit));            
            
            % Find where the second point started
            second_point = [traversal_2.X(path_segments(ith_hit)) traversal_2.Y(path_segments(ith_hit))];
            
            % Percent along length
            travel = sum((point - second_point).^2,2).^0.5;
            percent_along_length = travel/s_diff_2_this_hit;

            % Interpolate the s-coordinate
            s_coordinate_2_this_hit_interp = ...
                traversal_2.Station(path_segments(ith_hit)) + ...
                s_diff_2_this_hit * percent_along_length;
            
            %             % Do they match?
            %             if (s_coordinate_2_this_hit_interp ~= s_coordinate_2_this_hit)
            %                 warning('Interpolated and calculated s-coordinate distances did not match. This usually indicates a poorly formed traversal or a situation where a traversal loops back onto itself.');
            %             end
            
            % Update the s-coordinates for traversal 2
            raw_s_coordinates_in_traversal_2 = [...
                raw_s_coordinates_in_traversal_2;
                s_coordinate_2_this_hit_interp];

            
        end % Ends for loop through the points that were "hits"
        
    end % Ends check to see if distance is empty
        
end % Ends looping through segments

% Remove repeats
if ~isempty(raw_intersection_points)
    % USE THIS IF WE HAVE TO USE THE SEGMENT NUMBERS AS OUTPUT.
    % We don't do this yet, but might have a situation where we need to
    % know which segments are intersecting in each path. If that's the
    % case, we'll need the outputs of path_segments_1 and path_segments_2
    %     all_data = [raw_intersection_points ...
    %         raw_s_coordinates_in_traversal_1 ...
    %         raw_s_coordinates_in_traversal_2 ...
    %         raw_all_path_segments_1 ...
    %         raw_all_path_segments_2];
    %
    all_data = [raw_intersection_points ...
        raw_s_coordinates_in_traversal_1 ...
        raw_s_coordinates_in_traversal_2];
    
    % Keep only the data that is unique, and sort it too
    [cleaned_data,indices_original,~] = unique(all_data,'rows','sorted');
    %  [cleaned_data,~,~] = unique(all_data,'rows','sorted');
    
    % Pull out data
    intersection_points          = cleaned_data(:,1:2);
    s_coordinates_in_traversal_1 = cleaned_data(:,3);
    s_coordinates_in_traversal_2 = cleaned_data(:,4);
    all_path_segments_1          = raw_all_path_segments_1(indices_original,1);
    all_path_segments_2          = raw_all_path_segments_2(indices_original,1);
    
    %     % USE THIS IF KEEPING SEGMENT NUMBERS AS OUTPUTS
    %     all_path_segments_1          = cleaned_data(:,5);
    %     all_path_segments_2          = cleaned_data(:,6);
    %
else
    intersection_points          = [];
    s_coordinates_in_traversal_1 = [];
    s_coordinates_in_traversal_2 = [];
    all_path_segments_1          = [];
    all_path_segments_2          = [];
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
if flag_do_plot
    figure(fig_num);
    clf;
    hold on;
    grid on;
    
    % Plot the two traversals
    clear data
    data.traversal{1} = traversal_1;
    data.traversal{2} = traversal_2;       
    fcn_Path_plotTraversalsXY(data,fig_num);
    
    if ~isempty(intersection_points)
        % Plot the hit points
        plot(intersection_points(:,1),intersection_points(:,2),'ro','Markersize',20);
        
        % Label the points
        for ith_hit = 1:length(intersection_points(:,1))
            text(intersection_points(ith_hit,1),intersection_points(ith_hit,2),sprintf(...
                '%.0d: S1 is %.1f in segment %.0d, S2 is %.1f in segment %.0d',...
                ith_hit, ...
                s_coordinates_in_traversal_1(ith_hit),...
                all_path_segments_1(ith_hit),...
                s_coordinates_in_traversal_2(ith_hit),...
                all_path_segments_2(ith_hit)));
        end
        
        
        legend('Traversal 1','Traversal 2','Intersections');
        title(sprintf('%.0d intersections total',length(intersection_points(:,1))));
    else
        legend('Traversal 1','Traversal 2');
        title('No intersections detected');
    end
    
    % % Label the points with distances?
    % for i_point = 1:length(path(:,1))
    %     text(path(i_point,2),path(i_point,3),sprintf('%.2f',distances_point_to_path(i_point)));
    % end
    
    
    
    
end % Ends the flag_do_plot if statement

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends the function

