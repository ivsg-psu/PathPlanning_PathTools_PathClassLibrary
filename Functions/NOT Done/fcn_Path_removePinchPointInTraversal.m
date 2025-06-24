function [traversal_no_pinch_point] = ...
    fcn_Path_removePinchPointInTraversal(...
    traversal_with_pinch_point,...
    varargin)
% fcn_Path_removePinchPointInTraversal
% Given a traversal with a pinch point - an area where the traversal
% suddenly bends back on itself before continuing - this function removes
% the pinch point
%
% FORMAT:
%
%     [traversal_no_pinch_point] = ...
%         fcn_Path_removePinchPointInTraversal(...
%         traversal_with_pinch_point
%        (fig_num));
%
% INPUTS:
%
%      traversal_with_pinch_point: a traversal structure that specifies the
%      path, s-coordinates, etc of a traveral with a pinch point
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
%      traversal_no_pinch_point: a traversal structure that specifies the
%      path, s-coordinates, etc of a traveral with the pinch points
%      removed.
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_Path_convertPathToTraversalStructure
%      fcn_Path_plotTraversalsXY
%
% EXAMPLES:
%
% See the script: script_test_fcn_Path_removePinchPointInTraversal
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
    if nargin < 1 || nargin > 2
        error('Incorrect number of input arguments')
    end
    
    % Check the traversal_with_pinch_point input
    fcn_DebugTools_checkInputsToFunctions(traversal_with_pinch_point, 'traversal');
    
end


% Does user want to show the plots?
if 2 == nargin
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

Nsegments = length(traversal_with_pinch_point.X) - 1;

if Nsegments<2
    % There needs to be at least 2 segments to self-intersect
    traversal_no_pinch_point = traversal_with_pinch_point;
else    
    
    % Loop through all the segments in traversal 1, checking each for
    % intersections with traversal 2
    original_path_of_traversal = [traversal_with_pinch_point.X traversal_with_pinch_point.Y];
    
    % Set the flag to indicate that we are ONLY searching if the exact segment
    % crosses or not
    flag_search_type = 0;
    
    % Initialize values
    no_pinch_path = original_path_of_traversal(1,:);
    remaining_path = original_path_of_traversal(2:end,:);
    
    % Are there at least 3 points in the remaining path?
    while length(remaining_path(:,1))>=2
        % Define the sensor vector
        sensor_vector_start = no_pinch_path(end,:);
        sensor_vector_end = remaining_path(1,:);
        sensor_vector_length = sum((sensor_vector_end - sensor_vector_start).^2,2).^0.5;
                
        % Check to see if there are intersections
        [distance,hit_location,path_segments] = ...
            fcn_Path_findProjectionHitOntoPath(...
            remaining_path,...
            sensor_vector_start,sensor_vector_end,...
            flag_search_type);
        
        % Did we hit anything? If so, save it and set a flag that a pinch
        % point was hit!
        if isnan(distance) || (distance==sensor_vector_length(1,1) && path_segments==1)
            % Nothing hit
            no_pinch_path = [no_pinch_path; sensor_vector_end];  %#ok<AGROW>
            remaining_path = remaining_path(2:end,:);

        else % Hit something
            
            % Calculate how much s-distance is being "cut" to see if we
            % need to warn the user:

            % First, interpolate the s-coordinate after the 1st hit            
            travel = sum((sensor_vector_end - hit_location).^2,2).^0.5;
                       
            % Second: find the s-coordinate for traversal after this hit
            [~,s_coordinate_after_hit,~,~,~] = ...
                fcn_Path_snapPointToPathViaVectors(...
                hit_location, remaining_path); 
            
            % Third: add these values
            s_coordinates_trimmed = (travel + s_coordinate_after_hit);
            
            % Do we need to warn the user?
            if s_coordinates_trimmed > 10
                warning('10 meters or more were trimmed. This is a large amount!');
            end
            
            % Update the path                        
            remaining_path = [hit_location; remaining_path(path_segments+1:end,:)];
            
            % Make sure we didn't just repeat a point in the path (which
            % happens if the loop comes back onto itself)
            if isequal(remaining_path(1,:),remaining_path(2,:))
                remaining_path = remaining_path(2:end,:);
            end
            
        end % Ends check to see if distance is empty
        
    end % Ends while loop
    no_pinch_path = [no_pinch_path; remaining_path]; 
    
    % Clean up the path by removing repeats - these can occur when the path
    % loops back onto itself
    cleaned_path = no_pinch_path(1,:);
    for ith_row = 2:length(no_pinch_path)
        if ~isequal(no_pinch_path(ith_row,:),no_pinch_path(ith_row-1,:))
            cleaned_path = [cleaned_path; no_pinch_path(ith_row,:)]; %#ok<AGROW>
        end
    end % Ends the for loop over rows, to clean repeats
    no_pinch_path = cleaned_path;

    traversal_no_pinch_point = fcn_Path_convertPathToTraversalStructure(no_pinch_path);
end % Ends check to see if there are at least 3 segments




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
    data.traversal{1} = traversal_with_pinch_point;
    data.traversal{2} = traversal_no_pinch_point;
    fcn_Path_plotTraversalsXY(data,fig_num);
   
    legend('Original traversal with pinch point', 'Traversal with no pinch');
end % Ends the flag_do_plot if statement

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

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
