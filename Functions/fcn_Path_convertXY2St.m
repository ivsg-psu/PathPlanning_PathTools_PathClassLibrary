function St_points = fcn_Path_convertXY2St(referencePath,XY_points, varargin)
%% fcn_Path_convertXY2St
% Given a referencePath (N x 2) and a set of XY_points (N x 2), returns the
% Station-transverse coordinates that represent the XY points. Note: this
% function is essentially a streamlined version of
% fcn_Path_snapPointToPathViaVectors. For points snapped at verticies, the
% results may be a complex number where the real portion of the number is
% the distance along the flag_rounding_type projection vector, and the
% imaginary portion is the complex portion. The angle of the point relative
% to the vertex can be found using the complex vector angle.
% 
% FORMAT: 
%
%    St_points = fcn_Path_convertXY2St(referencePath,XY_points,
%    (flag_rounding_type), (fig_num));
%
% INPUTS:
%
%      referencePath: a Nx2 or Nx3 vector of [X Y (Z)] path points, where N
%      is the number of points the points on the path, with N >= 2.
%
%      XY_points: a Mx2 vector containing the [X Y] location of the points
%      or a Mx3 vector containing the [X Y Z] location of the points to
%      convert, where M is the number of rows of different points
%
%      (OPTIONAL INPUTS)
%      flag_rounding_type: a flag to indicate which type of projection is
%      used, especially when stations are located at the end-points of
%      segments within the nearby_traversal. When stations are at the
%      end-points of segments, the normal vector is undefined as it depends
%      on whether to use the prior or subsequent segment, or some
%      combination of these.
%
%      Note that the very first point always uses projections from the
%      following segement, and the very last point always uses the prior.
%      Otherwise, the flag determines behaviors for endpoints of internal
%      segments. The options include:
%
%          flag_rounding_type = 1;  % This is the default, and indicates
%          that the orthogonal projection of an endpoint is created by the
%          PRIOR segment leading up to each station query point.
% 
%          flag_rounding_type = 2;  % This indicates that the orthogonal
%          projection of an endpoint is created by the FOLLOWING segment
%          after each station query point.
% 
%          flag_rounding_type = 3;  % This indicates that the orthogonal
%          projection, ONLY if the station query falls at the joining point
%          between two segments (e.g. is on the "joint"), then the
%          projection is created by averaging the vector projections
%          created from the PRIOR segment and FOLLOWING segment.
% 
%          flag_rounding_type = 4;  % This indicates that the orthogonal
%          projections along segments should be calculated at the midpoints
%          of each segment, and then for each station qeuary, the vector
%          projections are interpolated from their prior and subsequent
%          vectors. For the endpoints - the start and end - the vectors are
%          aligned with the endpoints in the negative and positive
%          directions, respectively.
%
%      fig_num: figure number where results are plotted
%
% OUTPUTS:
%
%      St_points: a Mx2 or Mx3 vector containing the [S t] or [S t h]
%      coordinates corresponding to each [ X Y (Z)] input point.
%
% EXAMPLES:
%      
% See the script: script_test_fcn_Path_convertXY2St
% for a full test suite.
%
% DEPENDENCIES:
%
%     fcn_Path_snapPointToPathViaVectors
%     fcn_Path_checkInputsToFunctions
%     fcn_Path_convertPathToTraversalStructure
%
% This function was written on 2023_08_26 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2023_08_26 by S. Brennan
% -- first write of the code

flag_do_debug = 0; % Flag to plot the results for debugging
flag_do_plots = 0;
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

if flag_check_inputs == 1
    % Are there the right number of inputs?
    narginchk(2,4);
    
    % Check the data input
    fcn_Path_checkInputsToFunctions(referencePath, 'path2or3D');
    
    % Check that the dimension of the point and path match
    if length(XY_points(1,:)) ~= length(referencePath(1,:))
        error('The dimension of the XY_points points, in number of columns, must match the dimension of the referencePath, in number of columns');
    end
end


% Does user want to specify the rounding type?
flag_rounding_type = 1;
if 3 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        flag_rounding_type=temp;
    end
end

% Does user want to show the plots?
if 4 == nargin
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end

if flag_do_debug
    fig_debug = 888; %#ok<*UNRCH> 
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


[~,...
    s_coordinates,...
    first_path_point_indicies,...
    second_path_point_indicies,...
    percent_along_length,...
    distances_real,...
    distances_imaginary] = ...
    fcn_Path_snapPointToPathViaVectors(XY_points, referencePath, flag_rounding_type);

closest_point_on_path = referencePath(first_path_point_indicies,:) + percent_along_length.*(referencePath(second_path_point_indicies,:) - referencePath(first_path_point_indicies,:));
St_points = [s_coordinates distances_real+distances_imaginary*1i];

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
    axis equal;
    
    % Plot the path
    plot(referencePath(:,1),referencePath(:,2),'r.-','Linewidth',5,'Markersize',20);
    
    Npath = length(referencePath(:,1));
    % Label the path points
    for ith_point = 1:Npath
        text(referencePath(ith_point,1),referencePath(ith_point,2),sprintf('%.0d',ith_point),'Color',[1 0 0],'FontSize',12,'VerticalAlignment','bottom');
    end
    
    % Plot the query points
    plot(XY_points(:,1),XY_points(:,2),'k.','MarkerSize',20);
    
    Npoints = length(XY_points(:,1));
    % Label the query points
    for ith_point = 1:Npoints
        text(XY_points(ith_point,1),XY_points(ith_point,2),sprintf('%.0d',ith_point),'Color',[0 0 0],'FontSize',12,'VerticalAlignment','bottom');
    end
    
    
    
    % Plot the connecting vectors
    for ith_point = 1:Npoints
        connecting_vectors = XY_points(ith_point,:) - closest_point_on_path(ith_point,:);
        quiver(closest_point_on_path(ith_point,1),closest_point_on_path(ith_point,2),...
            connecting_vectors(1,1),connecting_vectors(1,2),...
            0,'Color',[0.5 0.5 0.5]);
        
        % Label the points with distances
        % labelpoint = (closest_path_points(ith_point,:) + XY_points(ith_point,:))/2;
        labelpoint = XY_points(ith_point,:);
        
        if distances_imaginary(ith_point)==0
            St_string = sprintf('St coordinates: [%.2f %.2f]',...
                St_points(ith_point,1), St_points(ith_point,2));
        else
            St_string = sprintf('St coordinates: [%.2f %.2f+%.2fi]',...
                St_points(ith_point,1), real(St_points(ith_point,2)), imag(St_points(ith_point,2)));
        end
        
        text(labelpoint(1,1),labelpoint(1,2),St_string,'VerticalAlignment','top');
        
        
        %         % Connect closest point on path to query point
        %         closest_point = closest_path_points(ith_point,:);
        %
        %         plot(...
        %             [XY_points(ith_point,1) closest_point(1,1)],...
        %             [XY_points(ith_point,2) closest_point(1,2)],'g-','Linewidth',2);
        %
        %         % Plot the closest point on path
        %         plot(closest_point(1,1),closest_point(1,2),'g.','Markersize',20);
        %         text(closest_point(1,1),closest_point(1,2),'Snap Point on Path');
    end
    
    % Make axis slightly larger
    temp = axis;
    axis_range_x = temp(2)-temp(1);
    axis_range_y = temp(4)-temp(3);
    percent_larger = 0.3;
    axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);

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
