function XY_points = fcn_Path_convertSt2XY(referencePath,St_input_points, varargin)
%% fcn_Path_convertSt2XY
% Given a referencePath (N x 2) and a set of XY_points (N x 2), returns the
% XY coordinates that represent the St points. 
% 
% FORMAT: 
%
%    XY_points = fcn_Path_convertSt2XY(referencePath,St_input_points,
%    (flag_rounding_type), (fig_num));
%
% INPUTS:
%
%      referencePath: a Nx2 or Nx3 vector of [X Y (Z)] path points, where N
%      is the number of points the points on the path, with N >= 2.
%
%      St_input_points: a Mx2 vector containing the [S t] location of the
%      points or a Mx3 vector containing the [S t h] location of the points
%      to convert, where M is the number of rows of different points
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
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%      XY_points: a Mx2 or Mx3 vector containing the [X Y] or [X Y Z]
%      coordinates corresponding to each [ S t (h)] input point.
%
% EXAMPLES:
%      
% See the script: script_test_fcn_Path_convertSt2XY
% for a full test suite.
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
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
    fcn_DebugTools_checkInputsToFunctions(referencePath, 'path2or3D');
    
    % Check that the dimension of the point and path match
    if length(St_input_points(1,:)) ~= length(referencePath(1,:))
        error('The dimension of the St_points, in number of columns, must match the dimension of the referencePath, in number of columns');
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

central_traversal = fcn_Path_convertPathToTraversalStructure(referencePath);

modified_stations = min(St_input_points(:,1),central_traversal.Station(end));
modified_stations = max(modified_stations,0);

[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalTraversalVectorsAtStations(modified_stations,central_traversal, flag_rounding_type);

unit_orthogonal_vectors = unit_normal_vector_end - unit_normal_vector_start;

point_dimension = length(St_input_points(1,:));
if point_dimension == 3
    rotation_matrix = [0 1 0; -1 0 0; 0 0 1];
else
    rotation_matrix = [0 1; -1 0];
end

unit_tangent_vectors = unit_orthogonal_vectors*rotation_matrix;

XY_points = unit_normal_vector_start + unit_orthogonal_vectors.*real(St_input_points(:,2)) + unit_tangent_vectors.*imag(St_input_points(:,2));


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
    
    % Plot the query stations
    plot(unit_normal_vector_start(:,1),unit_normal_vector_start(:,2),'k.','MarkerSize',20);
    
    
    % Plot the solution
    plot(XY_points(:,1),XY_points(:,2),'b.','MarkerSize',20);
    
    
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
