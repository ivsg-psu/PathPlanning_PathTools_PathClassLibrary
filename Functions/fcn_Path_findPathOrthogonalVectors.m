function [normal_unit_vectors_at_midpoints, normal_unit_vectors_at_joints] = ...
    fcn_Path_findPathOrthogonalVectors(path, varargin)

% fcn_Path_findPathOrthogonalVectors
% Given a central traversal, finds the orthogonal vectors of that
% traversal. Includes special flag input to define the meaning of
% orthogonal vectors at verticies and ends, as these may not match the
% midpoint vectors. If the input path is in 3D, the orthogonal vectors are
% not clearly defined and the projection is ONLY done in XY space, leaving
% the z-component of the projection as zero.
%
% FORMAT:
%
%      [normal_unit_vectors_at_midpoints, unit_joint_vectors] = ...
%        fcn_Path_findPathOrthogonalVectors(path,...
%        (flag_rounding_type),(fig_num));
%
% INPUTS:
%
%      path: a path type input (N x 2 or N x 3) that specifies the path
%      geometry.
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
%          vectors.
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%      normal_unit_vectors_at_midpoints: the unit Nx2 vector that is
%      the normal projection of the segements at the midpoints.
%
%      normal_unit_vectors_at_joints: the unit (N+1)x2 vector that is
%      the normal projection of the segements at the joints. There are N+1
%      joints because, for a N segment path, there are N-1 joints and the
%      starting and ending points are also joints.
%
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
% See the script: script_test_fcn_Path_findPathOrthogonalVectors
% for a full test suite.
%
% This function was written on 2023_08_27 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2023_08_27 - S. Brennan
% -- first write of the code via modification from
% fcn_Path_ findOrthogonalTraversalVectorsAtStations
% 2024_03_14 - S. Brennan
% -- fixed bug where snap breaks if path is passed in as a 3D vector
% 2025_06_23 - S. Brennan
% -- Updated debugging and input checks

% TO-DO
% Define search radius - need to let user define this as an input!

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 3; % The largest Number of argument inputs to the function
flag_max_speed = 0;
if (nargin==MAX_NARGIN && isequal(varargin{end},-1))
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
        narginchk(1,MAX_NARGIN);

        % Check the path input
        fcn_DebugTools_checkInputsToFunctions(path, 'path2or3D');
    end
end

% Does user want to specify the rounding type?
flag_rounding_type = 1;
if 2 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        flag_rounding_type=temp;
    end
end

% Does user want to show the plots?
flag_do_plots = 0; % Default is to NOT show plots
if (0==flag_max_speed) && (MAX_NARGIN == nargin)
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        fig_num = temp;
        figure(fig_num);
        flag_do_plots = 1;
    end
else
    if flag_do_debug
        fig_debug = 818181; %#ok<NASGU,*UNRCH>
        flag_do_plots = 1;
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

%% Find the midpoint tangent vectors for all segments in the central path
% These are used to calculate the midpoint vectors at each segment
% Convert these into vector form
tangent_vectors_at_midpoints = diff(path);
magnitudes = sum(tangent_vectors_at_midpoints.^2,2).^0.5;
tangent_unit_vectors_at_midpoints = tangent_vectors_at_midpoints./magnitudes;

if length(tangent_unit_vectors_at_midpoints(1,:))==2
    normal_unit_vectors_at_midpoints= tangent_unit_vectors_at_midpoints*[0 1; -1 0];
elseif  length(tangent_unit_vectors_at_midpoints(1,:))==3
    normal_unit_vectors_at_midpoints= tangent_unit_vectors_at_midpoints*[0 1 0; -1 0 0; 0 0 1];
else
    error('A 2D or 3D vector is expected');
end

%% Find the joint tangent vectors for all joints in the central path
% The joint projection vector for each joint will depend on the input
% options. For options where the orthogonal projection only uses the prior
% or subsequent values, we can use the previously calculated vectors. For
% options where averaging is occuring, we use averages of vectors.
%
% Note that the start and end points are considered "joints" also, so if
% there are N segments, then there are N+1 joints.
%
% Here's the rounding options:
%      flag_rounding_type = 1;  % This is the default, and indicates that
%      the orthogonal projection at a joint is created by the PRIOR
%      segment.
%
%      flag_rounding_type = 2;  % This indicates that the orthogonal
%      projection at a joint is created by the FOLLOWING segment.
%
%      flag_rounding_type = 3;  % This indicates that the orthogonal
%      projection, ONLY at a joint, is created by averaging both the
%      PRIOR segment and FOLLOWING segment.
%
%      flag_rounding_type = 4;  % This indicates that the orthogonal
%      projections at the joint is created by averaging. Additionally,
%      along the entire segment, vectors should be calculated by
%      interpolating between the midpoint and the joint vectors.

switch flag_rounding_type
    case 1  % Default - use orthogonal projection via prior segment
        normal_unit_vectors_at_joints = [normal_unit_vectors_at_midpoints(1,:); normal_unit_vectors_at_midpoints];
    case 2  % Use orthogonal projection of subsequent segment
        normal_unit_vectors_at_joints = [normal_unit_vectors_at_midpoints; normal_unit_vectors_at_midpoints(end,:)];
    case 3  % Average the two segments but only at the endpoints
        average_vectors = (normal_unit_vectors_at_midpoints(1:end-1,:)+normal_unit_vectors_at_midpoints(2:end,:))/2;
        magnitudes = sum(average_vectors.^2,2).^0.5;
        unit_average_vectors = average_vectors./magnitudes;
        normal_unit_vectors_at_joints = [normal_unit_vectors_at_midpoints(1,:); unit_average_vectors; normal_unit_vectors_at_midpoints(end,:)];
    case 4 % Average the two segments along all points
        average_vectors = (normal_unit_vectors_at_midpoints(1:end-1,:)+normal_unit_vectors_at_midpoints(2:end,:))/2;
        magnitudes = sum(average_vectors.^2,2).^0.5;
        unit_average_vectors = average_vectors./magnitudes;
        normal_unit_vectors_at_joints = [normal_unit_vectors_at_midpoints(1,:)*[0 1; -1 0]; unit_average_vectors; normal_unit_vectors_at_midpoints(end,:)*[0 -1; 1 0]];
    otherwise
        error('Unrecognized method in flag_rounding_type.');
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
    hold on;
    grid on;

    plot_color = [0 0 1];
    line_width = 3;

    % Plot the central traversal in blue, with black dots
    sizeOfMarkers = 40;
    plot(path(:,1),path(:,2), '.-','Color',plot_color,'Linewidth',line_width,'Markersize',sizeOfMarkers);
    plot(path(:,1),path(:,2), '.','Color',[0 1 0],'Linewidth',line_width,'Markersize',round(sizeOfMarkers*0.5));


    % Show the unit vectors at the midpoints in red
    midpoints = (path(1:end-1,:)+path(2:end,:))./2;
    quiver(midpoints(:,1),midpoints(:,2),normal_unit_vectors_at_midpoints(:,1),normal_unit_vectors_at_midpoints(:,2),0,'r','Linewidth',2);

    % Show the unit vectors at the joints in green
    quiver(path(:,1),path(:,2),normal_unit_vectors_at_joints(:,1),normal_unit_vectors_at_joints(:,2),0,'g','Linewidth',2);

    % Add a legend
    legend('Path','Path joints','Midpoint normal unit vectors','Joint normal unit vectors','Location','southwest');

    % Make axis slightly larger
    temp = axis;
    axis_range_x = temp(2)-temp(1);
    axis_range_y = temp(4)-temp(3);
    percent_larger = 0.3;
    axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);


end % Ends the flag_do_plots if statement

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

