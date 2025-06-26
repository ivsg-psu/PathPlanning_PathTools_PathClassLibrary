function [cut_path_before, cut_path_after] = fcn_Path_cutPathWithSegment(pathToCut, cutting_segment, varargin)
%% fcn_Path_cutPathWithSegment
% Given a pathToCut, an (N x 2) or (N x 3) matrix, and a path segment
% defined by two end-points, (N x 2) in form [start_x start_y, end_x,
% end_y], cuts the path at the point where it crosses the segment. If the
% cut does not pass through a point on the path, a point is inserted at the
% cut point so that the cut is exactly at the segment crossing. If there
% are more than one crossing, then the first crossing is used.
%
% FORMAT:
%
%    [cut_path_before, cut_path_after] = fcn_Path_cutPathWithSegment(pathToCut,cutting_segment,(fig_num));
%
% INPUTS:
%
%      pathToCut: a Nx2 or Nx3 vector of [X Y (Z)] path points, where N
%      is the number of points the points on the path, with N >= 2.
%
%      cutting_segment: a Mx2 vector containing the [X Y] location of the
%      start and end points of the cutting segment. Data should be in form
%      [start_x start_y, end_x, end_y]
%
%      (OPTIONAL INPUTS)
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%      cut_path_before: a Mx2 or Mx3 vector containing the points "before"
%      the segment, where before is the point of the path with station
%      coordinates before the cut.
%
%      cut_path_after: a Mx2 or Mx3 vector containing the points "after"
%      the segment, where after is the point of the path with station
%      coordinates after the cut.
%
% EXAMPLES:
%
% See the script: script_test_fcn_Path_cutPathWithSegment
% for a full test suite.
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%     fcn_Path_findProjectionHitOntoPath
%
% This function was written on 2023_09_26 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2023_09_26 by S. Brennan
% -- first write of the code
% 2025_06_23 - S. Brennan
% -- Updated debugging and input checks

% TO-DO
% (none)

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
        narginchk(2,MAX_NARGIN);

        % Check the data input
        fcn_DebugTools_checkInputsToFunctions(pathToCut, 'path2or3D');

        % Check that the dimension of the point and path match
        if (length(cutting_segment(1,:)) ~= 2) || (length(cutting_segment(:,1)) ~= 2)
            error('The cutting_segment definition must be a 2 x 2 matrix, formatted as: [start_x start_y, end_x, end_y]');
        end
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
end

if flag_do_debug
    fig_debug = 42545; %#ok<NASGU>
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

% Set the search type
% flag_search_type = 0; % 0 searches for exact crossing
% flag_search_type = 1; % 1 searches for any crossing, in any direction
flag_search_type = 2; % 2 searches for all crossings

% Initialize the outputs
cut_path_before = [];
cut_path_after  = [];


% Look for a crossing of the pathToCut with the cutting_segment
[distances,locations_of_crossing] = ...
    fcn_Path_findProjectionHitOntoPath(pathToCut,...
    cutting_segment(1,1:2),cutting_segment(2,1:2),...
    (flag_search_type), -1);

if all(isnan(distances)) || all(distances<0)
    % Return empty values - do nothing

else
    % For all the points that were hit, find the one with lowest
    % s-coordinate
    s_coordinate_of_cut = inf;
    correct_crossing = 0;
    for ith_crossing = 1:length(locations_of_crossing(:,1))
        location_of_crossing = locations_of_crossing(ith_crossing,:);
        [~,s_coordinate_of_crossing] = ...
            fcn_Path_snapPointToPathViaVectors(location_of_crossing, pathToCut, -1);
        if s_coordinate_of_crossing <= s_coordinate_of_cut
            s_coordinate_of_cut = s_coordinate_of_crossing;
            correct_crossing = ith_crossing;
        end
    end

    % Prep path for cutting
    traversalToCut = fcn_Path_convertPathToTraversalStructure(pathToCut, -1);
    trimmed_after_indicies  = traversalToCut.Station >=s_coordinate_of_cut;
    trimmed_before_indicies = traversalToCut.Station <=s_coordinate_of_cut;
    
    cut_path_before = pathToCut(trimmed_before_indicies,:);
    cut_path_after  = pathToCut(trimmed_after_indicies,:);

    % Do we need to insert a point?
    if ~any(traversalToCut.Station == s_coordinate_of_cut)
        if length(pathToCut(1,:))==3
            % Must find height at cut. Use interpolation to do this.

            index_before_cut = find(traversalToCut.Station<s_coordinate_of_cut,1,'last');
            station_before_cut = traversalToCut.Station(index_before_cut,1);
            height_before_cut  = pathToCut(index_before_cut,3);
            
            index_after_cut = find(traversalToCut.Station>s_coordinate_of_cut,1,'first');
            station_after_cut = traversalToCut.Station(index_after_cut,1);
            height_after_cut  = pathToCut(index_after_cut,3);

            stations = [station_before_cut; station_after_cut];
            heights  = [height_before_cut; height_after_cut];

            % Vq = interp1(X,V,Xq)
            height_at_cut = interp1(stations,heights,s_coordinate_of_cut);
            locations_of_crossing(1,3) = height_at_cut;

        end
        cut_path_after  = [locations_of_crossing(correct_crossing,:); cut_path_after];
        cut_path_before = [cut_path_before; locations_of_crossing(correct_crossing,:)];
    end
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
    axis equal;

    % Plot the pathToCut
    n_plots = 0;

    n_plots = n_plots+1;
    legend_strings{n_plots} = 'pathToCut';
    plot(pathToCut(:,1),pathToCut(:,2),'k.-','Linewidth',5,'Markersize',20);

    % Plot the cutting_segment
    n_plots = n_plots+1;
    legend_strings{n_plots} = 'cutting_segment';
    plot(cutting_segment(:,1),cutting_segment(:,2),'b.-','Linewidth',6,'MarkerSize',20);

    % Plot the cut_path_before
    if ~isempty(cut_path_before)
        n_plots = n_plots+1;
        legend_strings{n_plots} = 'cut_path_before';
        plot(cut_path_before(:,1),cut_path_before(:,2),'r-','Linewidth',4,'MarkerSize',20);
    end

    % Plot the cut_path_after
    if ~isempty(cut_path_after)
        n_plots = n_plots+1;
        legend_strings{n_plots} = 'cut_path_after';
        plot(cut_path_after(:,1),cut_path_after(:,2),'g-','Linewidth',2,'MarkerSize',20);
    end

    % h_legend = legend('pathToCut','cutting_segment','cut_path_before','cut_path_after');
    h_legend = legend(legend_strings);
    set(h_legend,'Interpreter','none');

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
