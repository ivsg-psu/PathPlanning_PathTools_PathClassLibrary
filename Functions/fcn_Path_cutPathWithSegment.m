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
%      fig_num: figure number where results are plotted
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
%     fcn_Path_checkInputsToFunctions
%     fcn_Path_findProjectionHitOntoPath
%
% This function was written on 2023_09_26 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2023_09_26 by S. Brennan
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
    narginchk(2,3);

    % Check the data input
    fcn_Path_checkInputsToFunctions(pathToCut, 'path2or3D');

    % Check that the dimension of the point and path match
    if (length(cutting_segment(1,:)) ~= 2) || (length(cutting_segment(:,1)) ~= 2)
        error('The cutting_segment definition must be a 2 x 2 matrix, formatted as: [start_x start_y, end_x, end_y]');
    end
end


% Does user want to show the plots?
if 3 == nargin
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

% Set the search type
flag_search_type = 0; % 0 searches for exact crossing
% flag_search_type = 1; % 1 searches for any crossing, in any direction

% Initialize the outputs
cut_path_before = [];
cut_path_after  = [];


% Look for a crossing of the pathToCut with the cutting_segment
[distance,location_of_crossing] = ...
    fcn_Path_findProjectionHitOntoPath(pathToCut,...
    cutting_segment(1,1:2),cutting_segment(2,1:2),...
    (flag_search_type));

if isnan(distance) || distance<0 
    % Return empty values - do nothing

else
    [~,s_coordinate_of_cut] = ...
        fcn_Path_snapPointToPathViaVectors(location_of_crossing, pathToCut);

    traversalToCut = fcn_Path_convertPathToTraversalStructure(pathToCut);

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
            location_of_crossing(1,3) = height_at_cut;

        end
        cut_path_after  = [location_of_crossing; cut_path_after];
        cut_path_before = [cut_path_before; location_of_crossing];
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
    plot(cutting_segment(:,1),cutting_segment(:,2),'b.-','Linewidth',5,'MarkerSize',20);

    % Plot the cut_path_before
    if ~isempty(cut_path_before)
        n_plots = n_plots+1;
        legend_strings{n_plots} = 'cut_path_before';
        plot(cut_path_before(:,1),cut_path_before(:,2),'r-','Linewidth',2,'MarkerSize',20);
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
