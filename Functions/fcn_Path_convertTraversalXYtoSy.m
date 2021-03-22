function [closestXs, closestYs, closestDistances] = ...
    fcn_Path_convertTraversalXYtoSy(...
    reference_station_points, reference_traversal, ...
    all_traversals, varargin)

% fcn_Path_convertTraversalXYtoSy
% Given a reference traversal and a set of stations along that traversal,
% finds the location on each nearby traversal that is closest to the central
% traveral at each station point. Closest is defined via an orthogonal
% projection (or modifications of orthogonal projections) from the central
% traversal outward toward nearby traversals. The results are then
% presented as an array of lateral offsets.
%
% NOTE: this is just a reformulation of the function:
% fcn_Path_findOrthoScatterFromTraversalToTraversals
% The main difference is the final plotting result, but otherwise it uses
% the same functionality.
%
% FORMAT:
%
%      [closestXs, closestYs, closestDistances] = ...
%        fcn_Path_convertTraversalXYtoSy(...
%        reference_station_points, reference_traversal, all_traversals,...
%        (flag_rounding_type),(search_radius),(fig_num));
%
% INPUTS:
%
%      reference_station_points: an N x 1 vector containing the stations on
%      the reference_traversal where the orthogonal projections should take
%      place
%
%      reference_traversal: a traversal structure that specifies the
%      traversal where projections to other traversals are taking place.
%
%      all_traversals: a traversals type data structure containing a cell
%      array of traversal structures that specifies the traversals where
%      intersections from the reference_traversal would hit.
%
%      (OPTIONAL INPUTS)
%      flag_rounding_type: a flag to indicate which type of projection is
%      used, especially when stations are located at the end-points of
%      segments within the nearby_traversal. Note that the very first point
%      always uses projections from the following segement, and the very
%      last point always uses the prior. The flag determines behaviors for
%      endpoints of internal segments. The options include:
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
%          vectors. Note that this can generate non-orthogonal results.
%
%      search_radius: how far in station distance for the search points to
%      look (default is 3 times the differences in station lengths within
%      all_traversals )
%
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      closestXs:  a NxM vector containing the [X] location of
%      the nearest points at the N stations projected orthogonally to the M
%      trajectories with all_traversals
%
%      closestYs:  a NxM vector containing the [Y] location of
%      the nearest points at the N stations projected orthogonally to the M
%      trajectories with all_traversals
%
%      closestDistancess:  a NxM vector containing the distance of
%      the nearest points at the N stations projected orthogonally to the M
%      trajectories with all_traversals. Note that positive distances are
%      those whose cross product from the reference_trajectory to the
%      intersection is positive, negative distances are in the opposite
%      direction
%
%
% EXAMPLES:
%
% See the script: script_test_fcn_Path_convertTraversalXYtoSy
% for a full test suite.
%
% DEPENDENCIES:
%
%      fcn_Path_FindOrthogonalHitFromTraversalToTraversals
%      fcn_Path_findOrthogonalTraversalVectorsAtStations
%
% This function was written on 2021_03_21 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
%      2021_03_21:
%      -- first write of the code moving this functionality out of
%      fcn_Path_findOrthoScatterFromTraversalToTraversals.m

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

if flag_check_inputs == 1
    % Are there the right number of inputs?
    if nargin < 3 || nargin > 6
        error('Incorrect number of input arguments')
    end
    
    % Check the station input
    fcn_Path_checkInputsToFunctions(reference_station_points, 'station');
    
    % Check the reference_traversal input
    fcn_Path_checkInputsToFunctions(reference_traversal, 'traversal');
    
    % Check the all_traversals input
    fcn_Path_checkInputsToFunctions(all_traversals, 'traversals');
    
end

% Fil in default values
Nstations = length(reference_station_points(:,1));
Ntraversals = length(all_traversals.traversal);
Station_central = reference_traversal.Station;

flag_rounding_type = 1;
% Does the user want to give a rounding type?
if 4 <= nargin
    flag_rounding_type = varargin{1};
end

% Define search radius?
search_radius = max(Station_central)*3;
if 5 <= nargin
    search_radius = varargin{2};
end


% Does user want to show the plots?
if 6 == nargin
    fig_num = varargin{3};
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

%% Initialize matrices
% These are used to save the current reference traversal (so we can
% check changes after this is updated) and initialize arrays for the loop

closestXs = zeros(Nstations,Ntraversals);
closestYs = zeros(Nstations,Ntraversals);
closestDistances = zeros(Nstations,Ntraversals);


%% For each traversal, project from reference orthogonally
% Search nearest points from reference path to each of the other
% traversals, saving the results as "closest" values

for ith_traversal = 1:Ntraversals
    nearby_traversal = all_traversals.traversal{ith_traversal};
    
    [closest_path_points,closest_distances] = ...
        fcn_Path_findOrthogonalHitFromTraversalToTraversals(...
        reference_station_points,...
        reference_traversal,...
        nearby_traversal,...
        flag_rounding_type,...
        search_radius);
    
    % Save final results as closest points
    closestXs(:,ith_traversal) = closest_path_points(:,1);
    closestYs(:,ith_traversal)  = closest_path_points(:,2);
    closestDistances(:,ith_traversal) = closest_distances;
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
    
    xlabel('Station [m]');
    ylabel('y-offset [m]');
    % axis equal;
    
    % Plot the results
    for ith_traversal = 1:length(closestDistances(1,:))        
        plot(reference_station_points, closestDistances(:, ith_traversal),'.-','Markersize',15,'Linewidth',3);
    end
    
    % For debugging    
    if 1== flag_do_debug  % The following plots the XY results
        figure(fig_num+1);
        clf;
        hold on;
        grid on;
        
        xlabel('X [m]');
        ylabel('Y [m]');
        axis equal;
        
        % Plot the central traversal
        plot(reference_traversal.X,reference_traversal.Y,'k.-','Linewidth',3,'Markersize',25);
        
        % Plot the paths
        fcn_Path_plotTraversalsXY(all_traversals,fig_num);
        
        % Plot the station points
        % Find the unit normal vectors at each of the station points
        [unit_normal_vector_start, unit_normal_vector_end] = ...
            fcn_Path_findOrthogonalTraversalVectorsAtStations(...
            reference_station_points,reference_traversal,flag_rounding_type);
        
        plot(unit_normal_vector_start(:,1),unit_normal_vector_start(:,2),'r.','Markersize',35);
        
        
        % Plot the results
        for i_point = 1:length(closestXs(:,1))
            plot(closestXs(i_point,:),closestYs(i_point,:),'bo-','Markersize',15,'Linewidth',3);
        end
        
        % Show the unit vectors
        normal_unit_vectors_at_stations = ...
            unit_normal_vector_end - unit_normal_vector_start;
        quiver(unit_normal_vector_start(:,1),unit_normal_vector_start(:,2),...
            normal_unit_vectors_at_stations(:,1),normal_unit_vectors_at_stations(:,2),0,'g','Linewidth',3);  % The 0 term is to prevent scaling
    end
    
end % Ends the flag_do_plot if statement

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends the function

