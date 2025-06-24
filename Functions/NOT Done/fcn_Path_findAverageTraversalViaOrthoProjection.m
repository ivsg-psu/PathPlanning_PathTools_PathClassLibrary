function [traversal_average, closestXs, closestYs, closestDistances] = ...
    fcn_Path_findAverageTraversalViaOrthoProjection(data,varargin)
% fcn_Path_findAverageTraversalViaOrthoProjection
% finds the average traversal of several traversals by taking a reference
% traversal or, if of a referemce traversal is not given, it uses as a
% reference the traversal with longest number of points. 
%
% As additional outputs, for each point in the traversal, this function
% also finds the intersection point in other traversals via orthogonal
% projection, and saves these points into arrays to denote the closest X
% and Y coordinates, and the distances. These X, Y, and distance arrays
% have M columns, one for each traversal, and N rows, one for each station
% in the reference traversal.
%
% FORMAT: 
%
%      [traversal_average, closestXs, closestYs, closestDistances] = ...
%      fcn_Path_findAverageTraversalViaOrthoProjection(...
%            data,
%            (reference_traversal),...
%            (num_iterations),
%            (weight_for_averaging),
%            (fig_num));
%
% INPUTS:
%
%      data: a traversals type data structure, namely a structure
%      containing a cell array of traversals, each with subfields of X, Y,
%      etc. in the following form
%           data.traversal{i_path}.X
%      Note that i_path denotes an index into a different traversal. Each
%      traversal will be compared separately. It is assumed there are M
%      traversals, with M >=1.
%
%      (OPTIONAL INPUTS)
%
%      reference_traversal: the traversal that is being used for comparison
%      for all the other traversals. If empty, it uses the longest
%      traversal (in number of stations) from the trajectories in the data
%      structure. 
%
%      num_iterations: an integer to specify how many iterations to perform
%      in the path averaging, where the solution of the first reference
%      traversal serves as the reference traversal for the second
%      iteration's average, etc. Usually, convergence occurs within 3 to 10
%      iterations (the default is 40). Note that, within the function, a debug
%      flag can be set to plot and analyze convergence.
%
%      weight_for_averaging: a value between 0 and 1 that keeps the prior
%      solution when a new solution is found, e.g. a value of 0.8 keeps 80%
%      of the prior average and 20% of the new average, during each
%      iteration. The effect is to low-pass filter the averaging process to
%      prevent numerical bouncing that often occurs (default is 0.8).
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%      traversal_average: the resulting traversal representing the average
%      path. Usually, it is decimated at an even spacing which is currently
%      set to 1 meter.
%
%      closestXs:  a N x M vector containing the [X] location of
%      the nearest points at the N average stations projected orthogonally
%      to the M trajectories within all_traversals
%
%      closestYs:  a N x M vector containing the [Y] location of
%      the nearest points at the N average stations projected orthogonally
%      to the M trajectories within all_traversals
%
%      closestDistancess:  a N x M vector containing the distance of
%      the nearest points at the N average stations projected orthogonally
%      to the M trajectories within all_traversals. Note that positive
%      distances are those whose cross product from the
%      reference_trajectory to the intersection is positive, negative
%      distances are in the opposite direction
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_Path_findTraversalWithMostData
%      fcn_Path_findOrthogonalTraversalVectorsAtStations
%      fcn_Path_findOrthoScatterFromTraversalToTraversals
%      fcn_Path_cleanPathFromForwardBackwardJogs
%      fcn_Path_plotTraversalsXY
%      fcn_Path_newTraversalByStationResampling
%      fcn_Path_convertPathToTraversalStructure
%
% EXAMPLES:
%      
%     See the script: 
%     script_test_fcn_Path_findAverageTraversalViaOrthoProjection
%     script_test_fcn_Path_findAverageTraversal
%     for a full test suite.
%
% This function was written on 2020_11_15 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%     
% 2020_11_15  - S. Brennan
% -- wrote the code originally - lots of bugs
% 2020_12_25  - S. Brennan
% -- added more comments
% 2021_01_01  - S. Brennan
% -- fixed the errors with interpolation
% -- fixed the bug with the end point truncating toward start
% 2021_01_06  - S. Brennan
% -- added functions for input checking
% 2021_01_07  - S. Brennan
% -- deleted unused functions
% -- renamed function to reflect the traversal output, not path
% 2021_01_09:  - S. Brennan
% -- corrected terminology in comments
% -- fixed the input argument notation to be traversals
% 2021_12_27:  - S. Brennan
% -- corrected dependencies in comments
% 2022_01_03:  - S. Brennan
% -- corrected typos in comments
% -- fixed a bug where the Z value is not defined in loop
% 2022_01_06:  - S. Brennan
% -- refactored code, added weighted averaging to prevent iteration
%     bouncing
% 2022_01_10:  - S. Brennan
% -- shut off debugging commentsclose
% 2024_03_14 - S. Brennan
% -- shut off traversal_average.Station calculation as it is giving wrong
% length

% TO DO
% Need to clean up the code - lots of code "lint"

flag_do_debug = 0; % Flag to show the results for debugging
flag_do_plots = 0; % Flag to plot the final results
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
    if nargin < 1 || nargin > 6
        error('Incorrect number of input arguments')
    end

    % Check the data input
    fcn_DebugTools_checkInputsToFunctions(data, 'traversals');
        
end

% Check to see if the reference_traversal was specified?
flag_calculate_reference_traversal = 1;
if nargin>=2
    temp = varargin{1};    
    if ~isempty(temp)
        reference_traversal = temp;
        if flag_check_inputs == 1           
            % Check the reference_traversal input
            fcn_DebugTools_checkInputsToFunctions(reference_traversal, 'traversal');
        end
        flag_calculate_reference_traversal = 0;
    end
end

% Check to see if the number of iterations was specified?
num_iterations = 40;  % the default number of iteration to find the average path 
if nargin >= 3
    temp = varargin{2};
    if ~isempty(temp)
        num_iterations = temp+1;
        if num_iterations<1
            error('Number of iterations must be at least 1');
        end
    end
end

% Check to see if the number of iterations was specified?
weight = 0.8;  % Default weight
if nargin >= 4
    temp = varargin{3};
    if ~isempty(temp)
        weight = temp;
    end
end
  

% Does user want to show the plots?
if 5 == nargin
    fig_num = varargin{4};
    figure(fig_num);
    flag_do_plots = 1;
else
    if flag_do_debug
        fig = figure; 
        fig_num = fig.Number;
        flag_do_plots = 1;
    end
end

%% Main code starts here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set decimation interval and averaging type
interval = 1; % meters
flag_rounding_type = 4; % Use smooth segment-to-segment rounding. See the help on fcn_Path_FindOrthogonalHitFromPathToPath for details
Ntraversals = length(data.traversal); % the number of traversals we will be averaging


%% Check to see if we need to calculate a reference traversal 
% or if it was given. If not given, use the longest traversal to start,
% with longest measured as the one with the longest station. 
if flag_calculate_reference_traversal
    index_of_longest = fcn_Path_findTraversalWithMostData(data);
    reference_traversal = data.traversal{index_of_longest}; %initial reference path
end

% Redecimate the reference traversal at stations corresponding to
% user-chosen intervals within the chosen reference traversal
reference_station_points    = (0:interval:reference_traversal.Station(end))';
redecimated_reference_traversal = ...
    fcn_Path_newTraversalByStationResampling(reference_traversal, reference_station_points);

% %% Indicate the starting reference stations
% % Do this by finding the length of each traversal, then finding the
% % standard deviation in length. Add this to the longest traversal just as a
% % safety amount, as the averaging process sometimes moves around a LOT>
% last_stations = zeros(Ntraversals,1);
% for ith_traversal = 1:Ntraversals
%     last_stations(ith_traversal,1) = data.traversal{ith_traversal}.Station(end);
% end
% std_end_station = std(last_stations);
% 
% extra_long_station = redecimated_reference_traversal.Station(end) + std_end_station;
% longer_reference_station_points    = (0:interval:extra_long_station)';
% 
% %% Set up all the traversals to be as long as these reference station points, plus a bit just in case
% for ith_traversal = 1:Ntraversals
%     current_X = data.traversal{ith_traversal}.X;
%     current_Y = data.traversal{ith_traversal}.Y;
%     current_Station = data.traversal{ith_traversal}.Station;
% 
%     
%     interp_X       = interp1(current_Station,current_X,longer_reference_station_points,'linear','extrap');
%     interp_Y       = interp1(current_Station,current_Y,longer_reference_station_points,'linear','extrap');
%     resampled_traversal = fcn_Path_convertPathToTraversalStructure([interp_X, interp_Y]);   
%     resampled_data.traversal{ith_traversal} = resampled_traversal;
% end
% 
% if flag_do_debug
%     fcn_Path_plotTraversalsXY(data,22888);
%     hold on;
%     fcn_Path_plotTraversalsXY(resampled_data,22888);
% end    
% resampled_data = data;

%% Define search radius in orthogonal direction based on standard deviation
std_deviation = fcn_Path_calcSingleTraversalStandardDeviation(redecimated_reference_traversal);
if std_deviation~=0
    search_radius = 5*std_deviation;
else
    search_radius = reference_traversal.Station(end);
end


%% Perform iterations to seek an average path
% * Project from the reference path to nearby trajectories to find projections
% * Average projections to find new average path. 
% * Clean up average and store results for next iteration or exit

% Initialize variables to hold results
average_traversals{num_iterations} = [];
total_error{num_iterations} = [];
mean_error = nan(num_iterations, 1);
iteration_error_X{num_iterations}  = 0;
iteration_error_Y{num_iterations}  = 0;

% Set iteration 1 values
average_traversals{1} = redecimated_reference_traversal;
total_error{1} = inf*redecimated_reference_traversal.X;

% Start iterations - at each iteration, find orthogonal projections,
% average results, then recalculate the reference traversal
for ith_iteration =2:num_iterations 

    % Show user what we are doing?        
    if flag_do_debug
        fprintf(1,'Averaging paths via iteration: %.0d / %.0d \n',ith_iteration,num_iterations);
        fprintf(1,'\tInitial maximum station: %.2f  \n',redecimated_reference_traversal.Station(end));    
        fprintf(1,'\tOld search radius: %.2f  \n',search_radius);
    end
    
    path_last  = redecimated_reference_traversal;   % Saves the prior path - used for error calculations later     
    
    %% For each traversal, project from reference orthogonally
    % For debugging: [closestXs, closestYs, closestDistances] = fcn_Path_findOrthoScatterFromTraversalToTraversals(reference_station_points, reference_traversal, data, flag_rounding_type,search_radius, 33304);
    [closestXs, closestYs, closestDistances] = fcn_Path_findOrthoScatterFromTraversalToTraversals(reference_station_points, redecimated_reference_traversal, data, flag_rounding_type,search_radius);
    
    %% Find and remove outliers due to distance jumps
    % Find the standard deviation in the distances
    sigma_distances = nanstd(closestDistances,0,'all');    
    search_radius = 5*sigma_distances;
    outliers = find(abs(closestDistances)>(5*sigma_distances));
    if ~isempty(outliers)
        closestXs(outliers) = NaN;
        closestYs(outliers) = NaN;
        closestDistances(outliers) = NaN;         %#ok<NASGU>
    end
    
    if flag_do_debug
        fprintf(1,'\tStandard deviation in distances: %.2f  \n',sigma_distances);
        fprintf(1,'\tNumber of outliers: %d  \n',length(outliers));
        fprintf(1,'\tNew search radius: %.2f  \n',search_radius);
    end

    
    
    %% Average distances at the projection points to generate a mean path    
    raw_path_mean = [nanmean(closestXs,2) nanmean(closestYs,2)];
    
    
    %% Call a special function to remove back-tracking behavior
    % Sometimes averaging can produce data that bounces forward and
    % backward. This function removes these
    path_mean = fcn_Path_cleanPathFromForwardBackwardJogs(raw_path_mean);        
    
    % Sometimes the path average has NaN values 
    % This will cause the interpolation step to fail.
    % This means that no paths were detected nearby the reference
    % traversal. In this case, we drop these points and the interpolation
    % will extrapolate via linear fit what the values should be
    path_mean_no_nan = path_mean(~isnan(path_mean(:,1)),:);
    path_mean_station_no_nan  = [0; cumsum(sqrt(sum(diff(path_mean_no_nan).^2,2)),'omitnan')];

    
    %% Interpolation of mean data to produce equal station intervals
    path_average_interp_X       = interp1(path_mean_station_no_nan,path_mean_no_nan(:,1),reference_station_points,'linear','extrap');
    path_average_interp_Y       = interp1(path_mean_station_no_nan,path_mean_no_nan(:,2),reference_station_points,'linear','extrap');
        
    %% Do weighted average of X and Y values
    path_new_X = weight*path_last.X + (1-weight)*path_average_interp_X;
    path_new_Y = weight*path_last.Y + (1-weight)*path_average_interp_Y;
    
    
    %% Convert path type back into traversal and save result
    redecimated_reference_traversal = fcn_Path_convertPathToTraversalStructure([path_new_X, path_new_Y]);  
    redecimated_reference_traversal.Station = reference_station_points; % The calculation of station is a bit off in the conversion, so fix it here.
    average_traversals{ith_iteration} = redecimated_reference_traversal;

    if flag_do_debug
        fprintf(1,'\tNew maximum station: %.2f  \n',redecimated_reference_traversal.Station(end));
    end
    
    %% Update error calculations for iterations 2 and onward    
    iteration_error_X{ith_iteration} = path_last.X - redecimated_reference_traversal.X;
    iteration_error_Y{ith_iteration} = path_last.Y - redecimated_reference_traversal.Y;

    % Calculate the total error
    X_error = iteration_error_X{ith_iteration};
    Y_error = iteration_error_Y{ith_iteration};
    total_error{ith_iteration} = (X_error.^2 + Y_error.^2).^0.5;
    mean_error(ith_iteration,1) = mean(total_error{ith_iteration});
    
    % Show results?
    if flag_do_debug
        fprintf(1,'\t Mean change from iteration %.0d to %.0d: %.3f \n',ith_iteration-1, ith_iteration,mean_error(ith_iteration,1));
    end

    
end

% Use final average path to define "true" s-coordinates of the original trajectories, using projection
traversal_average = fcn_Path_convertPathToTraversalStructure([path_mean_no_nan(:,1), path_mean_no_nan(:,2)]);  
% traversal_average.Station = reference_station_points; % The calculation of station is a bit off in the conversion, so fix it here.
 
% Calculate final results
[closestXs, closestYs, closestDistances] = ...
    fcn_Path_findOrthoScatterFromTraversalToTraversals(...
    traversal_average.Station, redecimated_reference_traversal, data,...
    flag_rounding_type,search_radius);


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
    
    % plot the final XY result        
    figure(fig_num); 
    clf;
    
    fcn_Path_plotTraversalsXY(data,fig_num);
    hold on;
    plot(traversal_average.X,traversal_average.Y,'b.-','Linewidth',4,'Markersize',20);
    title('original paths and final average path');
    xlabel('x [m]');
    ylabel('y [m]');
 
    
    %     % Plot the yaw results
    %     figure(111);
    %     clf
    %     hold on
    %     plot(reference_traversal.Station, reference_traversal.Yaw,'r-' ,'LineWidth',2)
    %     plot(reference_traversal.Station,fnc_yaw(reference_traversal.X,reference_traversal.Y),'b-' ,'LineWidth',4)
    %     for i_path= 1:length(data.traversal)
    %         plot(data.traversal{i_path}.Station,data.traversal{i_path}.Yaw)
    %     end
    %     legend('result of interpolation', 'result of final calculation')
    %     title('station vs yaw')
    %     xlabel('station [m]')
    %     ylabel('yaw [degree]')
    
    if flag_do_debug
        % Plot the path convergence
        figure(2255);
        clf;
        hold on;
        xlabel('Index')
        ylabel('Position change between iterations [m]')
        for ith_iteration = 1:length(iteration_error_X)
            title(sprintf('Iteration %.0d of %.0d',ith_iteration,length(iteration_error_X)));
            plot(average_traversals{ith_iteration}.X,average_traversals{ith_iteration}.Y,'-'); 
            drawnow;
            pause(0.1);
        end
        
        % Plot the error convergence
        figure(2233);
        clf;
        hold on;
        xlabel('Index')
        ylabel('Position change between iterations [m]')
        for ith_iteration = 1:length(iteration_error_X)
            title(sprintf('Iteration %.0d of %.0d',ith_iteration,length(iteration_error_X)));
            plot(total_error{ith_iteration});            
            pause(0.1);
        end
        
        
        % Plot the error convergence
        figure(2244);
        clf;
        
        semilogy(1:length(mean_error),mean_error,'.-','Markersize',30);
        xlabel('Iteration')
        ylabel('Mean change in distance')
        grid on
        
    end
    
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end
end % End of function


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
