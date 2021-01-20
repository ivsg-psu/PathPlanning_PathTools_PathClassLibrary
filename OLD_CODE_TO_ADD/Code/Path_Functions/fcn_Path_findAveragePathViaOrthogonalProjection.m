function [path_average_final, closestXs, closestYs, closestDistances] = fcn_Path_findAveragePathViaOrthogonalProjection(data,varargin)
% fcn_Path_findAveragePathViaOrthogonalProjection
% finds the average of several paths by taking a reference traversal (or,
% if one is not given, using the traversal with longest number of points)
% and for each point in the traversal finding the intersection point in other
% traversals via orthogonal projection. 
%
% FORMAT: 
%
%      path_average_final = ...
%      fcn_Path_findAveragePathViaOrthogonalProjection(...
%            data,
%            (reference_traversal),...
%            (flag_calculate_reference_traversal),(num_iterations),
%            (fig_num));
%
% INPUTS:
%
%      data: a structure containing i traversal fields, each with subfields
%      of X, Y, etc. in the following form
%           data.traversal{i_path}.X
%      Note that i_path denotes an array of paths. Each path will be
%      compared separately.
%
%      (OPTIONAL INPUTS)
%
%      reference_traversal: the traversal that is being used for comparison
%      for all the other traversals. If empty, sets
%      flag_calculate_reference_traversal to 1 and uses result of this.
%
%      flag_calculate_reference_traversal: a flag to calculate a reference
%      traveral using the longest traversal (in legth of vector) from the
%      trajectories in the data structure. If set to 1, this overrides the
%      reference_traversal path given as an input.
%
%      num_iterations: an integer to specify how many iterations to perform
%      in the path averaging, where the solution of the first reference
%      traversal serves as the reference traversal for the second
%      iteration's average, etc. Usually, convergence occurs within 3
%      iterations (the default). Note that, within the function, a debug
%      flag can be set to plot and analyze convergence.
%
%      fig_num: a figure number to plot results.
%
%
% OUTPUTS:
%
%      path_average_final: the resulting traversal representing the average
%      path. Usually, it is decimated at an even spacing which is currently
%      set to 1 meter.
%
%      closestXs:  a N x M vector containing the [X] location of
%      the nearest points at the N average stations projected orthogonally
%      to the M trajectories with all_traversals
%
%      closestYs:  a N x M vector containing the [Y] location of
%      the nearest points at the N average stations projected orthogonally
%      to the M trajectories with all_traversals
%
%      closestDistancess:  a N x M vector containing the distance of
%      the nearest points at the N average stations projected orthogonally
%      to the M trajectories with all_traversals. Note that positive
%      distances are those whose cross product from the
%      reference_trajectory to the intersection is positive, negative
%      distances are in the opposite direction
%
% DEPENDENCIES:
%
%      fcn_Path_findClosestPointsFromPath
%      fcn_Path_findTraversalWithMostData
%      fcn_Path_plotPathXY
%
% EXAMPLES:
%      
%     See the script: script_test_fcn_Path_findAveragePath
%     for a full test suite.
%
% This function was written on 2020_11_15 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%     2020_11_15: 
%     -- wrote the code originally - lots of bugs
%     2020_12_25:
%     -- added more comments
%     2021_01_01
%     -- fixed the errors with interpolation
%     -- fixed the bug with the end point truncating toward start
%     2021_01_06
%     -- added functions for input checking

% TO DO
% Pull out the last function that checks jogs

flag_do_debug = 0; % Flag to show the results for debugging
flag_do_plots = 0; % Flag to plot the final results
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'Starting function: %s, in file: %s\n',st(1).name,st(1).file);
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

% Set defaults
num_iterations = 3;  % the default number of iteration to find the average path 
flag_calculate_reference_traversal = 1;

if flag_check_inputs == 1
    % Are there the right number of inputs?
    if nargin < 1 || nargin > 5
        error('Incorrect number of input arguments')
    end

    % Check the data input
    fcn_Path_checkInputsToFunctions(data, 'traversals');
        
end

if nargin>=2
    temp = varargin{1};    
    if ~isempty(temp)
        reference_traversal = temp;
        if flag_check_inputs == 1           
            % Check the reference_traversal input
            fcn_Path_checkInputsToFunctions(reference_traversal, 'traversal');
        end 
    end
end

% Check to see if the reference traversal was given?
if nargin >=3
    temp = varargin{2};
    if ~isempty(temp)
        flag_calculate_reference_traversal = temp;
    end
end

% Check to see if the number of iterations was specified?
if nargin >= 4
    temp = varargin{3};
    if ~isempty(temp)
        num_iterations = temp;
    end
end
    
% Does user want to show the plots?
if 5 == nargin
    fig_num = varargin{3};
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

%% Check to see if we need to calculate a reference traversal 
% or if it was given
if flag_calculate_reference_traversal
    index_of_longest = fcn_Path_findTraversalWithMostData(data);
    reference_traversal = data.traversal{index_of_longest}; %initial reference path
end

%% Indicate the starting reference stations and looping variables
interval = 1; % meters
Ntraversals = length(data.traversal);
flag_rounding_type = 4; % See the help on fcn_Path_FindOrthogonalHitFromPathToPath for details


% Find a search radius to use - NEED TO AUTOMATE THIS TO PREVENT LOOPS
last_stations = zeros(Ntraversals,1);
for ith_traversal = 1:Ntraversals
    last_stations(ith_traversal,1) = data.traversal{ith_traversal}.Station(end);
end
std_end_station = std(last_stations);
search_radius = 5*std_end_station; % Expect adjacent paths to be within this many meters of the query point

% Initialize the iteration errors so we can record iteration convergence
iteration_error_X{num_iterations}  = 0;
iteration_error_Y{num_iterations}  = 0;

% Initialize the number of station points in this traversal
reference_station_points    = (0:interval:reference_traversal.Station(end))';

%% Perform iterations to seek an average path
% * Project from the reference path to nearby trajectories to find projections
% * Average projections to find new average path. 
% * Clean up average and store results for next iteration or exit

for ith_iteration =1:num_iterations 
    % Show user what we are doing?        
    if flag_do_debug
        fprintf(1,'Averaging paths via iteration: %.0d / %.0d \n',ith_iteration,num_iterations);
    end
    
    path_last  = reference_traversal;   % Saves the prior path - used for error calculations later     
    
    %% For each traversal, project from reference orthogonally
    [closestXs, closestYs, closestDistances] = ...
        fcn_Path_FindOrthogonalScatterFromPathToPaths(...
        reference_station_points, reference_traversal, data,...
        flag_rounding_type,search_radius);     
    
    %% Find and remove outliers due to distance jumps
    % Find the standard deviation in the distances
    sigma_distances = nanstd(closestDistances,0,'all');    
    outliers = find(abs(closestDistances)>(5*sigma_distances));
    if ~isempty(outliers)
        closestXs(outliers) = NaN;
        closestYs(outliers) = NaN;
        closestDistances(outliers) = NaN;         %#ok<NASGU>
    end
        
    %% For debugging: plot results thus far?
    if flag_do_debug        
        path_points_debug_fig = 22;
        figure(path_points_debug_fig); clf;
        fcn_Path_plotPathXY(data,path_points_debug_fig);
        hold on

        xlabel('x [m]');
        ylabel('y [m]');
        
        plot(reference_traversal.X,reference_traversal.Y,'r.','Markersize',25);

        for i_point = 1:length(closestXs(:,1))
            plot(closestXs(i_point,:),closestYs(i_point,:),'bo-');
        end
    end
    
    %% Average these projection points to generate an average path    
    raw_path_average = [nanmean(closestXs,2) nanmean(closestYs,2)];
    
    % Call a special function to remove back-tracking behavior
    path_average = fcn_Path_cleanForwardBackwardJogs(raw_path_average);        
    path_average_station  = [0; cumsum(sqrt(sum(diff(path_average).^2,2)))];
    
    % Show results?
    if flag_do_debug
        figure(path_points_debug_fig);        
        plot(path_average(:,1), path_average(:,2),'-','LineWidth',2)
        plot(path_average(:,1), path_average(:,2),'.','Markersize',15)
    end
    
    %% Redecimate the station points in this traversal
    reference_station_points    = [(0:interval:path_average_station(end))'; path_average_station(end)];   

    %% Interpolation of mean data to produce equal station intervals
    path_average_interp_X       = interp1(path_average_station,path_average(:,1),reference_station_points,'linear');
    path_average_interp_Y       = interp1(path_average_station,path_average(:,2),reference_station_points,'linear');
    path_average_interp_station = interp1(path_average_station,path_average_station,reference_station_points,'linear');

    
    %% Update path
    reference_traversal.X       = path_average_interp_X;
    reference_traversal.Y       = path_average_interp_Y;
    reference_traversal.Station = path_average_interp_station;

    %% Update error calculations for iterations 2 and onward
    if ith_iteration>=2
        shortest_length = min(length(path_last.X),length(reference_traversal.X));
        iteration_error_X{ith_iteration} = path_last.X(1:shortest_length) - reference_traversal.X(1:shortest_length);
        iteration_error_Y{ith_iteration} = path_last.Y(1:shortest_length) - reference_traversal.Y(1:shortest_length);
        
        % Show results?
        if flag_do_debug
            X_error = iteration_error_X{ith_iteration};
            Y_error = iteration_error_Y{ith_iteration};
            total_error = (X_error.^2 + Y_error.^2).^0.5;
            fprintf(1,'\t Mean change from iteration %.0d to %.0d: %.3f \n',ith_iteration-1, ith_iteration, mean(total_error));
        end
    end
    
    %% Plot the results?
    if flag_do_debug
        figure(path_points_debug_fig);
        plot(reference_traversal.X, reference_traversal.Y,'g.','Markersize',15);
        Nstations = length(reference_station_points);
        title(sprintf('Original paths and nearest points. Iteration %.0d of %.0d, with %.0d stations',ith_iteration,num_iterations,Nstations));
        pause(0.1)
    end
    
end

% Use final average path to define "true" s-coordinates of the original trajectories, using projection
path_average_final  = reference_traversal;

% Calculate final results
[closestXs, closestYs, closestDistances] = ...
    fcn_Path_FindOrthogonalScatterFromPathToPaths(...
    path_average_final.Station, reference_traversal, data,...
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
    fcn_Path_plotPathXY(data,path_points_debug_fig);
    hold on;
    plot(path_average_final.X,path_average_final.Y,'b.-','Linewidth',4,'Markersize',20);
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
        % Plot the error convergence
        figure(2233);
        clf;
        hold on;
        xlabel('Index')
        ylabel('Position change between iterations [m]')
        for ith_iteration = 1:length(iteration_error_X)
            title(sprintf('Iteration %.0d of %.0d',ith_iteration,length(iteration_error_X)));
            X_error = iteration_error_X{ith_iteration};
            Y_error = iteration_error_Y{ith_iteration};
            total_error = (X_error.^2 + Y_error.^2).^0.5;
            plot(total_error);
            
            pause(0.5);
        end
        
        
        % Plot the error convergence
        figure(2244);
        clf;
        hold on;
        xlabel('Iteration')
        ylabel('Mean change in distance')
        for ith_iteration = 2:length(iteration_error_X)
            X_error = iteration_error_X{ith_iteration};
            Y_error = iteration_error_Y{ith_iteration};
            total_error = (X_error.^2 + Y_error.^2).^0.5;
            plot(ith_iteration,mean(total_error),'.','Markersize',30);
        end
    end
    
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end
end

% ===============================================
%  nested function for yaw angle calculation (NOT needed for real data)
% ===============================================
function lane_yaw = fnc_yaw(X,Y) %#ok<*DEFNU>
% Yaw,north is zero, clockwise is positive direction, range 0-360 degrees

% Make sure that x and y are column vectors.
X=X(:);
Y=Y(:);
% calculate the atan2d and convert it to Yaw
yaw = 90 - atan2d(diff(Y), diff(X)); %
lane_yaw = [yaw(1); yaw];
n_yaw  = lane_yaw < 0;
lane_yaw(n_yaw) = lane_yaw(n_yaw)+360;
end



function unwrapped_angle = fcn_DataClean_unwrapAngles(wrapped)
initial_angle = wrapped(1,1);
change_in_angle = [0; diff(wrapped)];
index_jumps = find(change_in_angle>180);
change_in_angle(index_jumps) = change_in_angle(index_jumps)-360;
index_jumps = find(change_in_angle<-180);
change_in_angle(index_jumps) = change_in_angle(index_jumps)+360;
unwrapped_angle = cumsum(change_in_angle) + initial_angle;

% Shift all data up or down 
mean_angle = mean(unwrapped_angle);
good_mean = mod(mean_angle,360);
shift = mean_angle - good_mean;
unwrapped_angle = unwrapped_angle - shift;
end

function clean_path_average = fcn_Path_cleanForwardBackwardJogs(path_average)
% Find and remove outliers due to angle jumps

iteration_count = 1;
flag_average_is_good = 0;

while (0==flag_average_is_good)  && (iteration_count<=3)
    % Calculate angle changes between points
    diff_angles = fcn_Path_calcDiffAnglesBetweenPathSegments(path_average);
    
    % Find outliers
    outliers = find(abs(diff_angles)>pi/4);
        
    if ~isempty(outliers)
        % ID adjacent points
        outliers = unique([outliers;outliers+1;outliers-1]);
        outliers = min(outliers,length(path_average)-1);
        outliers = max(1,outliers);
        
        % Create a set of indices we will save
        indices = (1:length(path_average(:,1)))';
        indices(outliers+1) = 0;
        
        % Save the clean path
        clean_path_average = path_average(indices~=0,:);        
        
        % For debugging
        if 1==0
            figure(8888);
            clf;
            hold on;
            grid on;
            grid minor;
            
            x_indices = 1:length(diff_angles);
            plot(x_indices,diff_angles,'k-');
            plot(x_indices(outliers), diff_angles(outliers),'ro');
        end
        
    else
        flag_average_is_good = 1;
        clean_path_average = path_average;
    end
    
    % Show results for debugging?
    if 1==0
        figure(777777);
        clf;
        hold on;
        grid on;
        grid minor;
        axis equal;
        
        plot(path_average(:,1),path_average(:,2),'k.-');
        plot(path_average(outliers+1,1),path_average(outliers+1,2),'ro');
        plot(clean_path_average(:,1),clean_path_average(:,2),'b-');
    end
    
    % Increment the iteration count
    iteration_count = iteration_count + 1;
    
    % Reset the path average for the next round
    path_average = clean_path_average;
end

end
