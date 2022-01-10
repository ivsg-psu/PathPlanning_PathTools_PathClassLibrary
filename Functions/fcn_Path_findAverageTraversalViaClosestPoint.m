function average_traversal = fcn_Path_findAverageTraversalViaClosestPoint(data, varargin)
% fcn_Path_findAverageTraversalViaClosestPoint
% finds the average of several traversals by taking a reference traversal
% (or, if one is not given, using the traversal with longest number of
% points) and for each point in the traversal finding the nearest point in
% other traversals. The nearest point is determined by the "snap" of the
% reference traversal vertices to the closest location of each path. Thus,
% the resulting projection is orthogonal to each individual nearby path
% (and not usually orthogonal to the reference path).
%
% FORMAT: 
%
%      path_average_final = ...
%      fcn_Path_findAverageTraversalViaClosestPoint(...
%      data,...
%      reference_traversal,...
%      (flag_calculate_reference_traversal),(num_iterations))
%
% INPUTS:
%
%      data: a structure containing i traversal fields, each with subfields
%      of X, Y, etc. in the following form
%           data.traversal{i_path}.X
%      Note that i_path denotes an array of paths. Each path will be
%      compared separately
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
%      iteration's average, etc. Usually, convergence occurs within 3
%      iterations (the default is 10 for this function). Note that, within
%      the function, a debug flag can be set to plot and analyze
%      convergence.
%
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      path_average_final: the resulting traversal
%
% DEPENDENCIES:
%
%      fcn_Path_findClosestPointsToTraversal
%      fcn_Path_findTraversalWithMostData
%      fcn_Path_plotTraversalsXY
%
% EXAMPLES:
%      
%     See the scripts: 
%     script_test_fcn_Path_findAverageTraversalViaClosestPoint
%     script_test_fcn_Path_findAverageTraversal
%     for a full test suite.
%
% This function was written on 2020_11_13 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%      2020_11_13: 
%      -- wrote the code
%      2021_01_07
%      -- changed name and output to reflect traversal, not path
%      -- added figure number as optional input
%      -- added input checking
%     2022_01_03:
%     -- corrected typos in comments

% TO_DO:
%      -- Make the reference_traversal an optional input (see the
%      orthogonal projection code for an example)

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

% Set defaults
num_iterations = 10;  % the number of iteration to find the average path 
flag_calculate_reference_traversal = 1;

if flag_check_inputs == 1
    % Are there the right number of inputs?
    if nargin < 1 || nargin > 5
        error('Incorrect number of input arguments')
    end
    
    % Check the data input
    fcn_Path_checkInputsToFunctions(data, 'traversals');
    
end

% Check to see if the reference_traversal was specified?
if nargin>=2
    temp = varargin{1};    
    if ~isempty(temp)
        reference_traversal = temp;
        if flag_check_inputs == 1           
            % Check the reference_traversal input
            fcn_Path_checkInputsToFunctions(reference_traversal, 'traversal');
        end 
        flag_calculate_reference_traversal = 0;
    end
end

% Check to see if the number of iterations was specified?
if nargin >= 3
    temp = varargin{2};
    if ~isempty(temp)
        num_iterations = temp;
    end
end
        
% Does user want to show the plots?
if 4 == nargin
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

%% Solve for the circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check to see if we need to calculate a reference traversal or if it was
% given
if flag_calculate_reference_traversal
    index_of_longest = fcn_Path_findTraversalWithMostData(data);
    reference_traversal = data.traversal{index_of_longest}; %initial reference path
end

% Initialize the iteration errors
iteration_error_X{num_iterations}  = 0;
iteration_error_Y{num_iterations}  = 0;

% Perform a loop of the following:
% * Project from the reference path to nearby trajectories to find projections
% * Average projections to find new average path. Repeat until average
% path is not changing

for ith_iteration =1:num_iterations 
    % Show user what we are doing
    fprintf(1,'Averaging paths via iteration: %.0d / %.0d ',ith_iteration,num_iterations);
    
    % Save the current reference traversal (so we can check changes after
    % this is updated)
    path_last  = reference_traversal;

    % Search nearest points from reference path to each of the other
    % traversals, saving the results as "closest" values
    [closestXs,closestYs,closestZs,closestYaws] = ...
        fcn_Path_findClosestPointsToTraversal(reference_traversal, data,1);
    
    if flag_do_debug
        path_points_fig = 22;
        figure(path_points_fig); clf;
        fcn_Path_plotPathXY(data,path_points_fig)
        hold on
        title(sprintf('Original paths and nearest points. Iteration %.0d of %.0d',ith_iteration,num_iterations));
        xlabel('x [m]')
        ylabel('y [m]')

        for i_point = 1:length(closestXs(:,1))
            %for i_path= 1:length(data.traversal)
                % plot(closestXs(i_point,i_path),closestYs(i_point,i_path),'bo')
                plot(closestXs(i_point,:),closestYs(i_point,:),'bo-')
            %end
            pause(0.001);
        end
    end
    
    % Average these projection points to generate an average path    
    path_average = [mean(closestXs,2) mean(closestYs,2) mean(closestZs,2)];
    
    %     traversal = fcn_Path_convertXYtoTraversalStructure(path1(:,1),path1(:,2));
    %     data.traversal{1} = traversal;
    
    path_average_station  = [0; cumsum(sqrt(sum(diff(path_average).^2,2)))];
    path_average_yaw = mean(closestYaws,2);
    
    if flag_do_debug
        figure(path_points_fig);        
        plot(path_average(:,1), path_average(:,2),'LineWidth',3)
    end
    
    % TO DO - put the following into a path interpolation function
    % Interpolation of mean data by equal intervals 
    interval = 1; % meters
    reference_station_points    = (0:interval:path_average_station(end))';

    path_average_interp_X       = interp1(path_average_station,path_average(:,1),reference_station_points,'spline');
    path_average_interp_Y       = interp1(path_average_station,path_average(:,2),reference_station_points,'spline');
    path_average_interp_Z       = interp1(path_average_station,path_average(:,3),reference_station_points,'spline');
    path_average_interp_yaw     = interp1(path_average_station,path_average_yaw,reference_station_points,'spline');
    path_average_interp_station = interp1(path_average_station,path_average_station,reference_station_points,'spline');
            
    % update path 
    reference_traversal.X       = path_average_interp_X;
    reference_traversal.Y       = path_average_interp_Y;
    reference_traversal.Z       = path_average_interp_Z;
    reference_traversal.Yaw     = path_average_interp_yaw;
    reference_traversal.Station = path_average_interp_station;
    
    % END INSERTION INTO path interpolation function
    
    % Update error calculations for iterations 2 and onward
    if ith_iteration>=2
        shortest_length = min(length(path_last.X),length(reference_traversal.X));
        iteration_error_X{ith_iteration} = path_last.X(1:shortest_length) - reference_traversal.X(1:shortest_length);
        iteration_error_Y{ith_iteration} = path_last.Y(1:shortest_length) - reference_traversal.Y(1:shortest_length);
        
        X_error = iteration_error_X{ith_iteration};
        Y_error = iteration_error_Y{ith_iteration};
        total_error = (X_error.^2 + Y_error.^2).^0.5;
        fprintf(1,'\t Mean change: %.3f \n',mean(total_error));
    else
        fprintf(1,'\t Mean change: Inf \n');
    end

    % Plot the results?
    if flag_do_debug
        figure(path_points_fig);
        plot(path_average_interp(:,1), path_average_interp(:,2),'.','LineWidth',1)
    end
    
    pause(0.1)
end

% Use final average path to define "true" s-coordinates of the original trajectories, using projection
average_traversal  = reference_traversal;


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
    path_points_fig = 22222;
    figure(path_points_fig); clf;
    fcn_Path_plotTraversalsXY(data,path_points_fig)
    hold on
    plot(average_traversal.X,average_traversal.Y,'Linewidth',4);
    title('original path and final average path')
    xlabel('X [m]')
    ylabel('Y [m]')
    
    % Plot the yaw results
    figure(111);
    clf
    hold on
    plot(average_traversal.Station, average_traversal.Yaw*180/pi,'r-' ,'LineWidth',2)    
    for i_path= 1:length(data.traversal)
        plot(data.traversal{i_path}.Station(1:end-1),data.traversal{i_path}.Yaw*180/pi)
    end
    legend('average traversal', 'original traversals')
    title('station vs yaw')
    xlabel('Station [m]')
    ylabel('Yaw [degree]')
    
    
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

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end
end % Ends the function
