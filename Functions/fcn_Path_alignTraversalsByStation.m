function diff_angles = fcn_Path_alignTraversalsByStation(Path,varargin)
% fcn_Path_calcDiffAnglesBetweenSegments
% Calculates the change in angles between segments. If there are N points
% in the path, there are N-1 segments and thus N-2 angles between segments.
%
% Note that this method uses the dot product and the cross product to avoid
% errors caused by angle rollover that occur with arctan calculations.
%
% FORMAT: 
%
%       diff_angles = fcn_Path_calcDiffAnglesBetweenPathSegments(Path,(fig_num))
%
% INPUTS:
%
%      Path: an N x 2 vector with [X Y] data in each row. N must be >= 3.
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      diff_angles: tan (N-2) x 1 vector of the change in angles in radians
%
% EXAMPLES:
%      
%       See the script:
%       script_test_fcn_Path_calcDiffAnglesBetweenPathSegments.m for a full
%       test suite.
%
% This function was written on 2021_01_25 by Satya Prasad
% Questions or comments? szm888@psu.edu 

% Revision history:
%     2021_01_25
%     -- first writing of the code


flag_do_debug = 0; % Flag to plot the results for debugging
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
    if nargin < 1 || nargin > 2
        error('Incorrect number of input arguments')
    end
    
    % Check the Path variables        
    fcn_Path_checkInputsToFunctions(Path, 'paths');

end

% Does user want to show the plots?
if 2 == nargin
    fig_num = varargin{1};
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

% Note that this method uses the dot product and the cross product to avoid
% errors caused by angle rollover. One could use the arctan method of
% calculating angles, which works in general but fails for paths that point
% straight to the left due to the cross-over point for the atan2
% calculation. For example, the following code does NOT always work:
%
%     % angles = atan2(path_average(2:end,2)-path_average(1:end-1,2),path_average(2:end,1)-path_average(1:end-1,1));
%     % angles = [angles; angles(end)];  % Pad the last point twice
%     % diff_angles2 = abs(diff(angles));

            % initializing variables
            station_offsets = zeros(1,length(data.traversal)); 
            id_offset = zeros(1,length(data.traversal));
            station_offsets(1) = 0;
            id_offset(1) = 1;
            
            % Loop through all traversals other than the 1st one. For each
            % traversal, find the closest point on the first traversal, and
            % use this to calculate station offset and index offset of the
            % ith traversal relative to the first traversal
            
            for i = 2:length(data.traversal)
                % calculating x, y, and z position of one specific point on
                % each traversal given by index_to_match (a common point)
                X = data.traversal{i}.X(index_to_match);
                Y = data.traversal{i}.Y(index_to_match);
                Z = data.traversal{i}.Z(index_to_match);
                
                % finding the idx of the first traversal that is closet to
                % the current point
                [id, station, lateral_offset, heading_map] = obj.stationHere(X,Y,Z,0,0,0,data.traversal{1}.station,data.traversal{1}.X,data.traversal{1}.Y,data.traversal{1}.yaw,"xyz");
                station_offsets(i) = data.traversal{i}.station(index_to_match) - station; % distance offset of the current trajectory staiton from the first traj station 
                id_offset(i) = id(1); % index offset
                
                % correct the station by subtracting the offset
                data.traversal{i}.station = data.traversal{i}.station - station_offsets(i); 

            end
            
            if crop_to_common_station == 1
                
                % And we will chop off the beginning and ends of the data so they all have
                % the same length of data. First we need to find the maximum and minimum
                % station.
                start_stations = zeros(1,length(data.traversal));
                end_stations = zeros(1,length(data.traversal));
                for i = 1:length(data.traversal)
                    start_stations(i) = data.traversal{i}.station(1);
                    end_stations(i) = data.traversal{i}.station(end);
                end
                max_start_station = max(start_stations);
                min_end_station = min(end_stations);

                % Now crop the data by going through all the fields and
                % keeping only the indices that are within the bounds.
                for i = 1:length(data.traversal)

                   inds = find(data.traversal{i}.station >= max_start_station & data.traversal{i}.station <= min_end_station);

                    % Trim the data to the desired times
%                    data.traversal{i}.time = data.traversal{i}.time(inds) - data.traversal{i}.time(inds(1));
                    data.traversal{i}.station = data.traversal{i}.station(inds);
%                   data.traversal{i}.latitude = data.traversal{i}.latitude(inds);
%                   data.traversal{i}.longitude = data.traversal{i}.longitude(inds);
%                   data.traversal{i}.altitude = data.traversal{i}.altitude(inds);
%                   data.traversal{i}.east_velocity = data.traversal{i}.east_velocity(inds);
%                   data.traversal{i}.north_velocity = data.traversal{i}.north_velocity(inds);
%                   data.traversal{i}.U = data.traversal{i}.U(inds);
%                   data.traversal{i}.V = data.traversal{i}.V(inds);
                    data.traversal{i}.X = data.traversal{i}.X(inds);
                    data.traversal{i}.Y = data.traversal{i}.Y(inds);
                    data.traversal{i}.Z = data.traversal{i}.Z(inds);
                    data.traversal{i}.yaw = data.traversal{i}.yaw(inds);
%                   data.traversal{i}.yaw_rate = data.traversal{i}.yaw_rate(inds);

                end
                
            end % Ends if crop to common station
            
            if interpolate_station == 1 % Do we interpolate?
                
                % Interpolate the data to be on the same station
                % decimation. First we need to find the minimum and maximum
                % station.
                start_stations = zeros(1,length(data.traversal));
                end_stations = zeros(1,length(data.traversal));
                for i = 1:length(data.traversal)
                    start_stations(i) = data.traversal{i}.station(1);
                    end_stations(i) = data.traversal{i}.station(end);
                end
                min_start_station = min(start_stations);
                max_end_station = max(end_stations);

                new_station = min_start_station:station_decimation:max_end_station;
                for i = 1:length(data.traversal)

                    % It's possible that there are non-unique station
                    % measurements, so you cannot interpolate the data.
                    % This grabs only the unique data indices
                    [~,unique_inds] = unique(data.traversal{i}.station);

%                   data.traversal{i}.latitude = interp1(data.traversal{i}.station(unique_inds), data.traversal{i}.latitude(unique_inds), new_station,'linear','extrap');
%                   data.traversal{i}.longitude = interp1(data.traversal{i}.station(unique_inds), data.traversal{i}.longitude(unique_inds), new_station,'linear','extrap');
%                   data.traversal{i}.altitude = interp1(data.traversal{i}.station(unique_inds), data.traversal{i}.altitude(unique_inds), new_station,'linear','extrap');
%                   data.traversal{i}.east_velocity = interp1(data.traversal{i}.station(unique_inds), data.traversal{i}.east_velocity(unique_inds), new_station,'linear','extrap');
%                   data.traversal{i}.north_velocity = interp1(data.traversal{i}.station(unique_inds), data.traversal{i}.north_velocity(unique_inds), new_station,'linear','extrap');
%                   data.traversal{i}.U = interp1(data.traversal{i}.station(unique_inds), data.traversal{i}.U(unique_inds), new_station,'linear','extrap');
%                   data.traversal{i}.V = interp1(data.traversal{i}.station(unique_inds), data.traversal{i}.V(unique_inds), new_station,'linear','extrap');
                    data.traversal{i}.X = interp1(data.traversal{i}.station(unique_inds), data.traversal{i}.X(unique_inds), new_station,'linear','extrap');
                    data.traversal{i}.Y = interp1(data.traversal{i}.station(unique_inds), data.traversal{i}.Y(unique_inds), new_station,'linear','extrap');
                    data.traversal{i}.Z = interp1(data.traversal{i}.station(unique_inds), data.traversal{i}.Z(unique_inds), new_station,'linear','extrap');
                    data.traversal{i}.yaw = interp1(data.traversal{i}.station(unique_inds), data.traversal{i}.yaw(unique_inds), new_station,'linear','extrap');
%                   data.traversal{i}.yaw_rate = interp1(data.traversal{i}.station(unique_inds), data.traversal{i}.yaw_rate(unique_inds), new_station,'linear','extrap');
                    data.traversal{i}.station = new_station - min_start_station;

                end
                
            end


%% Any debugging?
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
    xlabel('Index');
    ylabel('Angle [deg]');

    % Plot the angle differences
    plot(diff_angles*180/pi,'k.-','Linewidth',3,'Markersize',25);

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end
end

%% Calculate cross products
function result = crossProduct(v,w)
result = v(:,1).*w(:,2)-v(:,2).*w(:,1);
end

