function [aligned_Data_ByStation,mean_Data] = fcn_Path_findAverageTraversalViaStationAlignment(data)
% fcn_Path_findAverageTraversalViaStationAlignment
% finds the average traversal of several traversals by averaging the
% lateral offsets at the same stations.
%
% FORMAT: 
%
%      [aligned_Data_ByStation,mean_Data] = fcn_Path_findAverageTraversalViaStationAlignment(data);
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
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%      (fill in)
%
% DEPENDENCIES:
%
%      (fill in)
%
% EXAMPLES:
%      
%     See the script: 
%     script_test_fcn_Path_findAverageTraversalViaStationAlignment
%     script_test_fcn_Path_findAveragePath
%     for a full test suite.
%
% This function was written on 2020_05_29 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%     2020_05_29:
%     -- fixed mod issue in yaw angle
%     2020_11_14
%     -- changed name to fcn_Path_findAveragePathViaStationAlignment
%     2021_01_02
%     -- changed input arguments for consistency
%     2022_01_03
%     -- updated header to more modern format

% for i = 1:numLaps % data structure convert
%     data.traversal{i} = lapData{i};      
% end

numLaps = length(data.traversal);

% Parameters for aligning each traversal by station
index_to_match = round(length(data.traversal{1}.X)/2); % A random point on each path to match to the first traversal
crop_to_common_station = 1;
interpolate_station = 1;
station_decimation = 0.5;

% Plot parameters
%font_size = 16;
%line_width = 2;
%sigma_multiplier = 2; % Specifies the number of sigmas to plot


% Create instance of Map class and process the GPS data to calculate
% station and project the latitude/longitude/altitude points into global XYZ
% We will set the first GPS reading as our 'base station' to calculate the
% global XYZ coordinates and station.

map_config = {};
%         map_config.reference_latitude = data.traversal{1}.latitude(1);
%        map_config.reference_longitude = data.traversal{1}.longitude(1);
%         map_config.reference_altitude = data.traversal{1}.altitude(1);

m = Map(map_config);
%start_offset_times = [0 0 0 0 0 0 0];
%end_offset_times = [0 0 0 0 0 0 0];

%% Define offset times to trim from the start and end respectively
% Then trim the data to only keep indices within those boundaries, and fill
% in the "data" structure

% start_offset_times = 0;  % Start time offset in seconds
% end_offset_times = 0;
for ith_Lap = 1:numLaps
       % data.traversal{i} = m.processGPS(data.traversal{i});% calculate, station, velocity,yaw and ENU
        % time_trip = data.traversal{i}.Clocks.targetTimeVector_GPS{5}; 
        %inds = find(time_trip >= start_offset_times(i) & time_trip <= time_trip(end) - end_offset_times(i));
        % inds = find(time_trip >= start_offset_times & time_trip <= time_trip(end) - end_offset_times);
        inds = (1:length(data.traversal{ith_Lap}.X))';
        % Trim the data to the desired times
        % data.traversal{i}.time = time_trip(inds);  
        data.traversal{ith_Lap}.station = data.traversal{ith_Lap}.Station(inds);
        %data.traversal{i}.X = data.traversal{i}.MergedGPS.xEast(inds); 
        %data.traversal{i}.Y = data.traversal{i}.MergedGPS.yNorth(inds);
        data.traversal{ith_Lap}.Z = data.traversal{ith_Lap}.Y*0;
        data.traversal{ith_Lap}.yaw = data.traversal{ith_Lap}.Y *0;
               
        % Zero each time
        %data.traversal{i}.time = data.traversal{i}.time - data.traversal{i}.time(1);
    
end

%% Align the data by station 
% this is within the Map.m subfunctions
aligned_Data_ByStation = m.alignDataByStation(data, index_to_match, crop_to_common_station, interpolate_station, station_decimation);     
fprintf('Data is aligned!\n');
    
%calculate the mean
all_xEast_measurements   = zeros(numLaps,length(aligned_Data_ByStation.traversal{1}.station));
all_yNorth_measurements  = all_xEast_measurements;
all_zUp_measurements     = all_xEast_measurements;
all_station_measurements = all_xEast_measurements;
all_yaw_measurements     = all_yNorth_measurements;

for ith_Lap = 1:numLaps
    all_xEast_measurements(ith_Lap,:) = aligned_Data_ByStation.traversal{ith_Lap}.X'; 
    all_yNorth_measurements(ith_Lap,:) = aligned_Data_ByStation.traversal{ith_Lap}.Y'; 
    all_zUp_measurements(ith_Lap,:) = aligned_Data_ByStation.traversal{ith_Lap}.Z'; 
    all_station_measurements(ith_Lap,:) = aligned_Data_ByStation.traversal{ith_Lap}.station';
    %all_yaw_measurements(i,:) = aligned_Data_ByStation.traversal{i}.yaw';
    all_yaw_measurements(ith_Lap,:) = aligned_Data_ByStation.traversal{ith_Lap}.yaw';  % Fixed this on 2020_05_29
    
    % Debugging plot to check for errors
    if 1==0
        figure(8945794);
        hold on;
        subplot(2,1,1);
        plot(all_xEast_measurements(ith_Lap,:),all_yNorth_measurements(ith_Lap,:));
        subplot(2,1,2);
        plot(all_yaw_measurements(ith_Lap,:));
    end

end

%average them 
mean_Data.mean_station = mean(all_station_measurements)';
mean_Data.mean_xEast = mean(all_xEast_measurements)';
mean_Data.mean_yNorth = mean(all_yNorth_measurements)';
mean_Data.mean_zUp = mean(all_zUp_measurements)';
mean_Data.mean_yaw = mean(all_yaw_measurements)';

% Place plotting code in main script
%{
legend_string = ''; % Initialize an empty string
for i_Laps = 1:numLaps
    plot(aligned_Data_ByStation.traversal{i_Laps}.X,aligned_Data_ByStation.traversal{i_Laps}.Y,'-','LineWidth', 1);
    hold on;
    legend_string{i_Laps} = sprintf('Lap %d',i_Laps);
end
grid on; 
xlabel('x [m]') %set  x label 
ylabel('y [y]') % set y label 

plot(mean_Data.mean_xEast ,mean_Data.mean_yNorth,'k','LineWidth', 1.2)
plot(RouteStructure.start_xEast ,RouteStructure.start_yNorth,'go','MarkerSize', 14)

title('Plot of Average Path'); 
legend_string{numLaps+1} =  'mean';
legend_string{numLaps+2} =  'startPoint';
legend(legend_string);

h_fig = figure(16700);
set(h_fig,'Name','mean_yaw_in_Laps');
legend_string = ''; % Initialize an empty string
for i_Laps = 1:numLaps
   plot(aligned_Data_ByStation.traversal{i_Laps}.station,aligned_Data_ByStation.traversal{i_Laps}.yaw,'-','LineWidth', 1.2);
   hold on;
   legend_string{i_Laps} = sprintf('Lap %d',i_Laps);
end
plot(aligned_Data_ByStation.traversal{1}.station ,mean_Data.mean_yaw,'r','LineWidth', 2.4)
grid on; 
xlabel('station [m]') %set  x label 
ylabel('yaw [deg]') % set y label 
title('Plot of mean\_yaw\_in\_Laps'); 
legend([legend_string, 'mean_route']);

%}
end


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
