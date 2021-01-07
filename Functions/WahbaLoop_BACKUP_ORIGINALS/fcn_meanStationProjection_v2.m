function [aligned_Data_ByStation,mean_Data] = fcn_meanStationProjection_v2(lapData,numLaps)
%UNTITLED4 Summary of this function goes here

%   Detailed explanation goes here
%%%%calculate mean station by laps for time filtered data (not completed)

for i = 1:numLaps % data structure convert
    data.traversal{i}=lapData{i};        
end

% Parameters for aligning each traversal by station
index_to_match = 500; % A random point on each path to match to the first traversal
crop_to_common_station = 1;
interpolate_station = 1;
station_decimation = 0.5;

% Plot parameters
font_size = 16;
line_width = 2;
sigma_multiplier = 2; % Specifies the number of sigmas to plot


            % Create instance of Map class and process the GPS data to calculate
            % station and project the latitude/longitude/altitude points into global XYZ
            % We will set the first GPS reading as our 'base station' to calculate the
            % global XYZ coordinates and station.
         map_config = {};
%         map_config.reference_latitude = data.traversal{1}.latitude(1);
 %        map_config.reference_longitude = data.traversal{1}.longitude(1);
%         map_config.reference_altitude = data.traversal{1}.altitude(1);

        m = Map(map_config);  %%==============??????
        start_offset_times = [0 0 0 0 0 0 0];
        end_offset_times = [0 0 0 0 0 0 0];
for i = 1:numLaps
       % data.traversal{i} = m.processGPS(data.traversal{i});% calculate, station, velocity,yaw and ENU
        
        inds = find(data.traversal{i}.GPS_Time >= start_offset_times(i) & data.traversal{i}.GPS_Time <= data.traversal{i}.GPS_Time(end) - end_offset_times(i));
        
        % Trim the data to the desired times
        data.traversal{i}.time = data.traversal{i}.GPS_Time(inds);
        data.traversal{i}.station = data.traversal{i}.station(inds);
        data.traversal{i}.latitude = data.traversal{i}.latitude(inds);
        data.traversal{i}.longitude = data.traversal{i}.longitude(inds);
        data.traversal{i}.altitude = data.traversal{i}.altitude(inds);
        data.traversal{i}.east_velocity = data.traversal{i}.velEast(inds);
        data.traversal{i}.north_velocity = data.traversal{i}.velNorth(inds);
        %data.traversal{i}.U = data.traversal{i}.U(inds);
        %data.traversal{i}.V = data.traversal{i}.V(inds);
        data.traversal{i}.X = data.traversal{i}.xEast(inds);
        data.traversal{i}.Y = data.traversal{i}.yNorth(inds);
        data.traversal{i}.Z = data.traversal{i}.zUp(inds);
        data.traversal{i}.yaw = data.traversal{i}.yaw_angle_from_velocity(inds);
       % data.traversal{i}.yaw_rate = data.traversal{i}.yaw_rate(inds);
        
        % Zero each time
        data.traversal{i}.time = data.traversal{i}.time - data.traversal{i}.time(1);
    
end

    aligned_Data_ByStation = m.alignDataByStation(data, index_to_match, crop_to_common_station, interpolate_station, station_decimation);
       
    fprintf('Data is aligned!\n');
    
    %calculate the mean
    all_xEast_measurements = zeros(numLaps,length(aligned_Data_ByStation.traversal{1}.station));
    all_yNorth_measurements= all_xEast_measurements;
    all_yaw_measurements = all_yNorth_measurements;
for i = 1:numLaps
    all_xEast_measurements(i,:) = aligned_Data_ByStation.traversal{i}.X'; 
    all_yNorth_measurements(i,:) = aligned_Data_ByStation.traversal{i}.Y'; 
    all_yaw_measurements(i,:) = aligned_Data_ByStation.traversal{i}.yaw';
end

mean_Data.mean_xEast = mean(all_xEast_measurements);
mean_Data.mean_yNorth = mean(all_yNorth_measurements);
mean_Data.mean_yaw = mean(all_yaw_measurements);
 %
 h_fig = figure(167);
set(h_fig,'Name','mean_ENU_in_Laps');
legend_string = ''; % Initialize an empty string
for i_Laps = 1:numLaps
    plot(aligned_Data_ByStation.traversal{i_Laps}.X,aligned_Data_ByStation.traversal{i_Laps}.Y,'-','LineWidth', 1.2);
    hold on;
    legend_string{i_Laps} = sprintf('Lap %d',i_Laps);
end
plot(mean_Data.mean_xEast ,mean_Data.mean_yNorth,'r','LineWidth', 2.4)
plot(mean_Data.mean_xEast(860) ,mean_Data.mean_yNorth(860),'go','MarkerSize', 14)
grid on; 
xlabel('time [s]') %set  x label 
ylabel('station [m]') % set y label 
title('Plot of mean ENU'); 
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
title('Plot of mean_yaw_in_Laps'); 
legend(legend_string, 'mean');


end

