function [lapData,numLaps] = fcn_breakDataIntoLaps_v1(clean_data,start_point)
%fcn_break_data_into_laps Summary of this function goes here
%    break given data into laps and calculate station 
%   Detailed explanation goes here
%% Find distances to start
start_xEast = start_point(3);
start_yNorth = start_point(4);
distances_to_start_in_meters = ((clean_data.xEast - start_xEast).^2 + (clean_data.yNorth - start_yNorth).^2).^0.5;
indices_data = (1:length(distances_to_start_in_meters))';

% Positions close to the start line will be within the distance threshold
distance_threshold = 5; % In meters

% Grab only the data that is close to the starting point
closeToStartLineIndices = find(distances_to_start_in_meters<distance_threshold);  

% Fill in the distances to these indices
close_locations = distances_to_start_in_meters*NaN;
close_locations(closeToStartLineIndices) = distances_to_start_in_meters(closeToStartLineIndices);

% Find the inflection point by finding crossover point in the derivative
changing_direction = diff(close_locations);
crossing_points = find(changing_direction(1:end-1).*changing_direction(2:end)<0);

% Plot the crossing points result (for debugging)
h_fig = figure(99947464);
set(h_fig,'Name','Distance_to_startPoint');
plot(indices_data,distances_to_start_in_meters,'b'); 
hold on;
plot(indices_data,close_locations,'r'); 
plot(indices_data(crossing_points),distances_to_start_in_meters(crossing_points),'co'); 
grid on
ylim([-1 25]);
hold on;
ylabel('Distance to start point [m]');
xlabel('Index of trip [unitless]');

%% Unpack the data into laps

numLaps = length(crossing_points) -1;
for i_Laps = 1:numLaps
    indices_for_lap = crossing_points(i_Laps):crossing_points(i_Laps+1);
    lapData{i_Laps}.GPS_Time = clean_data.GPS_Time(indices_for_lap); 
    lapData{i_Laps}.longitude = clean_data.Longitude(indices_for_lap); 
    lapData{i_Laps}.latitude = clean_data.Latitude(indices_for_lap); 
    lapData{i_Laps}.altitude = clean_data.Altitude(indices_for_lap);
    
    lapData{i_Laps}.xEast = clean_data.xEast(indices_for_lap); 
    lapData{i_Laps}.yNorth =clean_data. yNorth(indices_for_lap); 
    lapData{i_Laps}.zUp = clean_data.zUp(indices_for_lap); 
    
    lapData{i_Laps}.velNorth = clean_data.velNorth(indices_for_lap); 
    lapData{i_Laps}.velEast = clean_data.velEast(indices_for_lap); 
    lapData{i_Laps}.velUp = clean_data.velUp(indices_for_lap); 

    lapData{i_Laps}.yaw_angle_from_velocity =clean_data.Yaw_deg_from_velocity(indices_for_lap);
    
end

% calculate station 
for i_Laps = 1:numLaps
    lapData{i_Laps}.station(1) = 0; %sqrt((lapData{i_Laps}.clean_xEast(1))^2+ (lapData{i_Laps}.clean_yNorth(1) )^2); 
    for i = 2:length(lapData{i_Laps}.GPS_Time)
        delta_station= sqrt((lapData{i_Laps}.xEast(i) - lapData{i_Laps}.xEast(i-1))^2+ (lapData{i_Laps}.yNorth(i) - lapData{i_Laps}.yNorth(i-1))^2); 
        lapData{i_Laps}.station(i) =  lapData{i_Laps}.station(i-1) + delta_station;
%         
%                     station = sqrt(diff(X).^2 + diff(Y).^2 + diff(Z).^2);
%             station = cumsum(station);
%             station = [0; station];
    end
end

end

