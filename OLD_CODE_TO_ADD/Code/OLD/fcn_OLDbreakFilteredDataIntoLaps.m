function [lapData,numLaps] = fcn_breakFilteredDataIntoLaps(FilteredData,RouteStructure)
%fcn_break_data_into_laps Summary of this function goes here
%    break given data into laps and calculate station 

%% grad data
start_xEast =  RouteStructure.start_xEast;
start_yNorth = RouteStructure.start_yNorth;

xEast = FilteredData.MergedGPS.xEast;
yNorth = FilteredData.MergedGPS.yNorth;

%% Find distances to start

distances_to_start_in_meters = ((xEast - start_xEast).^2 + (yNorth - start_yNorth).^2).^0.5;

%indices_data = (1:length(distances_to_start_in_meters))';

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

%{
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
%}

%% Unpack the data into laps
names = fieldnames(FilteredData); 
numLaps = length(crossing_points) -1;

time_ROS = FilteredData.Clocks.targetTimeVector_ROS{5};
time_GPS = FilteredData.Clocks.targetTimeVector_GPS{5};

for i_Laps = 1:numLaps
    
    Indices_for_lap = crossing_points(i_Laps):crossing_points(i_Laps+1);
    laps_GPStime_start(i_Laps) = time_GPS(Indices_for_lap(1));
    laps_GPStime_end(i_Laps) = time_GPS(Indices_for_lap(end));
    laps_ROStime_start(i_Laps) = time_ROS (Indices_for_lap(1));
    laps_ROStime_end(i_Laps) = time_ROS(Indices_for_lap(end));
           
    for i_data = 1: length(names) %i_sensor 
        
        data_name = names{i_data};
        d = FilteredData.(data_name); % extract field  data 
        subfieldNames = fieldnames(d); % Grab all the subfields
        
        if 1==i_data  %break Clock structrue
            lapData{i_Laps}.(data_name).startTimeGPS = laps_GPStime_start(i_Laps);
            lapData{i_Laps}.(data_name).endTimeGPS =  laps_GPStime_end(i_Laps);
            lapData{i_Laps}.(data_name).startTimeROS = laps_ROStime_start(i_Laps);
            lapData{i_Laps}.(data_name).endTimeROS =laps_ROStime_end(i_Laps);
            num = 5;
                for i_time = 1: length(FilteredData.Clocks.targetTimeVector_ROS)
                    t1 = FilteredData.Clocks.targetTimeVector_ROS{i_time} ; %copy 
                    t2 = FilteredData.Clocks.targetTimeVector_GPS{i_time} ; %copy 
                    index= find(t2 >= laps_GPStime_start(i_Laps) & t2 <=  laps_GPStime_end(i_Laps));

                    lapData{i_Laps}.(data_name).(subfieldNames{num}){i_time} = t1(index);
                    lapData{i_Laps}.(data_name).(subfieldNames{num+1}){i_time} =t2(index);
                end
            continue
        end 
         %lapData{i_Laps}.(data_name) = lapData{i_Laps}.(names{1});
       
        time= FilteredData.Clocks.targetTimeVector_GPS{5};
        index_time_laps= time >= laps_GPStime_start(i_Laps) & time <=  laps_GPStime_end(i_Laps);
        
        for i_subField =1:length(subfieldNames) %setdiff([7:length(subfieldNames)],[8,9,29])
            % Grab the name of the ith subfield
            subFieldName = subfieldNames{i_subField};
            dout.(subFieldName) = d.(subFieldName);  % Copy over the field itself first
            data_mid = d.(subFieldName); 
            % Trim the data according to laps time 
            if  length(data_mid) == 1
                lapData{i_Laps}.(data_name).(subFieldName) =data_mid;
            else
                lapData{i_Laps}.(data_name).(subFieldName) =data_mid(index_time_laps);
            end
           % timeFilteredData.(data_name) = dout; % Save results to main structure
        end
        index_time_laps=[];
    end
end

%calculation station 
for i_Laps = 1:numLaps
    station = sqrt((diff(lapData{i_Laps}.MergedGPS.xEast)).^2 + (diff(lapData{i_Laps}.MergedGPS.yNorth)).^2);
    station = cumsum(station);
    station = [0; station];
    lapData{i_Laps}.station=  station;
    station = [];
end


% calculate station 
% for i_Laps = 1:numLaps
%     lapData{i_Laps}.station(1) = 0; %sqrt((lapData{i_Laps}.clean_xEast(1))^2+ (lapData{i_Laps}.clean_yNorth(1) )^2); 
%     lapData{i_Laps}.timeFilteredData.GPS_Hemisphere.velMagnitude
%     
%     for i = 2:length(lapData{i_Laps}.GPS_Time)
%         delta_station= sqrt((lapData{i_Laps}.xEast(i) - lapData{i_Laps}.xEast(i-1))^2+ (lapData{i_Laps}.yNorth(i) - lapData{i_Laps}.yNorth(i-1))^2); 
%         lapData{i_Laps}.station(i) =  lapData{i_Laps}.station(i-1) + delta_station;
% %         
% %                     station = sqrt(diff(X).^2 + diff(Y).^2 + diff(Z).^2);
% %             station = cumsum(station);
% %             station = [0; station];
%     end
% end

% Eliminate incorrectly identified laps which have very little data points in them
n = numel(lapData);
idx = 1;
for i = 1:n
    if numel(lapData{1,i}.station) > 100
        temp{idx} = lapData{1,i};
        %plot(lapData{1,i}.MergedGPS.xEast, lapData{1,i}.MergedGPS.yNorth)
        idx = idx + 1;
    end
end

lapData = temp;
numLaps = numel(lapData);

end

