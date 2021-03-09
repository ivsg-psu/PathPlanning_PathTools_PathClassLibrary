function [lapData,numLaps] = fcn_RouteSegments_breakTraversalIntoLaps(FilteredData,RouteStructure)
%fcn_RouteSegments_breakTraversalIntoLaps 
%     This function breaks traversal data into laps and calculates station.
%     It assumes that the input data has been filtered and cleaned, and has
%     the "FilteredData" structure type. See the data processing examples
%     for typical route data for examples.

% Revision history:
%     2020_05_20 - S. Brennan updated the code to have:
%     -- more comments, 
%     -- changed name of crossing_points to crossing_indices.
%     -- made the function more verbose
%     -- fixed bug with using fields that were not equi-sampled, which caused
%        averages to give wrong values later
%     -- fixed error where it was not using the end locations to calculate
%     segments
%     2020_05_29:
%     -- fixed mod issue in yaw angle
%     2020_11_10 
%     -- changed function name from fcn_breakFilteredDataIntoLaps to fcn_RouteSegments_breakTraversalIntoLaps in prep for DataClean class

flag_do_debug = 1;

%% Let the user know what we are doing
if flag_do_debug
    % Grab function name
    st = dbstack;
    namestr = st.name;

    % Show what we are doing
    fprintf(1,'\nWithin function: %s\n',namestr);
    fprintf(1,'\t Calculating the number of possible laps.\n');   
    fprintf(1,'\t Starting point is: \n');
    fprintf(1,'\t \t start_xEast:  %.4f meters in local coords.\n',RouteStructure.start_xEast);
    fprintf(1,'\t \t start_yNorth: %.4f meters in local coords.\n',RouteStructure.start_yNorth);
    fprintf(1,'\t \t end_xEast:  %.4f meters in local coords.\n',RouteStructure.end_xEast);
    fprintf(1,'\t \t end_yNorth: %.4f meters in local coords.\n',RouteStructure.end_yNorth);
end


%% grad data
% start and end points of the LAP are entered manually in defineRouteStructure
start_xEast  = RouteStructure.start_xEast;
start_yNorth = RouteStructure.start_yNorth;
end_xEast    = RouteStructure.end_xEast;
end_yNorth   = RouteStructure.end_yNorth;

xEast = FilteredData.MergedGPS.xEast;
yNorth = FilteredData.MergedGPS.yNorth;


%% Find distances to start

% TO DO - functionalize this section as it is repeated twice here!
% 
% This is done by calculating the distance of each point to the "start
% point" in a lap. This distance will get get smaller as we approach a
% lap's end, then bigger as we go away. So we take the derivative (e.g. the
% difference via "diff", and then find the inflection point in the
% derivative by checking the exact location where the derivative changes
% from positive to negative. At that point, the index is as close as
% possible to the "start point"

% Find Euclidean distance from lap's starting point to GPS data point
distances_to_start_in_meters = ((xEast - start_xEast).^2 + (yNorth - start_yNorth).^2).^0.5;

% Positions close to the start line will be within the following distance threshold
distance_threshold = 5; % In meters

% Grab only the data that is close to the starting point as measured by
% distance threshold
closeToStartLineIndices = find(distances_to_start_in_meters<distance_threshold);  

% Fill in the distances to these indices
close_start_locations = distances_to_start_in_meters*NaN; % Create an array with a bunch of NaN
close_start_locations(closeToStartLineIndices) = ...
    distances_to_start_in_meters(closeToStartLineIndices); % Fill the idx that are close to starting point

% Find the inflection point by finding crossover point in the derivative of
% distance. If approaching a crossing point, then distance will decrease
% and then increase. In this particular case, sequence of distance slopoes
% will go from negative to positive. We have to check for two cases
% however, and both occur when there are not enough points on approach to a
% starting point to enable detection of the derivative of distance. These
% are noted in the errors given below.

changing_direction = diff(close_start_locations);
if length(changing_direction(:,1))>=2
    start_crossing_points = find(changing_direction(1:end-1).*changing_direction(2:end)<0);
    if isempty(start_crossing_points)
        error('Insufficient points found within distance threshold. Unsure how to continue.');
    end
else
    error('Two or less points found within the distance threshold. Unsure how to continue as exact minimum requires three points to exist within the threshold');
end


%% Find distances to end
% This follows the same process as above start points
distances_to_end_in_meters = ((xEast - end_xEast).^2 + (yNorth - end_yNorth).^2).^0.5;

% Positions close to the end line will be within the distance threshold
distance_threshold = 5; % In meters

% Grab only the data that is close to the end point
closeToEndLineIndices = find(distances_to_end_in_meters<distance_threshold);  

% Fill in the distances to these indices
close_end_locations = distances_to_end_in_meters*NaN;
close_end_locations(closeToEndLineIndices) = distances_to_end_in_meters(closeToEndLineIndices);

% Find the inflection point by finding crossover point in the derivative
changing_direction = diff(close_end_locations);
if length(changing_direction(:,1))>=2
    end_crossing_points = find(changing_direction(1:end-1).*changing_direction(2:end)<0);
    if isempty(end_crossing_points)
        error('Insufficient points found within distance threshold. Unsure how to continue.');
    end
else
    error('Two or less points found within the distance threshold. Unsure how to continue as exact minimum requires three points to exist within the threshold');
end
    
%% Match starting and ending points up to each other
num_good_laps = 0;
start_end_index_pairs = []; % Contains the start and end idx of LAPS? NOT lane change areas?
for i=1:length(start_crossing_points)
    next_end_index = find(end_crossing_points>start_crossing_points(i),1);
    if ~isempty(next_end_index)
        num_good_laps = num_good_laps+1;
        start_end_index_pairs = [start_end_index_pairs; start_crossing_points(i) end_crossing_points(next_end_index)]; %#ok<AGROW>
    end
end

%% Debugging plots
%{
% Plot the crossing points result (for debugging)
h_fig = figure(99947464);
set(h_fig,'Name','Analysis of Start/End Point Calculation');
subplot(3,1,1);
hold on;
plot(xEast,yNorth,'b');
plot(start_xEast,start_yNorth,'go');
plot(end_xEast,end_yNorth,'ro');

subplot(3,1,2);
indices_data = (1:length(distances_to_start_in_meters))';
hold on;
plot(indices_data,distances_to_start_in_meters,'b'); 
plot(indices_data,close_start_locations,'r'); 
plot(indices_data(start_crossing_points),distances_to_start_in_meters(start_crossing_points),'co'); 
grid on
ylabel('Distance to start point [m]');
xlabel('Index of trip [unitless]');

subplot(3,1,3);
indices_data = (1:length(distances_to_end_in_meters))';
hold on;
plot(indices_data,distances_to_end_in_meters,'b'); 
plot(indices_data,close_end_locations,'r'); 
plot(indices_data(end_crossing_points),distances_to_end_in_meters(end_crossing_points),'co'); 
grid on
ylabel('Distance to start point [m]');
xlabel('Index of trip [unitless]');
%}

%% Unpack the data into laps
names = fieldnames(FilteredData); % Grab names of all the fields
numLaps = length(start_end_index_pairs(:,1));

if flag_do_debug
    fprintf(1,'\t Number of laps found is: %d\n',numLaps);
end

% Fill in the time vectors (used later in the loop. Target times represent
% the true GPS and ROS time stamps to use, for data that is sampled at 1,
% 2, 3, 4, and 5 hundreths of a second. There are 5 vectors of ROS and GPS
% target times for each situation, but may be more. The length of the
% targetTimeVector_ROS (or GPS) matches the number of points
time_ROS = FilteredData.Clocks.targetTimeVector_ROS; 
time_GPS = FilteredData.Clocks.targetTimeVector_GPS; 

% For laps, we need to pick a common time base - choose a 20 Hz sample rate
sample_every_deciseconds = 5;
time20Hz_ROS = time_ROS{sample_every_deciseconds};
time20Hz_GPS = time_GPS{sample_every_deciseconds};
num_timepoints = length(time20Hz_GPS(:,1));

%% Check to make sure filtered data is ALL at 20 Hz
if flag_do_debug
    fprintf(1,'\t Checking that all data is of correct sampling length before doing lap calculations.\n');
end
   
% Loop through all the data fields
for i_data = 1: length(names) %i_sensor
    
    % Grab the data field's name we are looking for
    data_name = names{i_data};
    fprintf(1,'\t \t Starting verification of data field: %s\n',data_name);
    
    d = FilteredData.(data_name); % extract field  data
    subfieldNames = fieldnames(d); % Grab all the subfields
    
    % Is the current data structure for time? If so, do NOT use it
    if ~strcmp(data_name,'Clocks')  % only enter beyond here if subfields OTHER than time
        for i_subField =1:length(subfieldNames) 
            % Grab the name of the ith subfield
            subFieldName = subfieldNames{i_subField};
            fprintf(1,'\t \t \t Verifying data subfield: %s\n',subFieldName);
            
            full_length_subdata = d.(subFieldName);

            % Trim the data according to GPS time, if data is a vector
            if  length(full_length_subdata) ~= 1 % Only look at vector data
                % Check if data needs to be resampled to 20 Hz?
                if length(full_length_subdata)~=num_timepoints
                    fprintf(1,'\t \t \t \t Subfield is of wrong length, resampling to 20 Hz\n');
                    fprintf(1,'\t \t \t \t so that it has a length target of: %d\n',num_timepoints);
                    deciseconds_used_by_sensor = 0;
                    for i_deciseconds = 1:length(FilteredData.Clocks.targetTimeVector_GPS)
                        if length(full_length_subdata) == length(FilteredData.Clocks.targetTimeVector_GPS{i_deciseconds})
                            deciseconds_used_by_sensor = i_deciseconds;
                        end
                    end % Ends for loop looking for decisecond value
                    if deciseconds_used_by_sensor == 0
                        error('Unable to match sensor to a time vector');
                    else
                        fprintf(1,'\t \t \t \t Subfield is originally sampled at %d deciseconds.\n', deciseconds_used_by_sensor);
                        full_length_subdata = interp1(...
                            FilteredData.Clocks.targetTimeVector_GPS{deciseconds_used_by_sensor},...
                            full_length_subdata,...
                            FilteredData.Clocks.targetTimeVector_GPS{sample_every_deciseconds},...
                            'nearest');
                        fprintf(1,'\t \t \t \t Subfield is now at %d deciseconds.\n', sample_every_deciseconds);
                        fprintf(1,'\t \t \t \t With length of %d indices.\n', length(full_length_subdata));
                    end % Ends check to see if deciseconds is zero
                end % Ends check to see if length of subdata is correct
                
                % Save data into FilteredData
                FilteredData.(data_name).(subFieldName) = full_length_subdata;
                
            end % Ends check to see if data is a vector
        end % Ends loop through subfields
    end % Ends if statement if the data is not clock data
end % Ends loop through the data structure




%% Fill in the lap data from filtered data
% Initialize the lapData structure
for i_Laps = 1:numLaps
    lapData{i_Laps} = FilteredData; %#ok<AGROW>
end

% At this point, lapData has 20 struct and EACH struct contains FULL LENGTH data

% Loop through the number of laps we've detected thus far
for i_Laps = 1:numLaps
    if flag_do_debug
        fprintf(1,'\t For lap: %d\n',i_Laps);
    end
    
    %Indices_for_lap = start_crossing_points(i_Laps):start_crossing_points(i_Laps+1);
    Indices_for_lap = start_end_index_pairs(i_Laps,1):start_end_index_pairs(i_Laps,2);
    if flag_do_debug
        fprintf(1,'\t \t Indices for lap go from: %d to %d\n',Indices_for_lap(1),Indices_for_lap(end));
    end
   

    % Loop through all the sensors, grabbing data just for the current lap       
    for i_data = 1: length(names) %i_sensor 
        
        % Grab the data field's name we are looking for
        data_name = names{i_data};
        fprintf(1,'\t \t Starting lap extraction for data field: %s\n',data_name);

        d = FilteredData.(data_name); % extract field  data 
        subfieldNames = fieldnames(d); % Grab all the subfields
        
        % Is the current data structure for time?
        if strcmp(data_name,'Clocks')  % update Clock structure for this specific lap
                      
            % Fill in the start/end times for this lap
            if flag_do_debug
                fprintf(1,'\t \t \t GPS time goes from: %.6f to %.6f\n',time20Hz_GPS(Indices_for_lap(1)),time20Hz_GPS(Indices_for_lap(end)));
                fprintf(1,'\t \t \t ROS time goes from: %.6f to %.6f\n',time20Hz_ROS(Indices_for_lap(1)),time20Hz_ROS(Indices_for_lap(end)));
                fprintf(1,'\t \t \t There are %d data points in this lap.\n',length(Indices_for_lap));                
            end
            
            lapData{i_Laps}.Clocks.startTimeGPS = time20Hz_GPS(Indices_for_lap(1));
            lapData{i_Laps}.Clocks.endTimeGPS   = time20Hz_GPS(Indices_for_lap(end));
            lapData{i_Laps}.Clocks.startTimeROS = time20Hz_ROS(Indices_for_lap(1));
            lapData{i_Laps}.Clocks.endTimeROS   = time20Hz_ROS(Indices_for_lap(end));
            
            % Need to copy the target time into each lap, for each time
            % vector that may be used for a target_time
            for i_time = 1:length(FilteredData.Clocks.targetTimeVector_GPS)
                lap_time_indices= find(time_GPS{i_time}>= time20Hz_GPS(Indices_for_lap(1)) & time_GPS{i_time} <=  time20Hz_GPS(Indices_for_lap(end)));                
                lapData{i_Laps}.Clocks.targetTimeVector_ROS{i_time} = lapData{i_Laps}.Clocks.targetTimeVector_ROS{i_time}(lap_time_indices);
                lapData{i_Laps}.Clocks.targetTimeVector_GPS{i_time} = lapData{i_Laps}.Clocks.targetTimeVector_GPS{i_time}(lap_time_indices);
            end       
            
            % Check that the lengths make sense
            if length(Indices_for_lap)~=length(lapData{i_Laps}.Clocks.targetTimeVector_GPS{sample_every_deciseconds})
                warning('Mismatch detected in GPS time length versus lap lengths.\n');
                warning('Expecting the 20Hz data to have length of: %d\n',length(lapData{i_Laps}.Clocks.targetTimeVector_GPS{sample_every_deciseconds}));
                warning('But data had a length of: %d\n',length(Indices_for_lap));
                error('Unable to continue...');                
            end
        else            
            % For all the non-time subfields OTHER than time, need to fill these in
            for i_subField =1:length(subfieldNames) %setdiff([7:length(subfieldNames)],[8,9,29])
                % Grab the name of the ith subfield
                subFieldName = subfieldNames{i_subField};
                fprintf(1,'\t \t \t Starting lap extraction for data subfield: %s\n',subFieldName);
                
                full_length_subdata = d.(subFieldName);
                % Trim the data according to laps time
                if  length(full_length_subdata) == 1
                    lapData{i_Laps}.(data_name).(subFieldName) =full_length_subdata;
                else
                   % Save data into lapData
                    lapData{i_Laps}.(data_name).(subFieldName) =full_length_subdata(Indices_for_lap);
                end
            end
        end
    end
end


%% calculate the station distance within each lap
for i_Laps = 1:numLaps
    station = sqrt((diff(lapData{i_Laps}.MergedGPS.xEast)).^2 + (diff(lapData{i_Laps}.MergedGPS.yNorth)).^2);
    station = cumsum(station);
    station = [0; station]; %#ok<AGROW>
    lapData{i_Laps}.station=  station;
    clear station; 
end

% calculate station (OLD METHOD)
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


%% calculate the yaw angle coordinates so that they align
for i_Laps = 1:numLaps
    % Find the starting yaw angle
    starting_angle = lapData{i_Laps}.MergedGPS.Yaw_deg(1);
    number_times_360 = floor(starting_angle/360);
    offset_yaw = number_times_360*360;
    lapData{i_Laps}.MergedGPS.Yaw_deg = lapData{i_Laps}.MergedGPS.Yaw_deg - offset_yaw;
end

%% Eliminate incorrectly identified laps which have very little data points in them
num_good_laps = 0; % Initialize number of laps
% Loop through each lap data
for i_Laps = 1:numel(lapData)
    % Check that number of data points is more than 100
    if numel(lapData{1,i_Laps}.station) > 100
        % If it is, increment number of good laps, and save into temp data
        % structure
        num_good_laps = num_good_laps + 1;
        temp{num_good_laps} = lapData{1,i_Laps}; %#ok<AGROW>
    end
end
% Copy temp data structure into final output, as well as number of good
% laps
lapData = temp;
numLaps = num_good_laps;

end

