% Authors: Dr. Brennan, Guangwei, Liming
% Revision history
% 2020_05_20 - fixed bug on the yaw angle plots
% 2020_06_20 - add raw data query functions
%% Clear the command window and workspace

clear all %#ok<CLALL>

% Add driver for database
javaaddpath('C:\Users\Guangwei Zhou\AppData\Roaming\MathWorks\MATLAB\R2020a\drivers\postgresql-42.2.9.jar')

%%
flag.DBquery = true; %set to true if you want to query raw data from database insteading of loading from default *.mat file
try
    fprintf('Starting the code using variable rawData of length: %d\n', length(rawData));
catch
    if flag.DBquery == true
        
        %---- Step1 :CONNECT TO  DATABASE ------------------------ %
        % choose different database name to connect to them
        database_name = 'mapping_van_raw';
        
        % Connect raw data database
        MDB = MapDatabase(database_name); % instance of MapDatabase Class
        % MDB.db.db_connection % show connection details
        
        % show all tables schema in the database
        tables = MDB.db.ShowTables();
        
        %---- Step2 :show all trips and pick the ones you want to query --------- %
        
        % query parameters
        MDB.zero_time = 0;  % offset the timestamps starting from zero
        MDB.verbose = 1;  % show the code processing details
        MDB.convert_GPS_to_ENU = 1; % you can choose the reference point by assign value to options.ENU_ref 
        MDB.separate_by_lap = 0; % default is 0
        
        % check all the trips
        trips = MDB.fetchTrips();
        
        % ----> pick trips you want to query(choose from trips table manually)
        % trip_names = {'Test Track Decision Points with Lane Change MappingVan 2020-03-13','Test Track MappingVan 2019-10-19'};
        trip_names = {'Test Track Decision Points with Lane Change MappingVan 2020-03-13'};
       
        prompt = ['Are you going to query the data of trip : \n' trip_names '\n?[y/n]'];
        User_input = input(strjoin(prompt),'s'); 
        
        if strcmpi(User_input,'y')
            fprintf(1,'Thanks. Let''s query it...\n');
        else
            fprintf(1,'Query is aborted. \nYou can Re-pick the trips name.\n');
            return
        end
        
        %---- Step3 :query data by trips ------------------------ %
        % quey trips id according to the trip name
        trip_id = [];
        for i = 1:length(trip_names)
            trip_id = cat(1,trip_id,trips.id(strcmp(trips.name, trip_names(i))));
        end
        
        % Pick sensors you want to query. 1 means query data from that sensor
        options = {};
        options.sensors.base_station = 1;  % default is 1
        options.sensors.hemisphere_gps = 1; % default is 1
        options.sensors.NovAtel_gps = 1; % default is 1
        options.sensors.garmin_gps = 1; % default is 1
        options.sensors.garmin_velocity = 1; % default is 0
        options.sensors.steering_angle =1; % default is 1
        options.sensors.NovAtel_imu = 1;% default is 1
        options.sensors.adis_imu = 1;% default is 1
        options.sensors.encoder_left_right = 1;% default is 1
        options.sensors.laser = 0; % default is 0
        options.sensors.front_left_camera = 0; % default is 0
        options.sensors.front_right_camera = 0; % default is 0
        options.sensors.front_center_camera = 0; % default is 0
        options.ENU_ref = 0; % 0 use default setting in database, 1 test track, 2 LTI, Larson  Transportation Institute
        
        % fetchByTripID
        result = MDB.fetchByTripID(trip_id,options);
        % disconnect with DB
        MDB.disconnect();
        %---- Step4 :pre process the data ------------------------ %
        % Notes: need Guangwei's help
        rawData = fcn_preProcessQueryResult(result);
        
    else
        % add the data path
        addpath '../data'
        
        % Load the raw data
        % This data will have outliers, be unevenly sampled, have multiple and inconsistent measurements of the same variable.
        filename  = 'MappingVan_DecisionMaking_03132020.mat';
        variable_names = 'MappingVan_DecisionMaking_03132020';
        rawData = fcn_loadRawData(filename,variable_names);

    end
    
    rawDataTimeFixed = fcn_removeTimeGapsFromRawData(rawData);
end

% Data clean and merge
% Fill in the sigma values for key fields. This just calculates the sigma values for key fields (velocities,
% accelerations, angular rates in particular), useful for doing outlier detection, etc. in steps that follow.
rawDataWithSigmas = fcn_loadSigmaValuesFromRawData(rawDataTimeFixed);


% NOTE: the following function changes the yaw angles to wind (correctly)
% up or down)

% Remove outliers on key fields via median filtering
% This removes outliers by median filtering key values.
rawDataWithSigmasAndMedianFiltered = fcn_medianFilterFromRawAndSigmaData(rawDataWithSigmas);

% PLOTS to show winding up or down:
% figure(2); plot(mod(rawDataWithSigmas.GPS_Novatel.Yaw_deg,360),'b')
% figure(3); plot(mod(rawDataWithSigmasAndMedianFiltered.GPS_Novatel.Yaw_deg,360),'k')
% figure(4); plot(rawDataWithSigmasAndMedianFiltered.GPS_Novatel.Yaw_deg,'r')

% Clean the raw data
cleanData = fcn_cleanRawDataBeforeTimeAlignment(rawDataWithSigmasAndMedianFiltered);

% Align all time vectors, and make time a "sensor" field
cleanAndTimeAlignedData = fcn_alignToGPSTimeAllData(cleanData);

% Time filter the signals
timeFilteredData = fcn_timeFilterData(cleanAndTimeAlignedData);

% Calculate merged data via Baysian averaging across same state
mergedData = fcn_mergeTimeAlignedData(timeFilteredData);

% Remove jumps from merged data caused by DGPS outages
mergedDataNoJumps = fcn_removeDGPSJumpsFromMergedData(mergedData,rawData);

% Calculate the KF fusion of single signals
mergedByKFData = mergedDataNoJumps;  % Initialize the structure with prior data

% KF the yawrate and yaw together
t_x1 = mergedByKFData.MergedGPS.GPS_Time;
x1 = mergedByKFData.MergedGPS.Yaw_deg;
x1_Sigma = mergedByKFData.MergedGPS.Yaw_deg_Sigma;
t_x1dot = mergedByKFData.MergedIMU.GPS_Time;
x1dot = mergedByKFData.MergedIMU.ZGyro*180/pi;
x1dot_Sigma = mergedByKFData.MergedIMU.ZGyro_Sigma*180/pi;
nameString = 'Yaw_deg';
[x_kf,sigma_x] = fcn_KFmergeStateAndStateDerivative(t_x1,x1,x1_Sigma,t_x1dot,x1dot,x1dot_Sigma,nameString);
mergedByKFData.MergedGPS.Yaw_deg = x_kf;
mergedByKFData.MergedGPS.Yaw_deg_Sigma = sigma_x;


% KF the xEast_increments and xEast together
t_x1 = mergedByKFData.MergedGPS.GPS_Time;
x1 = mergedByKFData.MergedGPS.xEast;
x1_Sigma = mergedByKFData.MergedGPS.xEast_Sigma;
t_x1dot = mergedByKFData.MergedGPS.GPS_Time;
x1dot = mergedByKFData.MergedGPS.xEast_increments/0.05;  % The increments are raw changes, not velocities. Have to divide by time step.
x1dot_Sigma = mergedByKFData.MergedGPS.xEast_increments_Sigma/0.05;
nameString = 'xEast';
[x_kf,sigma_x] = fcn_KFmergeStateAndStateDerivative(t_x1,x1,x1_Sigma,t_x1dot,x1dot,x1dot_Sigma,nameString);
mergedByKFData.MergedGPS.xEast = x_kf;
mergedByKFData.MergedGPS.xEast_Sigma = sigma_x;

% KF the yNorth_increments and yNorth together
t_x1 = mergedByKFData.MergedGPS.GPS_Time;
x1 = mergedByKFData.MergedGPS.yNorth;
x1_Sigma = mergedByKFData.MergedGPS.yNorth_Sigma;
t_x1dot = mergedByKFData.MergedGPS.GPS_Time;
x1dot = mergedByKFData.MergedGPS.yNorth_increments/0.05;  % The increments are raw changes, not velocities. Have to divide by time step.
x1dot_Sigma = mergedByKFData.MergedGPS.yNorth_increments_Sigma/0.05;
nameString = 'yNorth';
[x_kf,sigma_x] = fcn_KFmergeStateAndStateDerivative(t_x1,x1,x1_Sigma,t_x1dot,x1dot,x1dot_Sigma,nameString);
mergedByKFData.MergedGPS.yNorth = x_kf;
mergedByKFData.MergedGPS.yNorth_Sigma = sigma_x;


%% Merge data together by laps
lane_change.Clocks = mergedByKFData.Clocks;
lane_change.MergedGPS = mergedByKFData.MergedGPS;
loops = lane_change;
fn = fieldnames(lane_change.MergedGPS);
for j = 1:length(fn)
    fn_str = string(fn(j));
    lane_change.MergedGPS.(fn_str) = lane_change.MergedGPS.(fn_str)(1:end);
end

% define parameters of s-coordinate (start point and etc.)
% Deinfe the start and end point of one loop
RouteStructure = fcn_defineRouteStructure("lane_change_full_loop");


%% Break loop data into laps
[lap_data, num_laps] = fcn_breakFilteredDataIntoLaps(lane_change, RouteStructure);% lapData only contains lane-change segment?

% FOR DEBUGGING on 2020_05_29 - delete later if no longer needed
% figure(333332);
% clf;
% hold on;
% for i_Laps = 1:numLaps
%     plot(lapData{i_Laps}.MergedGPS.Yaw_deg,'-','LineWidth', 1);
% 
%     %legend_string{i_Laps} = sprintf('Lap %d',i_Laps);
% end

%% Apply the Algorithm in script_test_path_averaging
figure(2)
clf
hold on
for i_path= 1:num_laps
    plot(lap_data{i_path}.MergedGPS.xEast,lap_data{i_path}.MergedGPS.yNorth)
end

title('path')
xlabel('x [m]')
ylabel('y [m]')

% Pick the trajectory with the most data points so that the average path converges faster
data_length = zeros(1, num_laps);
for i_path = 1:num_laps
    data_length(i_path) = length(lap_data{i_path}.MergedGPS.xEast);
end
[~,initital_reference_path_id] = max(data_length);
path = lap_data{initital_reference_path_id}; %initial reference path

% run the simulation
fprintf(1,'\nAverage path is generating ...\n');
% profile off
average_start = tic;
% profile on -timer 'performance'
nbs_iteration = 2; % number of iterations to find the average path, at least 2
for i =1:nbs_iteration
    % Average these nearby points to generate an average path
    [closestXs, closestYs, closestZs, closestYaws] = fcn_findClosestPointsFromPath(path, lap_data, 1,0);
    
    path_average = [mean(closestXs,2) mean(closestYs,2) mean(closestZs,2)];
    path_average_yaw = mean(closestYaws,2);
    path_average_station  = [0; cumsum(sqrt(sum(diff(path_average).^2,2)))];
    
    % store last path 
    path_last  = path;
    
    % Interpolation of mean data, X,Y,Z
    % path_average_interp is same as path_average except the number of points
    interval = 0.5; % meters
    nb_points = round(path_average_station(end)/interval);
    
    path_average_interp = fcn_interparc(nb_points,path_average(:,1), path_average(:,2),path_average(:,3),'spline');
    
    % interplate the yaw
    path_average_interp_station = [0; cumsum(sqrt(sum(diff(path_average_interp).^2,2)))];
    path_average_interp_yaw = interp1(path_average_station, path_average_yaw, path_average_interp_station,'spline','extrap');
       
    % Project from the average path to nearby trajectories to find projections
    path.MergedGPS.xEast = path_average_interp(:,1);
    path.MergedGPS.yNorth = path_average_interp(:,2);
    path.MergedGPS.zUp = path_average_interp(:,3);
    path.station = [0; cumsum(sqrt(sum(diff(path_average_interp).^2,2)))];
    path.MergedGPS.Yaw_deg = path_average_interp_yaw;
    
    hold on
    plot(path_average_interp(:,1), path_average_interp(:,2),'-o','LineWidth',1)
    
    %pause(0.1)
end
% profile viewer
time_elapsed = toc(average_start);
fprintf(1,'Path averaging is done.\n\tWall time to run is: %.5f seconds. \n',time_elapsed);

%% plot the original path, averaged path, and the nearest points
figure(3)
clf
hold on
for i_path= 1:num_laps
    %plot(lap_data{i_path}.MergedGPS.xEast,lap_data{i_path}.MergedGPS.yNorth,'-')
    plot(closestXs(:,i_path), closestYs(:,i_path),'-s')
end
plot(path_last.MergedGPS.xEast,path_last.MergedGPS.yNorth,'r-o','LineWidth',3)
hold off
title('original path, averaged path, and the nearest points')
grid on
xlabel('x [m]')
ylabel('y [m]')
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use the data in path_last  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_mean = path_last.MergedGPS.xEast;
y_mean = path_last.MergedGPS.yNorth;
yaw_mean =  path_last.MergedGPS.Yaw_deg;

sigma = ((std(closestXs,0,2)).^2 + (std(closestYs,0,2)).^2).^0.5; % need confirm

%{
%% Average across all laps to find means, Notes: need edit fcn_meanStationProjection
% aligned_Data_ByStation and mean_Data only contains lane-change segment
[aligned_Data_ByStation, mean_Data] = fcn_meanStationProjection(lap_data,num_laps, RouteStructure); 

% obtain mean path data
x_mean = mean_Data.mean_xEast;
y_mean = mean_Data.mean_yNorth;
% calculate sigma (standard deviation)
sigma = fcn_calculateSigma(aligned_Data_ByStation,mean_Data,num_laps);


% Calculate lateral offset paths
s_coef = 1;
[x_offset_left, y_offset_left, x_offset_right, y_offset_right] = fcn_calculatePathOffsets(x_mean, y_mean, sigma, s_coef);
offset = 0.5*((x_offset_right - x_offset_left).^2 + (y_offset_right - y_offset_left).^2).^0.5;
%}

%% Test if closestSXY can be Used to Perform Histogram Algorithm
figure(8567)
hold on
test_idx = 200;
for i_laps = 1:num_laps
    plot(closestXs(:,i_laps), closestYs(:,i_laps))
    plot(closestXs(test_idx,i_laps), closestYs(test_idx,i_laps), 'ko')
end

plot(x_mean,y_mean,'k','LineWidth',2)
plot(x_mean(test_idx),y_mean(test_idx),'ro','LineWidth',2);

%% Get the error for yaw angle and position
% NOTES:
% Need to fix median filter - uses signal processing toolbox
% Need to fix the time error in length (see the -2 fix in the time
% conversion)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOT YET UPDATED BELOW
mean_projection_vector = [diff(x_mean) diff(y_mean)]; % [x(1)-x(2), y(1)-y(2); x(2)-x(3), y(2)-y(3);...]
mean_projection_vector = [mean_projection_vector(1,:); mean_projection_vector];  
mean_projection_vector = [mean_projection_vector zeros(length(mean_projection_vector(:,1)), 1)]; % Turn the vector into an xyz vector by adding a col of zeros

%all_dist = [];
all_error = [];
error_lap = [];
yaw_lap = [];
yaw_error_lap_unwrapped = [];
all_yaw_error = [];

for i_lap = 1:num_laps
    
    % Grab the deviation vector:[x-x_mean, y-y_mean, 0]
    deviation_vector = [closestXs(:,i_lap)-x_mean, closestYs(:,i_lap)-y_mean, 0*y_mean];
   
    % Grab the yaw angle for this run
    yaw_this_traversal = closestYaws(:, i_lap);

    % Calculate the yaw angle errors for this run
    yaw_error_this_traversal = closestYaws(:, i_lap) - yaw_mean;
    yaw_error_this_traversal_unwrapped = mod(yaw_error_this_traversal+180,360)-180;
    
    % Calculate the distance from the mean
    dist_from_mean = sum(deviation_vector.^2,2).^0.5;
    %all_dist = [all_dist, dist_from_mean];
    
    % Do the cross product so we can check which direction the error is in
    cross_product_projection_to_deviation = cross(mean_projection_vector,deviation_vector);
    
    % Keep just the z-component - the 3rd column
    sign_of_deviation = 2.0*(cross_product_projection_to_deviation(:,3)>0)-1;
    
    % Grab the error
    error_relative_to_path = dist_from_mean.*sign_of_deviation;
    
    % get the closest point from the mean on each lap?
    
    error_lap = [error_lap, error_relative_to_path];
    all_error = [all_error; error_relative_to_path]; %#ok<AGROW>
    all_yaw_error = [all_yaw_error; yaw_error_this_traversal]; %#ok<AGROW>
    yaw_lap = [yaw_lap, yaw_this_traversal]; %#ok<AGROW>
    yaw_error_lap_unwrapped = [yaw_error_lap_unwrapped, yaw_error_this_traversal_unwrapped]; %#ok<AGROW>
    
end

flag_eliminate_lane_changing = 1;
degree_threshold = 3; % Units are degrees
if flag_eliminate_lane_changing
    indices_not_changing_lanes = find(abs(mod(all_yaw_error+180,360)-180)<degree_threshold);
    all_error = all_error(indices_not_changing_lanes);
end

%% Classify all data points along the mean path to be either 1 lane area, 2 lane area, or transition area
yaw_error_diff = diff(yaw_error_lap_unwrapped);
yaw_error_std = std(yaw_error_lap_unwrapped')';

% Run a classifier on lane to force short segments to be ignored. The
% repelem function counts how many replecations of a designation occur in a
% sequence. We require 100 data points or more in order for a segment to be
% classified as 1 lane or 2 lanes.

is_one_lane = yaw_error_std < 1.2;
d1 = [true; diff(is_one_lane) ~= 0; true]; % TRUE if values change
num_rep1 = diff(find(d1));               
frequency_of_one = repelem(num_rep1, num_rep1); % if is_one_lane=[1,1,1,0,0,0,0,1], then frequency_of_one = [3,3,3,4,4,4,4,1]


is_two_lane = (yaw_error_std < 5) & (yaw_error_std > 2.5);
d2 = [true; diff(is_two_lane) ~= 0; true];  
num_rep2 = diff(find(d2));               
frequency_of_two = repelem(num_rep2, num_rep2);


% Eliminate '1's that only appear in short subarrays
for i = 1:length(yaw_error_lap_unwrapped)
    if is_one_lane(i) == 1 && frequency_of_one(i) < 100
        is_one_lane(i) = 0;
    end
    
    if is_two_lane(i) == 1 && frequency_of_two(i) < 100
        is_two_lane(i) = 0;
    end
end

lane_classification = zeros(length(x_mean), 1); % Identify the road segment to be 1 lane (1), 2 lanes (2), or transition area (0)
lane_classification(is_one_lane == 1) = 1;
lane_classification(is_two_lane == 1) = 2;

% Do the classifications here
one_lane_area = find(lane_classification == 1);
two_lane_area = find(lane_classification == 2);
transition_area = find(lane_classification == 0);

% Find the start and end location(index) of each area
one_lane_start = [];
one_lane_end = [];
two_lane_start = [];
two_lane_end = [];
transition_start = [];
transition_end = [];

% Store the classification into vectors, storing the "start" locations for
% one-lane, the "end" locations for one lane, etc.

for i = 1:length(lane_classification)
    % Store the start locations
    if i == 1 || lane_classification(i) ~= lane_classification(i-1)
        if lane_classification(i) == 1
            one_lane_start = [one_lane_start; i];
        elseif lane_classification(i) == 2
            two_lane_start = [two_lane_start; i];
        else
            transition_start = [transition_start; i];
        end
    end
    
    
    % Store the end locations
    if i == length(lane_classification) || lane_classification(i) ~= lane_classification(i+1)
        if lane_classification(i) == 1
            one_lane_end = [one_lane_end; i];
        elseif lane_classification(i) == 2
            two_lane_end = [two_lane_end; i];
        else
            transition_end = [transition_end; i];
        end
    end
end


%% Use Histogram Bincount to Determine Lane Positions (find the WIDTH of the lanes)
% For 2-lane area, there will be two peaks for the 2 histograms in this
% area. The peaks should correspond to the center of each lane.
% Select the data in all_error that corresponds to 2 Lane area
% The highest positive bin represents location of the right lane
% The highest negative bin represents location of the left lane

all_error_two_lane = error_lap(two_lane_start:two_lane_end, :);
all_error_two_lane = all_error_two_lane(:); % convert to 1D array

% Grab the histogram counts so that we can analyze the data
[counts, edge] = histcounts(all_error_two_lane, 75); % length of edge is always 1 more than count, fix by finding center of bin

% find the center location of each bin
bin_center =  0.5 * (edge(1:end-1) + edge(2:end));
right_lane_offset = bin_center(counts == max(counts(bin_center>0)));
left_lane_offset = bin_center(counts == max(counts(bin_center<0)));

lane_width = (abs(right_lane_offset)+abs(left_lane_offset))/2;
path_lane_width = zeros(length(x_mean),1);
path_lane_width(one_lane_area) = lane_width;
path_lane_width(transition_area) = 1.5*lane_width;
path_lane_width(two_lane_area) = 2*lane_width;

%% First Subplot: plot of lap data and average path
figure(1234);
subplot(2,3,1)
hold on
test_idx = 200;

% This is a plot of all the data separated into laps
legend_string = ''; % Initialize an empty string for legends
for i_laps = 1:num_laps
    plot(closestXs(:,i_laps), closestYs(:,i_laps))
    plot(closestXs(test_idx,i_laps), closestYs(test_idx,i_laps), 'ko')
    legend_string{i_laps} = sprintf('Lap %d',i_laps);
end

plot(x_mean,y_mean,'k','LineWidth',2)
plot(x_mean(test_idx),y_mean(test_idx),'ro','LineWidth',2);


xlabel('GPS x Coordinate [m]') %set  x label 
ylabel('GPS y Coordinate [m]') % set y label 
%legend(legend_string);
grid minor;

%% Second Subplot: Plot of all yaw error with respect to data index
subplot(2,3,2)
%xlim([0, length(all_yaw_error_unwrapped)]);
for i = 1:num_laps
    plot(yaw_error_lap_unwrapped(:, i));% All yaw angle plot
    hold on;
end

xlabel("Data Index");
ylabel("Yaw Angle Error from Mean [Degree]");


%% Third Subplot: plot standard deviation of the yaw error and lane classification boundaries
subplot(2,3,3)
plot(yaw_error_std)
xlim([0, length(yaw_error_std)]);
xlabel("Data Index");
ylabel("Standard Deviation of Yaw Angle Error [Degree]");
grid minor;
hold on;

%plot of Lane Classification
p(1) = plot(one_lane_area, zeros(length(one_lane_area),1), 'r.', 'MarkerSize',10);
p(2) = plot(two_lane_area, zeros(length(two_lane_area),1), 'g.', 'MarkerSize',10);
p(3) = plot(transition_area, zeros(length(transition_area),1), 'b.', 'MarkerSize',10);
legend(p(1:3), 'One Lane Area', 'Two Lane Area', 'Transition Area')

%% Forth Subplot: Lateral Offset histogram
subplot(2,3,4)
histogram(all_error_two_lane,75); 
title("75 Bins")
xlabel("Error Relative to Path [m]")
ylabel("Count")

%% Fifth Subplot: Polygons that indicate the lane classifications
subplot(2,3,5)
hold on
fcn_plotPathWithVarianceShading(x_mean(one_lane_area),y_mean(one_lane_area),'y', path_lane_width(one_lane_area),1);
fcn_plotPathWithVarianceShading(x_mean(two_lane_area),y_mean(two_lane_area),'r', path_lane_width(two_lane_area),1);
fcn_plotPathWithVarianceShading(x_mean(transition_area),y_mean(transition_area),'g', path_lane_width(transition_area),1);
for i_laps = 1:num_laps
    plot(closestXs(:,i_laps), closestYs(:,i_laps))
end
