%% Clear the command window
clc


%%
try
    fprintf('Starting the code using variable raawData of length: %d\n', length(rawData));
catch
    % add the data path
    addpath '../data';
    % Load the raw data
    filename  = 'MappingVan_DecisionMaking_02242020.mat';
    variable_names = 'MappingVan_DecisionMakingTestTrack_02242020';
    rawData = fcn_loadRawData(filename,variable_names);
    rawDataTimeFixed = fcn_removeTimeGapsFromRawData(rawData);
end

% Data clean and merge
% Fill in the sigma values for key fields. This just calculates the sigma values for key fields (velocities,
% accelerations, angular rates in particular), useful for doing outlier detection, etc. in steps that follow.
rawDataWithSigmas = fcn_loadSigmaValuesFromRawData(rawDataTimeFixed);

% Remove outliers on key fields via median filtering
% This removes outliers by median filtering key values.
rawDataWithSigmasAndMedianFiltered = fcn_medianFilterFromRawAndSigmaData(rawDataWithSigmas);

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


%% Separate Data into 5 Different Loops
% The start and end index are gotten from script_determineSeparationPoints
start_idx = [1,     15600, 42600, 56500, 68000];
end_idx   = [15370, 38550, 56500, 67500, 85950];

loop_1.Clocks = mergedByKFData.Clocks;
loop_1.MergedGPS = mergedByKFData.MergedGPS;
loop_2.Clocks = mergedByKFData.Clocks;
loop_2.MergedGPS = mergedByKFData.MergedGPS;
loop_3.Clocks = mergedByKFData.Clocks;
loop_3.MergedGPS = mergedByKFData.MergedGPS;
loop_4.Clocks = mergedByKFData.Clocks;
loop_4.MergedGPS = mergedByKFData.MergedGPS;
loop_5.Clocks = mergedByKFData.Clocks;
loop_5.MergedGPS = mergedByKFData.MergedGPS;


loops = [loop_1, loop_2, loop_3, loop_4, loop_5];
fields = ["loop_1","loop_2", "loop_3", "loop_4", "loop_5"];

main_branchPts = [];
for i = 1:length(loops)
    subplot(2,3,i)
    fn = fieldnames(loops(i).MergedGPS);
    for j = 1:length(fn)
        fn_str = string(fn(j));
        loops(i).MergedGPS.(fn_str) = loops(i).MergedGPS.(fn_str)(start_idx(i):end_idx(i));
    end
   
    % define parameters of s-coordinate (start point and etc.)
    loop_name = ['loop_', int2str(i)];
    RouteStructure = fcn_defineRouteStructure(loop_name);
    % break loop data into laps
    [lapData,numLaps] = fcn_OLDbreakFilteredDataIntoLaps(loops(i), RouteStructure);
    numLaps = 6;
    [aligned_Data_ByStation,mean_Data] = fcn_meanStationProjection(lapData,numLaps, RouteStructure);
    
    % obtain mean path data
    field_i = fields(i);
    x_mean.(field_i) = mean_Data.mean_xEast;
    y_mean.(field_i) = mean_Data.mean_yNorth;
    % calculate sigma
    sigma.(field_i) = fcn_OLDcalculateSigma(aligned_Data_ByStation,mean_Data,numLaps);
    
    %{
    % create field names for lapData (lap1, lap2, ...)
    lap_names = [];
    for k = 1:numLaps
        lap_names = [lap_names, string(['lap_',int2str(k)])];
    end
    %}
    loopsDataInLaps{i} = aligned_Data_ByStation.traversal;
    
    title(['Loop ', num2str(i), ' Data']);
end

%% Intersections 1 & 2
figure;
subplot(2,2,1);
[x_decision_12, y_decision_12] = fcn_calculateDecisionPoint(x_mean.loop_1, y_mean.loop_1, sigma.loop_1, x_mean.loop_2, y_mean.loop_2, sigma.loop_2, 2);
xlim([1120, 1280]);
ylim([6100, 6270])
grid on;
title('Y-Intersection 1 & 2');
box on

% Intersections 3
subplot(2,2,2)
[x_decision_3, y_decision_3] = fcn_calculateDecisionPoint(x_mean.loop_2, y_mean.loop_2, sigma.loop_2, x_mean.loop_3, y_mean.loop_3, sigma.loop_3)
grid on;
xlim([1035, 1085]);
ylim([6200, 6250]);
title('Y-Intersection 3');
box on

% Intersections 4 (cross road: straight & right)
subplot(2,2,3)
[x_decision_4sr, y_decision_4sr] = fcn_calculateDecisionPoint(x_mean.loop_2, y_mean.loop_2, sigma.loop_2, x_mean.loop_4, y_mean.loop_4, sigma.loop_4)
grid on;
xlim([1085, 1110]);
ylim([6320, 6345]);
title('Crossroad Straight & Right');
box on

% Intersections 4 (cross road: right & left)
subplot(2,2,4)
[x_decision_4sl, y_decision_4sl] = fcn_calculateDecisionPoint(x_mean.loop_2, y_mean.loop_2, sigma.loop_2, x_mean.loop_5(500:900), y_mean.loop_5(500:900), sigma.loop_5(500:900))
xlim([1085, 1110]);
ylim([6320, 6345]);
grid on;
title('Crossroad Left & Right');
box on


%% Save Branching Points Data
figure
hold on


main_branch_pts = [];

for i = 1:6
    [x_branch, y_branch] = fcn_branchingPointOfPaths(loopsDataInLaps{1}{i}.MergedGPS.xEast,loopsDataInLaps{1}{i}.MergedGPS.yNorth, ...
    loopsDataInLaps{2}{i}.MergedGPS.xEast, loopsDataInLaps{2}{i}.MergedGPS.yNorth);
    d1 = ((x_branch - x_decision_12(1)).^2 + (y_branch - y_decision_12(1)).^2).^0.5;
    d2 = ((x_branch - x_decision_12(2)).^2 + (y_branch - y_decision_12(2)).^2).^0.5;
    for j = 1:length(d1)
        if d1(j)< 25 | d2(j) <25
            main_branch_pts = [main_branch_pts; x_branch(j), y_branch(j)];
        end
    end
end
plot(x_mean.loop_1, y_mean.loop_1)
plot(x_mean.loop_2, y_mean.loop_2)


for i = 1:6
    [x_branch, y_branch] = fcn_branchingPointOfPaths(loopsDataInLaps{2}{i}.MergedGPS.xEast,loopsDataInLaps{2}{i}.MergedGPS.yNorth, ...
    loopsDataInLaps{3}{i}.MergedGPS.xEast, loopsDataInLaps{3}{i}.MergedGPS.yNorth);
    d1 = ((x_branch - x_decision_3(1)).^2 + (y_branch - y_decision_3(1)).^2).^0.5;
    for j = 1:length(d1)
        if d1(j)< 25
            main_branch_pts = [main_branch_pts; x_branch(j), y_branch(j)];
        end
    end
end
plot(x_mean.loop_2, y_mean.loop_2)
plot(x_mean.loop_3, y_mean.loop_3)

for i = 1:6
    [x_branch, y_branch] = fcn_branchingPointOfPaths(loopsDataInLaps{2}{i}.MergedGPS.xEast,loopsDataInLaps{2}{i}.MergedGPS.yNorth, ...
    loopsDataInLaps{4}{i}.MergedGPS.xEast, loopsDataInLaps{4}{i}.MergedGPS.yNorth);
    d1 = ((x_branch - x_decision_4sr(1)).^2 + (y_branch - y_decision_4sr(1)).^2).^0.5;
    for j = 1:length(d1)
        if d1(j)< 25
            main_branch_pts = [main_branch_pts; x_branch(j), y_branch(j)];
        end
    end
end
plot(x_mean.loop_2, y_mean.loop_2)
plot(x_mean.loop_4, y_mean.loop_4)



for i = 1:6
    [x_branch, y_branch] = fcn_branchingPointOfPaths(loopsDataInLaps{2}{i}.MergedGPS.xEast,loopsDataInLaps{2}{i}.MergedGPS.yNorth, ...
    loopsDataInLaps{5}{i}.MergedGPS.xEast, loopsDataInLaps{5}{i}.MergedGPS.yNorth);
    d1 = ((x_branch - x_decision_4sl(1)).^2 + (y_branch - y_decision_4sl(1)).^2).^0.5;
    for j = 1:length(d1)
        if d1(j)< 25
            main_branch_pts = [main_branch_pts; x_branch(j), y_branch(j)];
        end
    end
end
plot(x_mean.loop_2, y_mean.loop_2)
plot(x_mean.loop_5(500:900), y_mean.loop_5(500:900))
plot(main_branch_pts(:,1), main_branch_pts(:,2), 'ko')

save('main_branch_pts.mat', 'main_branch_pts')
%{
%% Clear the workspace
clc, clear

%% Load the raw data
% This data will have outliers, be unevenly sampled, have multiple and
% inconsistent measurements of the same variable. In other words, it is the
% raw data.
filename  = 'MappingVan_DecisionMakingTestTrack_02242020.mat';
variable_names = 'MappingVan_DecisionMakingTestTrack_02242020';
rawData = fcn_loadRawData(filename,variable_names);
rawDataTimeFixed = fcn_removeTimeGapsFromRawData(rawData);

%% Data clean and merge
% Fill in the sigma values for key fields. This just calculates the sigma values for key fields (velocities,
% accelerations, angular rates in particular), useful for doing outlier detection, etc. in steps that follow.
rawDataWithSigmas = fcn_loadSigmaValuesFromRawData(rawDataTimeFixed);

%% Remove outliers on key fields via median filtering
% This removes outliers by median filtering key values.
rawDataWithSigmasAndMedianFiltered = fcn_medianFilterFromRawAndSigmaData(rawDataWithSigmas);

%% Clean the raw data
cleanData = fcn_cleanRawDataBeforeTimeAlignment(rawDataWithSigmasAndMedianFiltered);

%% Align all time vectors, and make time a "sensor" field
cleanAndTimeAlignedData = fcn_alignToGPSTimeAllData(cleanData);

%% Time filter the signals
timeFilteredData = fcn_timeFilterData(cleanAndTimeAlignedData);

%% Calculate merged data via Baysian averaging across same state
mergedData = fcn_mergeTimeAlignedData(timeFilteredData);

%% Remove jumps from merged data caused by DGPS outages
mergedDataNoJumps = fcn_removeDGPSJumpsFromMergedData(mergedData,rawData);

%% Calculate the KF fusion of single signals
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


%% Separate Data into 5 Different Loops
% The start and end index are gotten from script_determineSeparationPoints
start_idx = [1,     15600, 42600, 57500, 68000, 91500];
end_idx   = [15370, 38550, 56500, 67500, 85950, 100000];

loop_1.Clocks = mergedByKFData.Clocks;
loop_1.MergedGPS = mergedByKFData.MergedGPS;
loop_2.Clocks = mergedByKFData.Clocks;
loop_2.MergedGPS = mergedByKFData.MergedGPS;
loop_3.Clocks = mergedByKFData.Clocks;
loop_3.MergedGPS = mergedByKFData.MergedGPS;
loop_4.Clocks = mergedByKFData.Clocks;
loop_4.MergedGPS = mergedByKFData.MergedGPS;
loop_5.Clocks = mergedByKFData.Clocks;
loop_5.MergedGPS = mergedByKFData.MergedGPS;

loops = [loop_1, loop_2, loop_3, loop_4, loop_5];
fields = ["loop_1","loop_2", "loop_3", "loop_4", "loop_5"];
loopsDataInLaps = {};

for i = 1:length(loops)
    subplot(2,3,i)
    fn = fieldnames(loops(i).MergedGPS);
    for j = 1:length(fn)
        fn_str = string(fn(j));
        loops(i).MergedGPS.(fn_str) = loops(i).MergedGPS.(fn_str)(start_idx(i):end_idx(i));
    end
    
    length(loops(i).MergedGPS.xEast)
    % define parameters of s-coordinate (start point and etc.)
    loop_name = ['loop_', int2str(i)];
    RouteStructure = fcn_defineRouteStructure(loop_name);
    % break loop data into laps
    [lapData,numLaps] = fcn_breakFilteredDataIntoLaps(loops(i), RouteStructure);
    
    [aligned_Data_ByStation,mean_Data] = fcn_meanStationProjection(lapData,numLaps, RouteStructure);
    % obtain mean path data
    field_i = fields(i);
    x_mean.(field_i) = mean_Data.mean_xEast;
    y_mean.(field_i) = mean_Data.mean_yNorth;
    % calculate sigma
    sigma.(field_i) = fcn_calculateSigma(aligned_Data_ByStation,mean_Data,numLaps);
    
    % create field names for lapData (lap1, lap2, ...)
    lap_names = [];
    for k = 1:numLaps
        lap_names = [lap_names, string(['lap_',int2str(k)])];
    end
    loopsDataInLaps{i} = lapData;
end

% put into subplots?
%% Intersections 1 & 2
figure
[x_decision_12, y_decision_12] = fcn_calculateDecisionPoint(x_mean.loop_1, y_mean.loop_1, sigma.loop_1, x_mean.loop_2, y_mean.loop_2, sigma.loop_2, 2)

%% Intersections 3
[x_decision_3, y_decision_3] = fcn_calculateDecisionPoint(x_mean.loop_2, y_mean.loop_2, sigma.loop_2, x_mean.loop_3, y_mean.loop_3, sigma.loop_3)

%% Intersections 4 (cross road: straight & right)
[x_decision_4sr, y_decision_4sr] = fcn_calculateDecisionPoint(x_mean.loop_2, y_mean.loop_2, sigma.loop_2, x_mean.loop_4, y_mean.loop_4, sigma.loop_4)

%% Intersections 4 (cross road: straight & left)
[x_decision_4sl, y_decision_4sl] = fcn_calculateDecisionPoint(x_mean.loop_2, y_mean.loop_2, sigma.loop_2, x_mean.loop_5(500:900), y_mean.loop_5(500:900), sigma.loop_5(500:900))

%% Intersection 5 PROBLEMATIC
[x_decision_4sl, y_decision_4sl] = fcn_calculateDecisionPoint(x_mean.loop_2, y_mean.loop_2, sigma.loop_2, x_mean.loop_5(1000:end), y_mean.loop_5(1000:end), sigma.loop_5(1000:end))
%}