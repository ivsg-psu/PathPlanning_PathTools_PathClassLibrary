% Authors: Dr. Brennan, Guangwei, Liming
% Revision history
% 2020_05_20 - fixed bug on the yaw angle plots
% 2020_06_20 - add raw data query functions



%% Prep the workspace

% Clear the command window and workspace
clc
clear all %#ok<CLALL>

% Make sure we can see the utilities folder 
addpath '../Utilities';


% Add driver for database
%javaaddpath('C:\Users\Guangwei Zhou\AppData\Roaming\MathWorks\MATLAB\R2020a\drivers\postgresql-42.2.9.jar')

%%
flag.DBquery = false; %set to true if you want to query raw data from database insteading of loading from default *.mat file
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
        addpath '../Data';
        
        
        % Load the raw data
        % This data will have outliers, be unevenly sampled, have multiple and inconsistent measurements of the same variable.
        filename  = 'MappingVan_DecisionMaking_03132020.mat';
        variable_names = 'MappingVan_DecisionMaking_03132020';
        rawData = fcn_DataClean_loadRawData(filename,variable_names);

    end
    
    rawDataTimeFixed = fcn_DataClean_removeTimeGapsFromRawData(rawData);
end

% Data clean and merge
% Fill in the sigma values for key fields. This just calculates the sigma values for key fields (velocities,
% accelerations, angular rates in particular), useful for doing outlier detection, etc. in steps that follow.
rawDataWithSigmas = fcn_DataClean_loadSigmaValuesFromRawData(rawDataTimeFixed);


% NOTE: the following function changes the yaw angles to wind (correctly)
% up or down)

% Remove outliers on key fields via median filtering
% This removes outliers by median filtering key values.
rawDataWithSigmasAndMedianFiltered = fcn_DataClean_medianFilterFromRawAndSigmaData(rawDataWithSigmas);

% PLOTS to show winding up or down:
% figure(2); plot(mod(rawDataWithSigmas.GPS_Novatel.Yaw_deg,360),'b')
% figure(3); plot(mod(rawDataWithSigmasAndMedianFiltered.GPS_Novatel.Yaw_deg,360),'k')
% figure(4); plot(rawDataWithSigmasAndMedianFiltered.GPS_Novatel.Yaw_deg,'r')

% Clean the raw data
cleanData = fcn_DataClean_cleanRawDataBeforeTimeAlignment(rawDataWithSigmasAndMedianFiltered);

% Align all time vectors, and make time a "sensor" field
cleanAndTimeAlignedData = fcn_DataClean_alignToGPSTimeAllData(cleanData);

% Time filter the signals
timeFilteredData = fcn_DataClean_timeFilterData(cleanAndTimeAlignedData);

% Calculate merged data via Baysian averaging across same state
mergedData = fcn_DataClean_mergeTimeAlignedData(timeFilteredData);

% Remove jumps from merged data caused by DGPS outages
mergedDataNoJumps = fcn_DataClean_removeDGPSJumpsFromMergedData(mergedData,rawData);

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
[x_kf,sigma_x] = fcn_DataClean_KFmergeStateAndStateDerivative(t_x1,x1,x1_Sigma,t_x1dot,x1dot,x1dot_Sigma,nameString);
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
[x_kf,sigma_x] = fcn_DataClean_KFmergeStateAndStateDerivative(t_x1,x1,x1_Sigma,t_x1dot,x1dot,x1dot_Sigma,nameString);
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
[x_kf,sigma_x] = fcn_DataClean_KFmergeStateAndStateDerivative(t_x1,x1,x1_Sigma,t_x1dot,x1dot,x1dot_Sigma,nameString);
mergedByKFData.MergedGPS.yNorth = x_kf;
mergedByKFData.MergedGPS.yNorth_Sigma = sigma_x;





if 1==0
    % The following shows that we should NOT use yaw angles to calculate yaw rate
    fcn_plotArtificialYawRateFromYaw(MergedData,timeFilteredData);
    
    % Now to check to see if raw integration of YawRate can recover the yaw
    % angle
    fcn_plotArtificialYawFromYawRate(MergedData,timeFilteredData);
    
    
    %fcn_plotArtificialVelocityFromXAccel(MergedData,timeFilteredData);
    fcn_plotArtificialPositionFromIncrementsAndVelocity(MergedData,cleanAndTimeAlignedData)
end


%% Update plotting flags to allow merged data to now appear hereafter

clear plottingFlags
plottingFlags.fields_to_plot = [...
    %     {'All_AllSensors_velMagnitude'}...
    %     {'All_AllSensors_ZGyro'},...
    %     {'All_AllSensors_yNorth_increments'}...
    %     {'All_AllSensors_xEast_increments'}...
    %     {'All_AllSensors_xEast'}...
    %     {'All_AllSensors_yNorth'}...
    %     {'xEast'}...
    %     {'yNorth'}...
    %     {'xEast_increments'}...
    %     {'yNorth_increments'}...
    %     {'All_AllSensors_Yaw_deg'},...
    %     {'Yaw_deg'},...
    %     {'ZGyro_merged'},...
    %     {'All_AllSensors_ZGyro_merged'},...
    {'XYplot'},...
    %     {'All_AllSensors_XYplot'},...
    ];
% Define what is plotted
plottingFlags.flag_plot_Garmin = 0;

%
% % THE TEMPLATE FOR ALL PLOTTING
% fieldOrdering = [...
%     {'Yaw_deg'},...                 % Yaw variables
%     {'Yaw_deg_from_position'},...
%     {'Yaw_deg_from_velocity'},...
%     {'All_SingleSensor_Yaw_deg'},...
%     {'All_AllSensors_Yaw_deg'},...
%     {'Yaw_deg_merged'},...
%     {'All_AllSensors_Yaw_deg_merged'},...
%     {'ZGyro'},...                   % Yawrate (ZGyro) variables
%     {'All_AllSensors_ZGyro'},...
%     {'velMagnitude'},...            % velMagnitude variables
%     {'All_AllSensors_velMagnitude'},...
%     {'XAccel'},...                  % XAccel variables
%     {'All_AllSensors_XAccel'},...
%     {'xEast_increments'},...        % Position increment variables
%     {'All_AllSensors_xEast_increments'},...
%     {'yNorth_increments'},...
%     {'All_AllSensors_yNorth_increments'},...
%     {'zUp_increments'},...
%     {'All_AllSensors_zUp_increments'},...
%     {'XYplot'},...                  % XY plots
%     {'All_AllSensors_XYplot'},...
%     {'xEast'},...                   % xEast and yNorth plots
%     {'All_AllSensors_xEast'},...
%     {'yNorth'},...
%     {'All_AllSensors_yNorth'},...
%     {'zUp'},...
%     {'All_AllSensors_zUp'},...
%     {'DGPS_is_active'},...
%     {'All_AllSensors_DGPS_is_active'},...
%     %     {'velNorth'},...                % Remaining are not yet plotted - just kept here for now as  placeholders
%     %     {'velEast'},...
%     %     {'velUp'},...
%     %     {'Roll_deg'},...
%     %     {'Pitch_deg'},...
%     %     {'xy_increments'}... % Confirmed
%     %     {'YAccel'},...
%     %     {'ZAccel'},...
%     %     {'XGyro'},...
%     %     {'YGyro'},...
%     %     {'VelocityR},...
%     ];

% Define which sensors to plot individually
plottingFlags.SensorsToPlotIndividually = [...
    %    {'GPS_Hemisphere'}...
    %    {'GPS_Novatel'}...
    {'MergedGPS'}...
    %    {'VelocityProjectedByYaw'}...
    %     {'GPS_Garmin'}...
    %     {'IMU_Novatel'}...
    %     {'IMU_ADIS'}...
    %     {'Input_Steering'}...
    %     {'Encoder_RearWheels'}...
    %     {'MergedIMU'}...   
]; 

% Define zoom points for plotting
% plottingFlags.XYZoomPoint = [-4426.14413504648 -4215.78947791467 1601.69022519862 1709.39208889317]; % This is the corner after Toftrees, where the DGPS lock is nearly always bad
% plottingFlags.TimeZoomPoint = [297.977909295872          418.685505549775];
% plottingFlags.TimeZoomPoint = [1434.33632953011          1441.17612419014];
% plottingFlags.TimeZoomPoint = [1380   1600];
% plottingFlags.TimeZoomPoint = [760 840];
% plottingFlags.TimeZoomPoint = [596 603];  % Shows a glitch in the Yaw_deg_all_sensors plot
% plottingFlags.TimeZoomPoint = [1360   1430];  % Shows lots of noise in the individual Yaw signals
% plottingFlags.TimeZoomPoint = [1226 1233]; % This is the point of time discontinuity in the raw dat for Hemisphere
% plottingFlags.TimeZoomPoint = [580 615];  % Shows a glitch in xEast_increments plot
% plottingFlags.TimeZoomPoint = [2110 2160]; % This is the point of discontinuity in xEast
% plottingFlags.TimeZoomPoint = [2119 2129]; % This is the location of a discontinuity produced by a variance change
% plottingFlags.TimeZoomPoint = [120 150]; % Strange jump in xEast data
plottingFlags.TimeZoomPoint = [185 185+30]; % Strange jump in xEast data


% if isfield(plottingFlags,'TimeZoomPoint')
%     plottingFlags = rmfield(plottingFlags,'TimeZoomPoint');
% end
 
 
 % These set common y limits on values
% plottingFlags.ylim.('xEast') = [-4500 500];
plottingFlags.ylim.('yNorth') = [500 2500];

plottingFlags.ylim.('xEast_increments') = [-1.5 1.5];
plottingFlags.ylim.('All_AllSensors_xEast_increments') = [-1.5 1.5]; 
plottingFlags.ylim.('yNorth_increments') = [-1.5 1.5];
plottingFlags.ylim.('All_AllSensors_yNorth_increments') = [-1.5 1.5]; 

plottingFlags.ylim.('velMagnitude') = [-5 35]; 
plottingFlags.ylim.('All_AllSensors_velMagnitude') = [-5 35]; 


plottingFlags.PlotDataDots = 0; % If set to 1, then the data is plotted as dots as well as lines. Useful to see data drops.

%% Plot the results
%fcn_plotStructureData(rawData,plottingFlags);
%fcn_plotStructureData(rawDataTimeFixed,plottingFlags);
%fcn_plotStructureData(rawDataWithSigmas,plottingFlags);
%fcn_plotStructureData(rawDataWithSigmasAndMedianFiltered,plottingFlags);
%fcn_plotStructureData(cleanData,plottingFlags);
%fcn_plotStructureData(cleanAndTimeAlignedData,plottingFlags);
%fcn_plotStructureData(timeFilteredData,plottingFlags);
%fcn_plotStructureData(mergedData,plottingFlags);
fcn_plotStructureData(mergedDataNoJumps,plottingFlags);
%fcn_plotStructureData(mergedByKFData,plottingFlags);

% The following function allows similar plots, made when there are repeated
% uncommented versions above, to all scroll/zoom in unison.
%fcn_plotAxesLinkedTogetherByField;


%% Export results to Google Earth?
%fcn_exportXYZ_to_GoogleKML(rawData.GPS_Hemisphere,'rawData_GPS_Hemisphere.kml');
%fcn_exportXYZ_to_GoogleKML(mergedData.MergedGPS,'mergedData_MergedGPS.kml');
fcn_exportXYZ_to_GoogleKML(mergedDataNoJumps.MergedGPS,[dir.datafiles 'mergedDataNoJumps_MergedGPS.kml']);


%% Save cleaned data to .mat file 
% Liming stopped here
newStr = regexprep(trip_name{1},'\s','_'); % replace whitespace with underscore
newStr = strrep(newStr,'-','_');
cleaned_fileName = [newStr,'_cleaned'];
eval([cleaned_fileName,'=mergedByKFData'])
save(strcat(dir.datafiles,cleaned_fileName,'.mat'),cleaned_fileName)

