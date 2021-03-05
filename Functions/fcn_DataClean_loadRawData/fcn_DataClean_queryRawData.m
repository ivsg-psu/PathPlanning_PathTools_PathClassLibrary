%%%%%%%%%%%%%%%%%%%%%  Function fcn_DataClean_queryRawData %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%      query raw data from a *.mat file or from mapping_van_raw database
%
% Input Variables:
%      DBqueryFlag = if we query data from database,format: true or false 
%      if true(query from batabase) 
%         varargin{1} = databaseName, format: string
%         varargin{2} = queryCondition, format: string
%      if false(load from *.mat file) 
%         varargin{1} = filename, format:string
%         varargin{2} = variable_names,format:string
%
% Returned Results:
%      rawData = structured data format, format:struct
%      varargout = Variable-length output argument list
% Example:
%      1. query from database 
%         rawData = fcn_DataClean_queryRawData(flag.DBquery,'mapping_van_raw','trip'); % more query condition can be set in the function 
%      2. load from file
%         filename  = 'MappingVan_DecisionMaking_03132020.mat';
%         variable_names = 'MappingVan_DecisionMaking_03132020';
%         rawData = fcn_DataClean_queryRawData(flag.DBquery,filename,variable_names); % more query condition can be set in the function 
%
% Processing Flow:
%      IVSG yaw defintion: north is zero, clockwise is positive direction, range 0-360 degrees
%
% Restrictions/Notes:
%
% The following functions are called:
%      This function has two class dependencies: MapDatabase,and Database
%
% Author:             Liming Gao, winstonglm@gamil.com
% Created Date:       2021-01-08:
% Revisions:
%
% To do list:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rawData,varargout] = fcn_DataClean_queryRawData(DBqueryFlag,varargin)


if DBqueryFlag == true
    
    database_name = varargin{1};
    % set some parameters before running the query
    % (1). query condition
    queryCondition = varargin{2};
    %queryCondition = 'trip'; % can be 'trip', 'date', or 'driver'
    % (2). data pre-processing
    queryFlag.zero_time = 0;  % offset the timestamps starting from zero
    queryFlag.verbose = 1;  % show the code processing details
    queryFlag.convert_GPS_to_ENU = 1; % you can choose the reference point by assign value to queryFlag.ENU_ref
    queryFlag.separate_by_lap = 0; % default is 0, do not set this to true for the data without laps
    % (3). choose the sensor you want to query
    queryFlag.ENU_ref = 0; % 0:use default setting in database, 1:test track, 2: LTI, Larson  Transportation Institute
    queryFlag.sensors.base_station = 1;  % default is 1
    queryFlag.sensors.hemisphere_gps = 1; % default is 1
    queryFlag.sensors.NovAtel_gps = 1; % default is 1
    queryFlag.sensors.garmin_gps = 1; % default is 1
    queryFlag.sensors.garmin_velocity = 1; % default is 0
    queryFlag.sensors.steering_angle =1; % default is 1
    queryFlag.sensors.NovAtel_imu = 1;% default is 1
    queryFlag.sensors.adis_imu = 1;% default is 1
    queryFlag.sensors.encoder_left_right = 1;% default is 1
    queryFlag.sensors.laser = 0; % default is 0
    queryFlag.sensors.front_left_camera = 0; % default is 0
    queryFlag.sensors.front_center_camera = 0; % default is 0
    queryFlag.sensors.front_right_camera = 0; % default is 0
    
    % run query function
    [result,trip_name,~,trip_id_cleaned] = fcn_DataClean_mappingVanRawDataQuery(database_name,queryCondition,queryFlag);
    
    % load and pre-process the data
    % rawData = fcn_preProcessQueryResult(result); % no longer use
    flag_do_debug = 1;
    rawData = fcn_DataClean_loadRawData(flag_do_debug,result); %% try to use the anonymous function https://www.youtube.com/watch?v=kE4SUA392_I&feature=youtu.be
    varargout{1} = trip_name; % output trip name
    varargout{2} = trip_id_cleaned; % trip id will be used at the data database  
    varargout{3} = result.base_station; % used for LLA and ENU conversion
else
    
    % Load the raw data from file
    % This data will have outliers, be unevenly sampled, have multiple and inconsistent measurements of the same variable.
    
    filename = varargin{1};
    variable_names = varargin{2};
    % load and pre-process the data
    % rawData = fcn_loadRawData(filename,variable_names); % no longer use
    flag_do_debug = 1;
    rawData = fcn_DataClean_loadRawData(flag_do_debug,filename,variable_names); %
    
end

end

