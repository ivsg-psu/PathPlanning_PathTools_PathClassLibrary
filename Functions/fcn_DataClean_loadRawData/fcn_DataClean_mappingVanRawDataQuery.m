%%%%%%%%%%%%%%%%%%%%%  Function fcn_mappingVanRawQuery %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%      query data from database mapping_van_raw
%
% Input Variables:
%      queryCondition = query method,format: string
%      queryFlag = same flag parameters for any query, format: struct
%
% Returned Results:
%      result = queried results data, format:struct
%      varargout = Variable-length output argument list
% Example:
%      [result,trip_name,trip_ids] = fcn_mappingVanRawQuery(queryCondition,queryFlag);
%
% Processing Flow:
%      IVSG yaw defintion: north is zero, clockwise is positive direction, range 0-360 degrees
%
% Restrictions/Notes:
% For more details about setup the Matlab-database , check https://github.com/ivsg-psu/Databases_Connections_ConnectMATLABToPostgreSQLServer/blob/master/README_to_Start.txt
% 
% The following functions are called:
%      This function has two class dependencies: MapDatabase,and Database
%
% Author:             Liming Gao, winstonglm@gamil.com
% Created Date:       2020-10-19:
% 
% Revisions:
%
% To do list:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [result,varargout] = fcn_DataClean_mappingVanRawDataQuery(database_name,queryCondition,queryFlag)

%---- Step1 :CONNECT TO  DATABASE ------------------------ %
% 
% javaaddpath('~\MATLAB\R2020a\drivers\postgresql-42.2.9.jar') % Add driver
% for database, this is not necessary if you setup the matlab-databse connector correctly. 

% choose different database name to connect to them
% Connect raw data database
MDB = MapDatabase(database_name); % instance of MapDatabase Class
% MDB.db.db_connection % show connection details

% show all tables schema in the database
% tables = MDB.db.ShowTables();

%---- Step2 : queryby method vrification

% method may be any of {'trip', 'bag', 'driver', 'date'} any simple contraction thereof.
valid = {'trip', 'date', 'bag_id', 'driver'};
[queryCondition,errstr] = validstring(queryCondition,valid);
if ~isempty(errstr)
    error('INTERPARC:incorrectmethod',errstr)
end

%---- Step3 : query data by specific condition

% show all trips and pick the ones you want to query --------- %

if strcmp(queryCondition(1:4),'trip') % query by trip
    %step1: set query parameters
    MDB.zero_time = queryFlag.zero_time;  % offset the timestamps starting from zero
    MDB.verbose = queryFlag.verbose;  % show the code processing details
    MDB.convert_GPS_to_ENU = queryFlag.convert_GPS_to_ENU; % you can choose the reference point by assign value to options.ENU_ref
    MDB.separate_by_lap = queryFlag.separate_by_lap; % default is 0
    
    % step2: pick trips you want to query(choose from trips table manually)
    % check all the trips files, each bag is trip here
    trips = MDB.fetchTrips();
    bagfiles = MDB.fetchBagFileIDs();
    trip_names = unique(trips.name,'stable'); % find trips
    
    % find daterange for each trip
    trip_withDate = trip_names;
    for i_trip= 1:length(trip_names)
        rows = strcmp(trips.name, trip_names{i_trip});
        id_trips = trips.id(rows); %find trip_ids for each trips
        bags_ofTrip = bagfiles(ismember(bagfiles.trips_id, id_trips),:); %find bagfiles_ids for each trips
        
        trip_startTime = bags_ofTrip.datetime(1); % this works because the data were sorted by datetime
        trip_endTime = bags_ofTrip.datetime_end(end);
        trip_withDate{i_trip} = [' ' num2str(i_trip) '.TripName:' trip_names{i_trip},'. Start_time:',trip_startTime{1},'. End_time:',trip_endTime{1}];
    end
    
    prompt_choose_trip = ['Which trip are you going to process? ' newline strjoin(trip_withDate,'\n') newline 'Input trip number: '];
    
    User_input = input(prompt_choose_trip);
    
    % pick the trip name
    trip_name = trip_names(User_input);
    
    prompt = ['Are you going to query the data of trip : ' trip_name{1} '?[y/n]'];
    User_input = input(prompt,'s');
    
    if strcmpi(User_input,'y')
        fprintf(1,'Thanks. Let''s query it...\n');
    else
        fprintf(1,'Query is aborted. \nYou can Re-pick the trips name.\n');
        return
    end
    
    % Step3: query data by trips ------------------------ %
    % query trips id according to the trip name
    trip_id = [];
    for i = 1:length(trip_name)
        trip_id = cat(1,trip_id,trips.id(strcmp(trips.name, trip_name{i})));
    end
    
    % Pick sensors you want to query. 1 means query data from that sensor
    options = {};
    options.sensors.base_station = queryFlag.sensors.base_station;  % default is 1
    options.sensors.hemisphere_gps = queryFlag.sensors.hemisphere_gps; % default is 1
    options.sensors.NovAtel_gps = queryFlag.sensors.NovAtel_gps; % default is 1
    options.sensors.garmin_gps = queryFlag.sensors.garmin_gps; % default is 1
    options.sensors.garmin_velocity = queryFlag.sensors.garmin_velocity; % default is 0
    options.sensors.steering_angle =queryFlag.sensors.steering_angle; % default is 1
    options.sensors.NovAtel_imu = queryFlag.sensors.NovAtel_imu;% default is 1
    options.sensors.adis_imu = queryFlag.sensors.adis_imu;% default is 1
    options.sensors.encoder_left_right = queryFlag.sensors.encoder_left_right;% default is 1
    options.sensors.laser = queryFlag.sensors.laser; % default is 0
    options.sensors.front_left_camera = queryFlag.sensors.front_left_camera; % default is 0
    options.sensors.front_center_camera = queryFlag.sensors.front_center_camera; % default is 0
    options.sensors.front_right_camera = queryFlag.sensors.front_right_camera; % default is 0
    
    options.ENU_ref = queryFlag.ENU_ref; % 0 use default setting in database, 1 test track, 2 LTI, Larson  Transportation Institute
    
    % fetchByTripID
    result = MDB.fetchByTripID(trip_id,options);
    
    %
    varargout{1} = trip_name; % output trip name
    varargout{2} = trip_id; % output trip id
    % disconnect with DB
    MDB.disconnect();
    
elseif strcmp(queryCondition(1:4),'date') % query by 'date', 'bag_id', 'driver'
    %step1: set query parameters
    MDB.zero_time = queryFlag.zero_time;  % offset the timestamps starting from zero
    MDB.verbose = queryFlag.verbose;  % show the code processing details
    MDB.convert_GPS_to_ENU = queryFlag.convert_GPS_to_ENU; % you can choose the reference point by assign value to options.ENU_ref
    MDB.separate_by_lap = queryFlag.separate_by_lap; % default is 0
    
    %step2: check all the trips and date
    trips = MDB.fetchTrips();
    bagfiles = MDB.fetchBagFileIDs();
    trip_names = unique(trips.name,'stable'); % find trips
    
    % find daterange for each trip
    trip_withDate = trip_names;
    for i_trip= 1:length(trip_names)
        rows = strcmp(trips.name, trip_names{i_trip});
        id_trips = trips.id(rows); %find trip_ids for each trips
        bags_ofTrip = bagfiles(ismember(bagfiles.trips_id, id_trips),:); %find bagfiles_ids for each trips
        
        trip_startTime = bags_ofTrip.datetime(1); % this works because the data were sorted by datetime
        trip_endTime = bags_ofTrip.datetime_end(end);
        trip_withDate{i_trip} = [' ' num2str(i_trip) '.TripName:' trip_names{i_trip},'. Start_time:',trip_startTime{1},'. End_time:',trip_endTime{1}];
    end
    
    %%step3: pick time range you want to query
    %%show available trips and time range
    fprintf([newline 'Available trips and time range:' newline strjoin(trip_withDate,'\n') newline]);
    
    % input start and end time
    prompt_InputStartTime =[newline 'Input Strat Time(Format: yyyy-mm-dd HH:MM:ss): '];
    User_InputStartTime = input(prompt_InputStartTime,'s');
    
    prompt_InputEndTime =[newline 'Input End Time(Format: yyyy-mm-dd HH:MM:ss): '];
    User_InputEndTime = input(prompt_InputEndTime,'s');
    
    prompt = [newline 'Are you going to query the data from : ' User_InputStartTime ' to ' User_InputEndTime '?[y/n]'];
    User_input = input(prompt,'s');
    
    % dateformat
    try
        datetime(User_InputStartTime,'InputFormat','yyyy-MM-dd HH:mm:ss');
        datetime(User_InputEndTime,'InputFormat','yyyy-MM-dd HH:mm:ss');
    catch
        fprintf(1,'Your Input date format is wrong.Please try again. \nQuery is aborted..\n');
        return
    end
    % confrim the information with user
    if strcmpi(User_input,'y')
        fprintf(1,'Thanks. Let''s query it...\n');
    else
        fprintf(1,'Query is aborted. \nYou can Re-pick the time range.\n');
        return
    end
    
    % Step4: query data by time range ------------------------ %
    
    % Pick sensors you want to query. 1 means query data from that sensor
    options = {};
    options.sensors.base_station = queryFlag.sensors.base_station;  % default is 1
    options.sensors.hemisphere_gps = queryFlag.sensors.hemisphere_gps; % default is 1
    options.sensors.NovAtel_gps = queryFlag.sensors.NovAtel_gps; % default is 1
    options.sensors.garmin_gps = queryFlag.sensors.garmin_gps; % default is 1
    options.sensors.garmin_velocity = queryFlag.sensors.garmin_velocity; % default is 0
    options.sensors.steering_angle =queryFlag.sensors.steering_angle; % default is 1
    options.sensors.NovAtel_imu = queryFlag.sensors.NovAtel_imu;% default is 1
    options.sensors.adis_imu = queryFlag.sensors.adis_imu;% default is 1
    options.sensors.encoder_left_right = queryFlag.sensors.encoder_left_right;% default is 1
    options.sensors.laser = queryFlag.sensors.laser; % default is 0
    options.sensors.front_left_camera = queryFlag.sensors.front_left_camera; % default is 0
    options.sensors.front_center_camera = queryFlag.sensors.front_center_camera; % default is 0
    options.sensors.front_right_camera = queryFlag.sensors.front_right_camera; % default is 0
    
    options.ENU_ref = queryFlag.ENU_ref; % 0 use default setting in database, 1 test track, 2 LTI, Larson  Transportation Institute
    
    % fetchByTimeRange
    result = MDB.fetchByTimeRange(User_InputStartTime,User_InputEndTime,options); %result of query, format is struct
    
    varargout{1} = User_InputStartTime; % output User_InputStartTime
    varargout{2} = User_InputEndTime; % output User_InputEndTime
    % disconnect with DB
    MDB.disconnect();
    
elseif strcmp(queryCondition(1:6),'driver') % query by driver
    % Note: need more edit
    
elseif strcmp(queryCondition(1:6),'location') % query by location
    % Note: need more edit
    
else
    error('no query condition');
    
    
end

end

% ===============================================
%  nested function for string vrification
% ===============================================
function [str,errorclass] = validstring(arg,valid)
% validstring: compares a string against a set of valid options
% usage: [str,errorclass] = validstring(arg,valid)
%
% If a direct hit, or any unambiguous shortening is found, that
% string is returned. Capitalization is ignored.
%
% arguments: (input)
%  arg - character string, to be tested against a list
%        of valid choices. Capitalization is ignored.
%
%  valid - cellstring array of alternative choices
%
% Arguments: (output)
%  str - string - resulting choice resolved from the
%        list of valid arguments. If no unambiguous
%        choice can be resolved, then str will be empty.
%
%  errorclass - string - A string argument that explains
%        the error. It will be one of the following
%        possibilities:
%
%        ''  --> No error. An unambiguous match for arg
%                was found among the choices.
%
%        'No match found' --> No match was found among
%                the choices provided in valid.
%
%        'Ambiguous argument' --> At least two ambiguous
%                matches were found among those provided
%                in valid.
%
%
% Example:
%  valid = {'off' 'on' 'The sky is falling'}
%
%
% See also: parse_pv_pairs, strmatch, strcmpi
%

ind = find(strncmpi(arg,valid,numel(arg)));
if isempty(ind)
    % No hit found
    errorclass = 'No match found';
    str = '';
elseif (length(ind) > 1)
    % Ambiguous arg, hitting more than one of the valid options
    errorclass = 'Ambiguous argument';
    str = '';
    return
else
    errorclass = '';
    str = valid{ind};
end

end % function validstring
