%%%%%%%%%%%%%%%%%%%%%  class MapDatabase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%      wrapper class of query functions
%
% Input Variables:
%      database_name = the name of database you want to connect,format:string, eg:database_name = 'mapping_van_raw';
%
% Returned Results:
%
% Example:
% m = MapDatabase(database_name);
%
% Processing Flow:
%
%
% Restrictions/Notes:
%   need Databse class
%
% The following functions(methods) are included:
%      obj = MapDatabase(database_name, ip_address, port, image_directory, username, password)
%      function result = fetchTrips(obj)
%
% Author:             Liming Gao
% Created Date:       2020-02-28
% Revisions:
%           2020-02-29:
%
% To do list:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MapDatabase < handle
    
    % Class properties and variables
    properties
        
        % Parameters that the user must define when instantiating the class
        image_directory
        database_name
        
        % General parameters
        zero_time = 0;
        verbose = 0;
        interpolate = 0; % Not used right now
        convert_GPS_to_ENU = 1;
        separate_by_lap = 0;
        
        % Default settings
        ip_address = '130.203.223.234' %  'localhost'; %Ip address of server host
        port = 5432;  % port number
        username = 'ivsg_db_user'   % user name for the server
        password = 'ivsg@DB320' % password
        
        spheroid = referenceEllipsoid('wgs84');
        %instance of DB
        db
        
    end
    
    properties (Access = private)
        
        minimum_time
        
    end
    
    methods
        
        % Constructor for class. This uses the internal connect function to
        % create a connection to the database.
        function obj = MapDatabase(database_name, ip_address, port, image_directory, username, password)
            
            if nargin == 1
                obj.database_name = database_name;
                obj.db = Database(database_name);
                
            elseif nargin == 5
                obj.database_name = database_name;
                obj.image_directory = image_directory;
                obj.ip_address = ip_address;
                obj.port = port;
                obj.username = username;
                obj.password = password;
                
                obj.db = Database(database_name, ip_address, port, username, password);
                
            else
                
                error('MapDatabase requires 6 or 1 inputs: database name,image_directory, ip address, port number,  database username and password');
                
            end
            
        end
        
        % =============================== function fetchSensor===============================
        % purpose:      query data from specific table
        % intput:       table = table name, format: string
        %               fields = fields, format: cell array
        %               where = where, format: cell array
        %               orderby = orderby, format: cell array
        % output:       result = struct format of query result with table format element
        function result = fetchSensor(obj, table, fields, where, orderby)
            
            [result,~, column_names] = obj.db.select(table, fields, where, orderby);
            
            result = obj.db.convertFromTableToStruct(column_names, result);
            
            if strcmp(table, 'laser')
                
                sscanf_format = repmat('%f',1,1141);
                for i = 1:size(result.time,1)
                    result.ranges{i} = sscanf(result.ranges{i}, sscanf_format)';
                    result.intensities{i} = sscanf(result.intensities{i}, sscanf_format)';
                end
                
            end
            
        end
        
        % =============================== function fetchAll===============================
        % purpose:      query data from BD given trip_id
        % intput:       trip_id = trip id, format: numeric array
        % output:       result_table = table format of query result
        %               result_struct = struct format of query result with array format element
        %               result = struct format of query result with table format element
        % ===========================================================================
        function result = fetchAll(obj, where, options)
            
            sensor_pose_trip_id = 1;
            
            result = {};
            
            orderby = 'time'; % sort the query by time
            
            if obj.verbose == 1
                fprintf(1,'Start loading sensors data... \n')
                fprintf(1,'----------------------- \n')
                tic
            end
            
            % 0. Base Station
            if nargin == 2 || (nargin == 3 && isfield(options.sensors,'base_station') && options.sensors.base_station == 1)
                
                if obj.verbose == 1
                    fprintf('Loading base station data...\n')
                end
                
                % determine base station id (ENU reference id)
                if options.ENU_ref == 0
                    if strncmp(where{1},'bag_files_id in',15)
                        sqlquery =['SELECT DISTINCT base_stations_id FROM trips where id in ( SELECT trips_id FROM bag_files where '...
                            replace(where{1},'bag_files_id','id')  ');']; % find base_station_id given bagfilesid, be carefule with the space
                        base_station_id = fetch(obj.db.db_connection,sqlquery);
                        id = base_station_id.base_stations_id;
                        
                        if length(id) > 1 % if More than one base stations are found, there are some issue in the LLA ENU transformation
                            warning('More than one base stations are found!');
                        end
                        sqlquery =['SELECT name FROM base_stations where id in (' num2str(id)  ');']; % find base station name given id
                        base_station_id = fetch(obj.db.db_connection,sqlquery);
                        fprintf('Base Station: %s\n', base_station_id.name{:})
                    else
                        id = 1;
                        fprintf('Base Station: Test Track\n')
                    end
                elseif options.ENU_ref == 1
                    id = 1;
                    fprintf('Base Station: Test Track\n')
                elseif options.ENU_ref == 2
                    id = 2;
                    fprintf('Base Station: LTI, Larson  Transportation Institute\n')
                else
                    error('No base station')
                end
                
                table = 'base_stations';
                fields = {'id', 'latitude', 'longitude', 'altitude', 'latitude_std', 'longitude_std', 'altitude_std'};
                where_base_station = {cat(2,'id = ', num2str(id))};
                result_base_station = fetchSensor(obj, table, fields, where_base_station, '');
                result.base_station = result_base_station;
                
                % REPLACE it with matlab built-in function
                % g = GPS(result.base_station.latitude, result.base_station.longitude, result.base_station.altitude);
                if obj.verbose == 1
                    fprintf(1,'Load base station data Done! \n\n')
                end
            end
            
            % 1. GPS (Hemisphere)
            if nargin == 2 || (nargin == 3 && isfield(options.sensors,'hemisphere_gps') && options.sensors.hemisphere_gps == 1)
                
                if obj.verbose == 1
                    fprintf('Loading GPS (Hemisphere) data...\n')
                end
                
                % sensor_id = 0;
                table = 'hemisphere_gps';
                %{'id','sensors_id','bag_files_id','navmode','latitude','longitude','altitude','geography','vnorth','veast','vup',
                % 'stddevresid','roll','pitch','yaw','gps_week','gps_seconds','numofsats','manualmark','ageofdiff',
                % 'extendedageofdiff','seconds','nanoseconds','time','timestamp','date_added'}
                fields = {'id', 'ageofdiff','gps_week','gps_seconds','stddevresid','extendedageofdiff','time',...
                    'latitude', 'longitude', 'altitude','navmode', 'vnorth','veast','vup','numofsats','manualmark'};
                result_gps = fetchSensor(obj, table, fields, where, orderby);
                result.Hemisphere_DGPS = result_gps;
                
                % result_pose = fetchSensorPose(obj, sensor_pose_trip_id, sensor_id);
                % result.gps.pose = result_pose;
                
                if obj.convert_GPS_to_ENU == 1 && (nargin == 3 && isfield(options.sensors,'base_station') && options.sensors.base_station == 1)
                    
                    [xEast, yNorth,zUp] = geodetic2enu(result.Hemisphere_DGPS.latitude,result.Hemisphere_DGPS.longitude ,result.Hemisphere_DGPS.altitude,...
                        result.base_station.latitude, result.base_station.longitude, result.base_station.altitude,obj.spheroid);
                    
                    % ENU = g.WGSLLA2ENU(result.gps.latitude, result.gps.longitude, result.gps.altitude);
                    station = obj.calculateStation(xEast, yNorth,zUp);
                    
                    result.Hemisphere_DGPS.xeast = xEast;
                    result.Hemisphere_DGPS.ynorth = yNorth;
                    result.Hemisphere_DGPS.zup = zUp;
                    result.Hemisphere_DGPS.station = station;
                    
                end
                
                if obj.verbose == 1
                    fprintf(1,'Load GPS (Hemisphere) data Done! \n\n')
                end
                
            end
            
            % 2. GPS (Novatel)
            if nargin == 2 || (nargin == 3 && isfield(options.sensors,'NovAtel_gps') && options.sensors.NovAtel_gps == 1)
                
                if obj.verbose == 1
                    fprintf('Loading GPS (Novatel) data...\n')
                end
                
                % sensor_id = 0;
                table = 'novatel_gps';
                % {'id','sensors_id','bag_files_id','status','latitude','longitude','altitude','geography','east_velocity','north_velocity','up_velocity','roll','pitch','yaw','gps_week','gps_seconds','seconds','nanoseconds','time','timestamp','date_added'}
                fields = {'id', 'latitude', 'longitude', 'altitude', 'gps_week','gps_seconds','time','east_velocity','north_velocity','up_velocity','roll', 'pitch', 'yaw','status'};
                result_gps = fetchSensor(obj, table, fields, where, orderby);
                result.GPS_Novatel = result_gps;
                
                % result_pose = fetchSensorPose(obj, sensor_pose_trip_id, sensor_id);
                % result.gps.pose = result_pose;
                
                % Calculate lateral and longitudinal velocity from INS
                result.GPS_Novatel.U = result.GPS_Novatel.east_velocity .* cos(result.GPS_Novatel.yaw) + result.GPS_Novatel.north_velocity .* sin(result.GPS_Novatel.yaw);
                result.GPS_Novatel.V = -result.GPS_Novatel.east_velocity .* sin(result.GPS_Novatel.yaw) + result.GPS_Novatel.north_velocity .* cos(result.GPS_Novatel.yaw);
                
                if obj.convert_GPS_to_ENU == 1 && (nargin == 3 && isfield(options.sensors,'base_station') && options.sensors.base_station == 1)
                    
                    [xEast, yNorth,zUp] = geodetic2enu(result.GPS_Novatel.latitude,result.GPS_Novatel.longitude ,result.GPS_Novatel.altitude,...
                        result.base_station.latitude, result.base_station.longitude, result.base_station.altitude,obj.spheroid);
                    
                    % ENU = g.WGSLLA2ENU(result.gps.latitude, result.gps.longitude, result.gps.altitude);
                    station = obj.calculateStation(xEast, yNorth,zUp);
                    
                    result.GPS_Novatel.xEast = xEast;
                    result.GPS_Novatel.yNorth = yNorth;
                    result.GPS_Novatel.zUp = zUp;
                    result.GPS_Novatel.station = station;
                    
                end
                
                if obj.verbose == 1
                    fprintf(1,'Load GPS (Novatel) data Done! \n\n')
                end
            end
            
            % 3. GPS (Garmin)
            if nargin == 2 || (nargin == 3 && isfield(options.sensors,'garmin_gps') && options.sensors.garmin_gps == 1)
                
                if obj.verbose == 1
                    fprintf('Loading GPS (Garmin) data...\n')
                end
                
                %sensor_id = 8;
                table = 'garmin_gps';
                % {'id','sensors_id','bag_files_id','status','service','latitude','longitude','altitude','geography',
                %  'roll','pitch','yaw','seconds','nanoseconds','position_covariance','position_covariance_type','time','timestamp','date_added'}
                fields = {'id', 'time', 'latitude', 'longitude', 'altitude','position_covariance','position_covariance_type'};
                result_gps = fetchSensor(obj, table, fields, where, orderby);
                result.Garmin_GPS = result_gps;
                
                %                 result_pose = fetchSensorPose(obj, sensor_pose_trip_id, sensor_id);
                %                 result.Garmin_GPS.pose = result_pose;
                
                if obj.convert_GPS_to_ENU == 1 && (nargin == 3 && isfield(options.sensors,'base_station') && options.sensors.base_station == 1)
                    
                    [xEast, yNorth,zUp] = geodetic2enu(result.Garmin_GPS.latitude,result.Garmin_GPS.longitude ,result.Garmin_GPS.altitude,...
                        result.base_station.latitude,result.base_station.longitude,result.base_station.altitude,obj.spheroid);
                    station = obj.calculateStation(xEast, yNorth,zUp);
                    
                    result.Garmin_GPS.xEast = xEast;
                    result.Garmin_GPS.yNorth = yNorth;
                    result.Garmin_GPS.zUp = zUp;
                    result.Garmin_GPS.station = station;
                    
                end
                if obj.verbose == 1
                    fprintf(1,'Load GPS (Garmin) data Done! \n\n')
                end
            end
            
            % 4. Garmin Velocity
            if nargin == 2 || (nargin == 3 && isfield(options.sensors,'garmin_velocity') && options.sensors.garmin_velocity == 1)
                
                if obj.verbose == 1
                    fprintf('Loading Garmin velocity data...\n')
                end
                
                % sensor_id = 8;
                table = 'garmin_velocity';
                fields = {'id', 'time', 'latitude', 'longitude', 'altitude', 'east_velocity', 'north_velocity'};
                result_gps = fetchSensor(obj, table, fields, where, orderby);
                result.garmin_velocity = result_gps;
                
                %                 result_pose = fetchSensorPose(obj, sensor_pose_trip_id, sensor_id);
                %                 result.garmin_velocity.pose = result_pose;
                
                if obj.convert_GPS_to_ENU == 1 && (nargin == 3 && isfield(options.sensors,'base_station') && options.sensors.base_station == 1)
                    
                    [xEast, yNorth,zUp] = geodetic2enu(result.garmin_velocity.latitude,result.garmin_velocity.longitude ,result.garmin_velocity.altitude,...
                        result.base_station.latitude,result.base_station.longitude,result.base_station.altitude,obj.spheroid);
                    station = obj.calculateStation(xEast, yNorth,zUp);
                    
                    result.garmin_velocity.xEast = xEast;
                    result.garmin_velocity.yNorth = yNorth;
                    result.garmin_velocity.zUp = zUp;
                    result.garmin_velocity.station = station;
                    
                end
                
                if obj.verbose == 1
                    fprintf(1,'Load Garmin velocity data Done! \n\n')
                end
                
            end
            
            % 5. steering_angle
            if nargin == 2 || (nargin == 3 && isfield(options.sensors,'steering_angle') && options.sensors.steering_angle == 1)
                
                if obj.verbose == 1
                    fprintf('Loading Steering_angle data...\n')
                end
                
                % sensor_id = 1;
                table = 'steering_angle';
                % {'id','bag_files_id','sensors_id','left_counts','right_counts','left_counts_filtered','right_counts_filtered','left_angle','right_angle','angle',
                % 'latitude','longitude','altitude','geography','seconds','nanoseconds','time','timestamp','date_added'}
                fields = {'id', 'time', 'left_counts','right_counts','left_counts_filtered','right_counts_filtered','left_angle','right_angle','angle'};
                result_steering = fetchSensor(obj, table, fields, where, orderby);
                result.Steering_angle = result_steering;
                
                % result_pose = fetchSensorPose(obj, sensor_pose_trip_id, sensor_id);
                % result.imu.pose = result_pose;
                
                if obj.verbose == 1
                    fprintf('Load Steering_angle data done! \n\n')
                end
            end
            
            % 6. IMU (Novatel_IMU)
            if nargin == 2 || (nargin == 3 && isfield(options.sensors,'NovAtel_imu') && options.sensors.NovAtel_imu == 1)
                
                if obj.verbose == 1
                    fprintf('Loading IMU (Novatel) data...\n')
                end
                
                % sensor_id = 1;
                table = 'novatel_imu';
                % {'id','bag_files_id','sensors_id','status','x_acceleration','y_acceleration','z_acceleration','x_angular_velocity','y_angular_velocity','z_angular_velocity','latitude','longitude','altitude',
                % 'roll','pitch','yaw','geography','gps_week','gps_seconds','seconds','nanoseconds','time','timestamp','date_added'}
                fields = {'id', 'time', 'gps_week','gps_seconds','status','x_acceleration', 'y_acceleration', 'z_acceleration', 'x_angular_velocity', 'y_angular_velocity', 'z_angular_velocity'};
                result_imu = fetchSensor(obj, table, fields, where, orderby);
                result.Novatel_IMU = result_imu;
                
                %                 result_pose = fetchSensorPose(obj, sensor_pose_trip_id, sensor_id);
                %                 result.imu.pose = result_pose;
                
                if obj.verbose == 1
                    fprintf('Load IMU (Novatel) data done! \n\n')
                end
            end
            
            % 7. IMU (Adis)
            if nargin == 2 || (nargin == 3 && isfield(options.sensors,'adis_imu') && options.sensors.adis_imu == 1)
                
                if obj.verbose == 1
                    fprintf('Loading IMU (Adis) data...\n')
                end
                
                % sensor_id = 1;
                table = 'adis_imu';
                % {'id','sensors_id','bag_files_id','status','x_acceleration','y_acceleration','z_acceleration','x_angular_velocity','y_angular_velocity','z_angular_velocity',
                % 'geography','roll','pitch','yaw','magnetic_x','magnetic_y','magnetic_z','temperature','pressure','seconds','nanoseconds','time','timestamp','date_added'}
                fields = {'id', 'time', 'x_acceleration', 'y_acceleration', 'z_acceleration', 'x_angular_velocity', 'y_angular_velocity', 'z_angular_velocity'};
                result_imu = fetchSensor(obj, table, fields, where, orderby);
                result.adis_IMU_data = result_imu;
                
                % result_pose = fetchSensorPose(obj, sensor_pose_trip_id, sensor_id);
                % result.imu.pose = result_pose;
                
                table = 'adis16407_parameters';
                % {'id','sensors_id','accel_scale','gyro_scale','mag_scale','bar_scale','temp_scale','gravity','accel_cov','gyro_cov','mag_cov','bar_cov','date_added'}
                fields = {'id','accel_scale','gyro_scale','mag_scale','bar_scale','temp_scale','gravity','accel_cov','gyro_cov','mag_cov','bar_cov'};
                [result_imu_parameters,~, column_names] = obj.db.select(table, fields);
                result_imu_parameters = obj.db.convertFromTableToStruct(column_names, result_imu_parameters);
                result.adis_IMU_parameters = result_imu_parameters;
                
                if obj.verbose == 1
                    fprintf('Load IMU (Adis) data done! \n\n')
                end
            end
            
            % 8. Encoders
            if nargin == 2 || (nargin == 3 && isfield(options.sensors,'encoder_left_right') && options.sensors.encoder_left_right == 1)
                
                if obj.verbose == 1
                    fprintf('Loading encoders data...\n')
                end
                
                % Left encoder and right encoder
                % sensor_id = 6;
                table = 'encoder';
                % {'id','bag_files_id','sensors_id','revolutions','left_counts','right_counts','left_delta_counts','right_delta_counts','left_angular_velocity',
                % 'right_angular_velocity','latitude','longitude','altitude','geography','seconds','nanoseconds','time','timestamp','date_added'}
                fields = {'id', 'time','revolutions','left_counts','right_counts','left_delta_counts','right_delta_counts','left_angular_velocity','right_angular_velocity'};
                % where_encoder = {where, cat(2,'sensor_id = ',num2str(sensor_id))};
                where_encoder = where;
                result_encoder = fetchSensor(obj, table, fields, where_encoder, orderby);
                result.Raw_encoder = result_encoder;
                %               result.encoder_left.delta_counts = -result.encoder_left.delta_counts;
                %               result.encoder_left.angular_velocity = -result.encoder_left.angular_velocity;
                
                table = 'encoder_parameters';
                % {'id','sensors_id','counts_per_revolution','date_added'}
                fields = {'id','counts_per_revolution'};
                [result_encoder_parameters,~, column_names] = obj.db.select(table, fields);
                result_encoder_parameters = obj.db.convertFromTableToStruct(column_names, result_encoder_parameters);
                result.encoder_parameters = result_encoder_parameters;
                
                if obj.verbose == 1
                    fprintf(1,'Load encoders data Done! \n\n')
                end
                
            end
            
            
            %Notes:(continue edit from here )
            % Laser
            if nargin == 2 || (nargin == 3 && isfield(options.sensors,'laser') && options.sensors.laser == 1)
                
                if obj.verbose == 1
                    fprintf('Loading laser data...\n')
                end
                
                sensor_id = 2;
                table = 'laser';
                fields = {'id', 'time', 'scan_time', 'ranges', 'intensities', 'latitude', 'longitude', 'altitude'};
                where_laser = {where, cat(2,'sensor_id = ',num2str(sensor_id))};
                result_laser = fetchSensor(obj, table, fields, where_laser, orderby);
                result.laser = result_laser;
                
                result_pose = fetchSensorPose(obj, sensor_pose_trip_id, sensor_id);
                result.laser.pose = result_pose;
                
                table = 'laser_parameters';
                fields = {'id', 'angle_min', 'angle_max', 'angle_increment', 'time_increment', 'range_min', 'range_max'};
                where_laser = {cat(2,'sensor_id = ', num2str(sensor_id)), cat(2,'trip_id = ', num2str(trip_id))};
                result_laser_parameters = fetchSensor(obj, table, fields, where_laser, '');
                result.laser.parameters = result_laser_parameters;
                
                if obj.convert_GPS_to_ENU == 1 && (nargin == 3 && isfield(options.sensors,'base_station') && options.sensors.base_station == 1)
                    
                    ENU = g.WGSLLA2ENU(result.laser.latitude, result.laser.longitude, result.laser.altitude);
                    station = obj.calculateStation(ENU');
                    
                    result.laser.station = station;
                    result.laser.X = ENU(1,:)';
                    result.laser.Y = ENU(2,:)';
                    result.laser.Z = ENU(3,:)';
                    
                end
                
            end
            
            %%% Camera %%%
            
            % Rear left camera
            if nargin == 2 || (nargin == 3 && isfield(options.sensors,'rear_left_camera') && options.sensors.rear_left_camera == 1)
                
                if obj.verbose == 1
                    fprintf('Loading rear left camera data...\n')
                end
                
                sensor_id = 3;
                table = 'camera';
                fields = {'id', 'time', 'seconds_triggered', 'nanoseconds_triggered', 'file_name', 'latitude', 'longitude', 'altitude', 'roll', 'pitch', 'yaw'};
                where_camera = {where, cat(2,'sensor_id = ',num2str(sensor_id))};
                result_camera = fetchSensor(obj, table, fields, where_camera, orderby);
                result.rear_left_camera = result_camera;
                
                result.rear_left_camera.time_triggered = result.rear_left_camera.seconds_triggered + result.rear_left_camera.nanoseconds_triggered * 10 ^ (-9);
                result.rear_left_camera = rmfield(result.rear_left_camera, {'seconds_triggered', 'nanoseconds_triggered'});
                
                result_pose = fetchSensorPose(obj, sensor_pose_trip_id, sensor_id);
                result.rear_left_camera.pose = result_pose;
                
                table = 'camera_parameters';
                fields = {'id', 'focal_x', 'focal_y', 'center_x', 'center_y', 'skew', 'image_width', 'image_height', 'distortion_k1', 'distortion_k2', 'distortion_p1', 'distortion_p2', 'distortion_k3'};
                where_camera = {cat(2,'sensor_id = ', num2str(sensor_id)), cat(2,'trip_id = ', num2str(trip_id))};
                result_camera_parameters = fetchSensor(obj, table, fields, where_camera, '');
                result.rear_left_camera.parameters = result_camera_parameters;
                
                for i = 1:length(result.rear_left_camera.file_name)
                    result.rear_left_camera.file_name{i} = MapDatabase.createPathToImage(obj, result.rear_left_camera.file_name{i});
                end
                
                if obj.convert_GPS_to_ENU == 1 && (nargin == 3 && isfield(options.sensors,'base_station') && options.sensors.base_station == 1)
                    
                    ENU = g.WGSLLA2ENU(result.rear_left_camera.latitude, result.rear_left_camera.longitude, result.rear_left_camera.altitude);
                    station = obj.calculateStation(ENU');
                    
                    result.rear_left_camera.station = station;
                    result.rear_left_camera.X = ENU(1,:)';
                    result.rear_left_camera.Y = ENU(2,:)';
                    result.rear_left_camera.Z = ENU(3,:)';
                    
                end
                
            end
            
            % Rear right camera
            if nargin == 2 || (nargin == 3 && isfield(options.sensors,'rear_right_camera') && options.sensors.rear_right_camera == 1)
                
                if obj.verbose == 1
                    fprintf('Loading rear right camera data...\n')
                end
                
                sensor_id = 4;
                table = 'camera';
                fields = {'id', 'time', 'seconds_triggered', 'nanoseconds_triggered', 'file_name', 'latitude', 'longitude', 'altitude', 'roll', 'pitch', 'yaw'};
                where_camera = {where, cat(2,'sensor_id = ',num2str(sensor_id))};
                result_camera = fetchSensor(obj, table, fields, where_camera, orderby);
                result.rear_right_camera = result_camera;
                
                result.rear_right_camera.time_triggered = result.rear_right_camera.seconds_triggered + result.rear_right_camera.nanoseconds_triggered * 10 ^ (-9);
                result.rear_right_camera = rmfield(result.rear_right_camera, {'seconds_triggered', 'nanoseconds_triggered'});
                
                result_pose = fetchSensorPose(obj, sensor_pose_trip_id, sensor_id);
                result.rear_right_camera.pose = result_pose;
                
                table = 'camera_parameters';
                fields = {'id', 'focal_x', 'focal_y', 'center_x', 'center_y', 'skew', 'image_width', 'image_height', 'distortion_k1', 'distortion_k2', 'distortion_p1', 'distortion_p2', 'distortion_k3'};
                where_camera = {cat(2,'sensor_id = ', num2str(sensor_id)), cat(2,'trip_id = ', num2str(trip_id))};
                result_camera_parameters = fetchSensor(obj, table, fields, where_camera, '');
                result.rear_right_camera.parameters = result_camera_parameters;
                
                for i = 1:length(result.rear_right_camera.file_name)
                    result.rear_right_camera.file_name{i} = MapDatabase.createPathToImage(obj, result.rear_right_camera.file_name{i});
                end
                
                if obj.convert_GPS_to_ENU == 1 && (nargin == 3 && isfield(options.sensors,'base_station') && options.sensors.base_station == 1)
                    
                    ENU = g.WGSLLA2ENU(result.rear_right_camera.latitude, result.rear_right_camera.longitude, result.rear_right_camera.altitude);
                    station = obj.calculateStation(ENU');
                    
                    result.rear_right_camera.station = station;
                    result.rear_right_camera.X = ENU(1,:)';
                    result.rear_right_camera.Y = ENU(2,:)';
                    result.rear_right_camera.Z = ENU(3,:)';
                    
                end
                
            end
            
            % Front camera
            if nargin == 2 || (nargin == 3 && isfield(options.sensors,'front_camera') && options.sensors.front_camera == 1)
                
                if obj.verbose == 1
                    fprintf('Loading front camera data...\n')
                end
                
                sensor_id = 5;
                table = 'camera';
                fields = {'id', 'time', 'seconds_triggered', 'nanoseconds_triggered', 'file_name', 'latitude', 'longitude', 'altitude', 'roll', 'pitch', 'yaw'};
                where_camera = {where, cat(2,'sensor_id = ',num2str(sensor_id))};
                result_camera = fetchSensor(obj, table, fields, where_camera, orderby);
                result.front_camera = result_camera;
                
                result.front_camera.time_triggered = result.front_camera.seconds_triggered + result.front_camera.nanoseconds_triggered * 10 ^ (-9);
                result.front_camera = rmfield(result.front_camera, {'seconds_triggered', 'nanoseconds_triggered'});
                
                result_pose = fetchSensorPose(obj, sensor_pose_trip_id, sensor_id);
                result.front_camera.pose = result_pose;
                
                table = 'camera_parameters';
                fields = {'id', 'focal_x', 'focal_y', 'center_x', 'center_y', 'skew', 'image_width', 'image_height', 'distortion_k1', 'distortion_k2', 'distortion_p1', 'distortion_p2', 'distortion_k3'};
                where_camera = {cat(2,'sensor_id = ', num2str(sensor_id)), cat(2,'trip_id = ', num2str(trip_id))};
                result_camera_parameters = fetchSensor(obj, table, fields, where_camera, '');
                result.front_camera.parameters = result_camera_parameters;
                
                for i = 1:length(result.front_camera.file_name)
                    result.front_camera.file_name{i} = MapDatabase.createPathToImage(obj, result.front_camera.file_name{i});
                end
                
                if obj.convert_GPS_to_ENU == 1 && (nargin == 3 && isfield(options.sensors,'base_station') && options.sensors.base_station == 1)
                    
                    ENU = g.WGSLLA2ENU(result.front_camera.latitude, result.front_camera.longitude, result.front_camera.altitude);
                    station = obj.calculateStation(ENU');
                    
                    result.front_camera.station = station;
                    result.front_camera.X = ENU(1,:)';
                    result.front_camera.Y = ENU(2,:)';
                    result.front_camera.Z = ENU(3,:)';
                    
                end
                
            end
            
            % Zero all the timestamps.
            if obj.zero_time == 1
                
                if nargin == 2 || isempty(options.sensors)
                    sensors = {'gps', 'imu', 'laser', 'encoder_left', 'encoder_right', 'rear_left_camera', 'rear_right_camera', 'front_camera'};
                else
                    sensors = fieldnames(options.sensors);
                end
                
                start_times = [];
                for i = 1:length(sensors)
                    
                    if nargin == 2 || options.sensors.(sensors{i}) == 1
                        if ~strcmp(sensors{i}, 'base_station')
                            
                            if ~strcmp(sensors{i}, 'front_camera') && ~strcmp(sensors{i}, 'rear_left_camera') && ~strcmp(sensors{i}, 'rear_right_camera')
                                start_time = result.(sensors{i}).time(1);
                            else
                                start_time = result.(sensors{i}).time_triggered(1);
                            end
                            start_times = [start_times start_time];
                        end
                    end
                    
                end
                
                [obj.minimum_time,~] = min(start_times);
                
                for i = 1:length(sensors)
                    
                    if nargin == 2 || options.sensors.(sensors{i}) == 1
                        if ~strcmp(sensors{i}, 'base_station')
                            if strcmp(sensors{i}, 'front_camera') || strcmp(sensors{i}, 'rear_left_camera') || strcmp(sensors{i}, 'rear_right_camera')
                                result.(sensors{i}).time_triggered = result.(sensors{i}).time_triggered - obj.minimum_time;
                            end
                            result.(sensors{i}).time = result.(sensors{i}).time - obj.minimum_time;
                        end
                    end
                    
                end
                
            end
            
            % Separate the data by laps
            if obj.separate_by_lap == 1
                
                sensors = fieldnames(result);
                
                result.laps = {};
                for i = 1:length(sensors)
                    
                    if ~strcmp(sensors{i}, 'base_station')
                        
                        lap_ranges = obj.findLoops(result.(sensors{i}).X, result.(sensors{i}).Y, result.(sensors{i}).Z);
                        
                        data_fields = fieldnames(result.(sensors{i}));
                        
                        for j = 1:length(data_fields)
                            
                            if ~strcmp(data_fields{j}, 'pose') && ~strcmp(data_fields{j}, 'parameters')
                                
                                data = result.(sensors{i}).(data_fields{j});
                                
                                for k = 1:size(lap_ranges,1)
                                    
                                    result.laps{k}.(sensors{i}).(data_fields{j}) = data(lap_ranges(k,1):lap_ranges(k,2));
                                    
                                end
                                
                                result.(sensors{i}) = rmfield(result.(sensors{i}), data_fields{j});
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
            
            if obj.verbose == 1
                
                fprintf('Finish load the sensors data!\n');
                toc
                fprintf('=============================\n');
                
            end
            
        end
        
        function result = fetchSensorPose(obj, trip_id, sensor_id)
            
            table = 'sensor_poses';
            fields = {'x', 'y', 'z', 'roll', 'pitch', 'yaw'};
            where = {cat(2,'trip_id = ',num2str(trip_id)), cat(2,'sensor_id = ',num2str(sensor_id))};
            [result, column_names] = obj.db.select(table, fields, where);
            
            result = obj.db.convertFromTableToArray(column_names, result);
            
        end
        
        function result = fetchByDate(obj, date, options)
            
            where_date = MapDatabase.createDateWhereQuery(date);
            
            if nargin == 2
                result = fetchAll(obj, where_date);
            elseif nargin == 3
                result = fetchAll(obj, where_date, options);
            end
            
        end
        
        
        function result = fetchByLatitudeLongitudeWithRadius(obj, latitude, longitude, radius, options)
            
            where = MapDatabase.createLatitudeLongitudeRadiusWhereQuery(latitude, longitude, radius);
            
            if nargin == 4
                
                result = fetchAll(obj, where);
                
            elseif nargin == 5
                
                result = fetchAll(obj, where, options);
                
            end
            
        end
        
        function result = fetchByBagFileID(obj,bag_file_id,options)
            
            where = {};
            for i = 1:length(bag_file_id)
                where{i} = cat(2,'bag_file_id = ', num2str(bag_file_id(i)));
            end
            
            where = strjoin(where, ' OR ');
            where = {cat(2,'(', where, ')')};
            
            if nargin == 2
                result = fetchAll(obj, where);
            elseif nargin == 3
                result = fetchAll(obj, where, options);
            end
            
        end
        
        function result = fetchByBagFileName(obj,names,options)
            
            where = {};
            for i = 1:length(names)
                where{i} = cat(2,'name = ''', names{i}, '''');
            end
            
            where = {strjoin(where, ' OR ')};
            
            table = 'bag_files';
            fields = {'id'};
            orderby = 'datetime';
            [result, ~] = obj.db.select(table, fields, where, orderby);
            
            bag_file_ids = table2array(result.id);
            
            where = {};
            for i = 1:length(bag_file_ids)
                where{i} = cat(2,'bag_file_id = ', num2str(bag_file_ids(i)));
            end
            
            where = strjoin(where, ' OR ');
            where = {cat(2,'(',where,')')};
            
            if nargin == 2
                result = fetchAll(obj, where);
            elseif nargin == 3
                result = fetchAll(obj, where, options);
            end
            
        end
        
        % =============================== function fetchBagFileIDs ===========================
        % purpose:      query the all the data from bag_files table of DB
        % intput:       none,
        % output:       result_table = table format of query result
        %               result_struct = struct format of query result with array format element
        %               result = struct format of query result with table format element
        % =====================================================================================
        function [result_table,result_struct,result] = fetchBagFileIDs(obj)
            
            table = 'bag_files';
            fields = {'id','trips_id','name','datetime','datetime_end','parsed','date_added'};
            where = '';
            orderby = 'datetime';
            [result_table,result,column_names] = obj.db.select(table, fields, where, orderby);
            result_table = sortrows(result_table,{'datetime','id'}); % sort by date and id
            result_struct = obj.db.convertFromTableToStruct(column_names, result_table);
            
            if obj.verbose == 1
                fprintf('BAG FILES\n------\n')
                for i = 1:length(result_struct.id)
                    fprintf('ID: %i, Trips_id: %i, Name: %s, Datetime_start: %s, Datetime_end: %s, Parsed: %i, Date Added: %s',...
                        result_struct.id(i), result_struct.trips_id(i),result_struct.name{i}, result_struct.datetime{i},result_struct.datetime_end{i}, result_struct.parsed(i), result_struct.date_added{i})
                    fprintf('\n')
                end
                fprintf('\n')
            end
            
        end
        
        % =============================== function fetchTrips =================================
        % purpose:      query the all the data from trips table of DB
        % intput:       none,
        % output:       result_table = table format of query result
        %               result_struct = struct format of query result with array format element
        %               result = struct format of query result with table format element
        % =====================================================================================
        function [result_table,result_struct,result] = fetchTrips(obj)
            
            table = 'trips';
            fields = {'id', 'name', 'date','base_stations_id', 'description', 'driver','passengers', 'notes', 'date_added'};
            where = '';
            orderby = 'date';
            [result_table,result,column_names] = obj.db.select(table, fields, where, orderby);
            result_table = sortrows(result_table,{'date','date_added','id'}); % sort by date and id
            result_struct = obj.db.convertFromTableToStruct(column_names, result_table);
            
            if obj.verbose == 1
                fprintf('TRIPS\n------\n')
                for i = 1:length(result_struct.id)
                    fprintf('ID: %i, Name: %s, Driver: %s, Date: %s, Notes: %s', result_struct.id(i),result_struct.name{i}, result_struct.driver{i}, result_struct.date{i}, result_struct.notes{i})
                    fprintf('\n')
                end
                fprintf('\n')
            end
        end
        
        % =============================== function fetchByTripID===============================
        % purpose:      query data from BD given trip_id
        % intput:       trip_id = trip id, format: numeric array
        % output:       result_table = table format of query result
        %               result_struct = struct format of query result with array format element
        %               result = struct format of query result with table format element
        % processing flow: step1: find the bag_files_id given trips_id,
        %                  step2: query each table by bag_files_id and options
        % =====================================================================================
        function result = fetchByTripID(obj,trip_id,options)
            
            % step1: find the bag_files_id given trips_id,
            where = strjoin(cellstr(num2str(trip_id)), ', '); % char
            where = {cat(2,'trips_id in (', where, ')')}; % 1*1 cell
            
            %where = join(cellstr(num2str(trip_id)),' , ');
            
            table = 'bag_files';
            fields = {'id','datetime'};
            [result,~, ~] = obj.db.select(table,fields,where);
            result = sortrows(result,{'datetime','id'}); % sort by date and id
            
            bag_file_ids = result.id; % bag_files_id
            
            % step2:  query each table by bag_files_id and options
            %{
            OR is slow
            where = {};
            for i = 1:length(bag_file_ids.id)
                where{i} = cat(2,'bag_files_id = ', num2str(bag_file_ids(i)));
            end
            
            where = strjoin(where, ' OR ');
            where = {cat(2,'(', where, ')')};
            %}
            where = strjoin(cellstr(num2str(bag_file_ids)), ', '); % char
            where = {cat(2,'bag_files_id in (', where, ')')}; % 1*1 cell
            
            if nargin == 2
                result = fetchAll(obj, where);
            elseif nargin == 3
                result = fetchAll(obj, where, options);
            end
            
        end
        
        % =============================== function fetchByTimeRange===============================
        % purpose:      query data from BD by time range
        % intput:       start_time = start time, format:'yyyy-MM-dd HH:mm:ss'
        %               end_time = start time, format:'yyyy-MM-dd HH:mm:ss'
        % output:       result_table = table format of query result
        %               result_struct = struct format of query result with array format element
        %               result = struct format of query result with table format element
        % processing flow: step1: find the bag_files_id given trips_id,
        %                  step2: query each table by bag_files_id and options
        % =====================================================================================
        function result = fetchByTimeRange(obj,start_time, end_time,options)
            
            % step1: generate the where condition by start date and end
            % date
            where_date = MapDatabase.createDateWhereQuery(start_time, end_time);
            
            % step2:  query each table by bag_files_id and options
           
            if nargin == 3
                result = fetchAll(obj, where_date);
            elseif nargin == 4
                result = fetchAll(obj, where_date, options);
            end
            
        end
                
        %==========================
        function interpolated_result = fetchInterpolated(obj,table,fields,timestamp)
            
            if obj.zero_time == 1
                time_diff = cat(2,'ABS(', num2str(timestamp) ,' - (time - ', num2str(obj.minimum_time), ')) AS t_diff');
            else
                time_diff = cat(2,'ABS(', num2str(timestamp) ,' - time) AS t_diff');
            end
            fields = [time_diff, 'time', fields];
            orderby = 't_diff';
            limit = 2;
            
            [result, column_names] = obj.db.select(table, fields, '', orderby, limit);
            
            result = obj.db.convertFromTableToArray(column_names, result);
            
            [t, ind] = sort(result.time);
            
            % Remove the time and the difference in time
            result = rmfield(result, 'time');
            result = rmfield(result, 't_diff');
            
            fields = fields(3:end);
            interpolated_result = {};
            for i = 1:length(fieldnames(result))
                
                d = result.(fields{i});
                d = d(ind);
                
                if obj.zero_time == 1
                    t = t - obj.minimum_time;
                end
                
                interpolated_result.(fields{i}) = d(1) + (timestamp - t(1)) * diff(d) / diff(t);
                
            end
            
        end
        
        % Disconnect from the database.
        function disconnect(obj)
            
            close(obj.db.db_connection)
            
            fprintf('Disconnected\n')
            
        end
        
    end
    
    methods (Static)
        
        function query_string = createDateWhereQuery(date_start, date_end)
            
            if nargin == 1
                
                query_string = cat(2,'timestamp::DATE = ''', date_start, '''::DATE');
                
            elseif nargin == 2
                
                % query_string = cat(2,'timestamp::DATE >= ''', date_start, '''::DATE AND timestamp::DATE <= ''', date_end, '''::DATE');
                query_string = cat(2,'timestamp >= ''', date_start, ''' AND timestamp <= ''', date_end, '''');
            end
            
            % Join everything into one where query string
            query_string = strjoin(string(query_string),'');
            
        end
        
        function query_string = createLatitudeLongitudeRadiusWhereQuery(latitude, longitude, radius)
            
            query_string = cat(2,'ST_DWithin(ST_MakePoint(', num2str(longitude,16), ',', num2str(latitude,16), '), geography(geography),', num2str(radius), ')');
            
            % Join everything into one where query string
            query_string = strjoin(string(query_string),'');
            
        end
        
        function query_string = createTimestampSelectQuery()
            
            query_string = 'seconds + nanoseconds * 10^(-9) AS time';
            
        end
        
        function path = createPathToImage(obj,file_name)
            
            path = cat(2,obj.image_directory, file_name(1:2), '/', file_name(3:4), '/', file_name, '.jpg');
            
        end
        
        function station = calculateStation(X, Y, Z)
            
            if nargin == 1
                ENU = X;
                X = ENU(:,1);
                Y = ENU(:,2);
                Z = ENU(:,3);
            end
            
            station = sqrt(diff(X).^2 + diff(Y).^2 + diff(Z).^2);
            station = cumsum(station);
            station = [0; station];
            
        end
        
        function lap_ranges = findLoops(X, Y, Z)
            
            distance_between_points = sqrt((X - X(1)).^2 + (Y - Y(1)).^2 + (Z - Z(1)).^2);
            
            potential_loop_end_points = find(distance_between_points < 1);
            
            d_ind = find(diff(potential_loop_end_points) ~= 1) + 1;
            number_of_loops = length(d_ind) + 1;
            
            lap_ranges = zeros(number_of_loops,2);
            
            if number_of_loops == 1
                
                lap_ranges = [1 length(X)];
                
            else
                
                for i = 1:number_of_loops
                    
                    if i == 1
                        
                        if length(d_ind) == 1
                            
                            [~,j] = min(distance_between_points(potential_loop_end_points(d_ind(1):end)));
                            
                            test = potential_loop_end_points(d_ind(1):end);
                            
                            lap_ranges(i,:) = [1 test(j)];
                            %                             lap_ranges(i,:) = [1 potential_loop_end_points(j+d_ind(1))];
                            
                        else
                            
                            [~,j] = min(distance_between_points(potential_loop_end_points(d_ind(1):d_ind(2))));
                            
                            lap_ranges(i,:) = [1 potential_loop_end_points(j+d_ind(1))];
                            
                        end
                        
                    elseif i == number_of_loops
                        
                        loop_range = potential_loop_end_points(d_ind(i-1):end);
                        [~,index_of_loop_end] = min(distance_between_points(loop_range));
                        lap_ranges(i,:) = [loop_range(index_of_loop_end)+1 length(distance_between_points)];
                        
                    else
                        
                        if i == number_of_loops-1
                            
                            loop_range = potential_loop_end_points(d_ind(i):end-1);
                            
                        else
                            
                            loop_range = potential_loop_end_points(d_ind(i):d_ind(i+1)-1);
                            
                        end
                        
                        [~,index_of_loop_end] = min(distance_between_points(loop_range));
                        
                        lap_ranges(i,:) = [lap_ranges(i-1,2)+1 loop_range(index_of_loop_end)];
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
end