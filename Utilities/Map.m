classdef Map < handle
   
    % Class properties and variables
    properties
        config
        gps 
    end
    
    methods
        
        % Constructor for class.
        function obj = Map(config)
            obj.config = config;
            
            if ~isempty(obj.config) && ~isempty(obj.config.reference_latitude)
                % initialize the parameters for converting LLA to ENU
       
                obj.gps = GPS(obj.config.reference_latitude, obj.config.reference_longitude, obj.config.reference_altitude);
            else
                obj.gps = GPS();
            end
            
        end
        
        function data = processGPS(obj, data)
            %convert LLA to ENU 
            ENU = obj.gps.WGSLLA2ENU(data.latitude, data.longitude, data.altitude);
            station = obj.calculateStation(ENU');

            data.station = station;
            data.X = ENU(1,:)';
            data.Y = ENU(2,:)';
            data.Z = ENU(3,:)';
            
            data.yaw = atan2(diff(data.Y), diff(data.X));
            data.yaw = [data.yaw; data.yaw(end)];
            data.yaw = unwrap(data.yaw);
            
            % Calculate lateral and longitudinal velocity from INS
            data.U = data.east_velocity .* cos(data.yaw) + data.north_velocity .* sin(data.yaw);
            data.V = -data.east_velocity .* sin(data.yaw) + data.north_velocity .* cos(data.yaw);
                
            
        end
        
        function [id, station, lateral_offset, heading_map] = stationHere(obj,latitude_x,longitude_y,altitude_z,origin_x,origin_y,origin_z,S,X_map,Y_map,psi_map,type)

            if type == "lla"

                enu = obj.g.WGSLLA2ENU(latitude_x,longitude_y,altitude_z,origin_x,origin_y,origin_z);
                x = enu(1);
                y = enu(2);

            elseif type == "xyz"

                x = latitude_x;
                y = longitude_y;

            else

                error("fcn_station_here: Parameter 'type' must be lla or xyz\n")

            end

%             %nearest_station_ind
%             id = knnsearch([X_map Y_map],[x y]);
% 
%             s = S(id);
%             x_map = X_map(id);
%             y_map = Y_map(id);
%             heading_map = psi_map(id);
% 
%             delta_X = -x_map + x;
%             delta_Y = -y_map + y;
%             station = delta_X * cos(-heading_map) - delta_Y * sin(-heading_map) + s;
%             lateral_offset = delta_X * sin(-heading_map) + delta_Y * cos(-heading_map);

            id = knnsearch([X_map Y_map],[x y],'k',2);  %Find k-nearest neighbors using object
            
            id = sort(id);
            
            map1 = [X_map(id(1)); Y_map(id(1))];
            map2 = [X_map(id(2)); Y_map(id(2))];
            location = [x; y];  
            
            s = S(id(1));  % The station of the nearest point on the map
            ab = map2 - map1; % The vector of the 2 nearest points on the map to query point
            
            ab_squared = dot(ab,ab);%%%==========?????square root
            ap = location - map1; % Define a vector from nearest map point to query point
            
            t = dot(ap,ab) ./ ab_squared;  %%the normalized factor - this is the fraction of vector ap that is aligned in the direction of ab
            
            point_on_line = map1 + ab * diag(t); %projection point  of second station on the reference path (first travel)
            lateral_offset = sqrt((point_on_line(1)-location(1))^2+(point_on_line(2)-location(2))^2);  %
            heading_map = psi_map(id(1));% atan2(map2(2)-map1(2),map2(1)-map1(1));
            station = s + sqrt((point_on_line(1)-map1(1))^2+(point_on_line(2)-map1(2))^2);

        end
        
        % This function aligns data by station
        function data = alignDataByStation(obj, data, index_to_match, crop_to_common_station, interpolate_station, station_decimation)
            % initializing variables
            station_offsets = zeros(1,length(data.traversal)); 
            id_offset = zeros(1,length(data.traversal));
            station_offsets(1) = 0;
            id_offset(1) = 1;
            
            % Loop through all traversals other than the 1st one. For each
            % traversal, find the closest point on the first traversal, and
            % use this to calculate station offset and index offset of the
            % ith traversal relative to the first traversal
            
            for i = 2:length(data.traversal)
                % calculating x, y, and z position of one specific point on
                % each traversal given by index_to_match (a common point)
                X = data.traversal{i}.X(index_to_match);
                Y = data.traversal{i}.Y(index_to_match);
                Z = data.traversal{i}.Z(index_to_match);
                
                % finding the idx of the first traversal that is closet to
                % the current point
                [id, station, lateral_offset, heading_map] = obj.stationHere(X,Y,Z,0,0,0,data.traversal{1}.station,data.traversal{1}.X,data.traversal{1}.Y,data.traversal{1}.yaw,"xyz");
                station_offsets(i) = data.traversal{i}.station(index_to_match) - station; % distance offset of the current trajectory staiton from the first traj station 
                id_offset(i) = id(1); % index offset
                
                % correct the station by subtracting the offset
                data.traversal{i}.station = data.traversal{i}.station - station_offsets(i); 

            end
            
            if crop_to_common_station == 1
                
                % And we will chop off the beginning and ends of the data so they all have
                % the same length of data. First we need to find the maximum and minimum
                % station.
                start_stations = zeros(1,length(data.traversal));
                end_stations = zeros(1,length(data.traversal));
                for i = 1:length(data.traversal)
                    start_stations(i) = data.traversal{i}.station(1);
                    end_stations(i) = data.traversal{i}.station(end);
                end
                max_start_station = max(start_stations);
                min_end_station = min(end_stations);

                % Now crop the data by going through all the fields and
                % keeping only the indices that are within the bounds.
                for i = 1:length(data.traversal)

                   inds = find(data.traversal{i}.station >= max_start_station & data.traversal{i}.station <= min_end_station);

                    % Trim the data to the desired times
                    data.traversal{i}.time = data.traversal{i}.time(inds) - data.traversal{i}.time(inds(1));
                    data.traversal{i}.station = data.traversal{i}.station(inds);
%                   data.traversal{i}.latitude = data.traversal{i}.latitude(inds);
%                   data.traversal{i}.longitude = data.traversal{i}.longitude(inds);
%                   data.traversal{i}.altitude = data.traversal{i}.altitude(inds);
%                   data.traversal{i}.east_velocity = data.traversal{i}.east_velocity(inds);
%                   data.traversal{i}.north_velocity = data.traversal{i}.north_velocity(inds);
%                   data.traversal{i}.U = data.traversal{i}.U(inds);
%                   data.traversal{i}.V = data.traversal{i}.V(inds);
                    data.traversal{i}.X = data.traversal{i}.X(inds);
                    data.traversal{i}.Y = data.traversal{i}.Y(inds);
                    data.traversal{i}.Z = data.traversal{i}.Z(inds);
                    data.traversal{i}.yaw = data.traversal{i}.yaw(inds);
%                   data.traversal{i}.yaw_rate = data.traversal{i}.yaw_rate(inds);

                end
                
            end % Ends if crop to common station
            
            if interpolate_station == 1 % Do we interpolate?
                
                % Interpolate the data to be on the same station
                % decimation. First we need to find the minimum and maximum
                % station.
                start_stations = zeros(1,length(data.traversal));
                end_stations = zeros(1,length(data.traversal));
                for i = 1:length(data.traversal)
                    start_stations(i) = data.traversal{i}.station(1);
                    end_stations(i) = data.traversal{i}.station(end);
                end
                min_start_station = min(start_stations);
                max_end_station = max(end_stations);

                new_station = min_start_station:station_decimation:max_end_station;
                for i = 1:length(data.traversal)

                    % It's possible that there are non-unique station
                    % measurements, so you cannot interpolate the data.
                    % This grabs only the unique data indices
                    [~,unique_inds] = unique(data.traversal{i}.station);

%                   data.traversal{i}.latitude = interp1(data.traversal{i}.station(unique_inds), data.traversal{i}.latitude(unique_inds), new_station,'linear','extrap');
%                   data.traversal{i}.longitude = interp1(data.traversal{i}.station(unique_inds), data.traversal{i}.longitude(unique_inds), new_station,'linear','extrap');
%                   data.traversal{i}.altitude = interp1(data.traversal{i}.station(unique_inds), data.traversal{i}.altitude(unique_inds), new_station,'linear','extrap');
%                   data.traversal{i}.east_velocity = interp1(data.traversal{i}.station(unique_inds), data.traversal{i}.east_velocity(unique_inds), new_station,'linear','extrap');
%                   data.traversal{i}.north_velocity = interp1(data.traversal{i}.station(unique_inds), data.traversal{i}.north_velocity(unique_inds), new_station,'linear','extrap');
%                   data.traversal{i}.U = interp1(data.traversal{i}.station(unique_inds), data.traversal{i}.U(unique_inds), new_station,'linear','extrap');
%                   data.traversal{i}.V = interp1(data.traversal{i}.station(unique_inds), data.traversal{i}.V(unique_inds), new_station,'linear','extrap');
                    data.traversal{i}.X = interp1(data.traversal{i}.station(unique_inds), data.traversal{i}.X(unique_inds), new_station,'linear','extrap');
                    data.traversal{i}.Y = interp1(data.traversal{i}.station(unique_inds), data.traversal{i}.Y(unique_inds), new_station,'linear','extrap');
                    data.traversal{i}.Z = interp1(data.traversal{i}.station(unique_inds), data.traversal{i}.Z(unique_inds), new_station,'linear','extrap');
                    data.traversal{i}.yaw = interp1(data.traversal{i}.station(unique_inds), data.traversal{i}.yaw(unique_inds), new_station,'linear','extrap');
%                   data.traversal{i}.yaw_rate = interp1(data.traversal{i}.station(unique_inds), data.traversal{i}.yaw_rate(unique_inds), new_station,'linear','extrap');
                    data.traversal{i}.station = new_station - min_start_station;

                end
                
            end
            
        end
        
    end
    
    methods (Static)
        
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
        
        function [left_lane_edge, right_lane_edge] = createRoadBoundaries(lane_center, road_yaw, lane_width)
                                                  
            road_yaw = road_yaw - road_yaw(1);
            left_lane_edge = [lane_center(1,:) - (lane_width/2) * cos(road_yaw - pi/2);
                              lane_center(2,:) - (lane_width/2) * sin(road_yaw - pi/2)];

            right_lane_edge = [lane_center(1,:) + (lane_width/2) * cos(road_yaw - pi/2);
                               lane_center(2,:) + (lane_width/2) * sin(road_yaw - pi/2)];
            
        end
        
    end

end