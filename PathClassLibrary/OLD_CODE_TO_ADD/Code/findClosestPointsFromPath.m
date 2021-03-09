%%%%%%%%%%%%%%%%%%%%%  Function findClosestPointsFromPath %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  
%      find the projection point on each traversals from reference path 
%
% Input Variables:
%      path          = reference path, format:struct with X,Y,Z,Station,and Yaw(optional)
%      data          = traversals data, format:struct with a traversal cell,
%                    and each attribute of the cell is a struct data with X,Y,Z,Station,and Yaw(optional)
%      flag_yaw      = set to 1 if you want to interplate the yaw too

% Returned Results:
%      [closestXs,closestYs,closestZs,closestYaws]  = the lane type of each lane, format:table
%
% Example:
% 
% 
% Processing Flow:
% 
% Restrictions/Notes:
%
% The following functions are called:
%      none
%
% Author:             Liming Gao
% Created Date:       2020-07-21
% Revisions:
%           2020-02-22: take yaw into account 
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [closestXs,closestYs,closestZs,closestYaws] = findClosestPointsFromPath(path,data,flag_yaw)
% path: the reference
if nargin < 3
    flag_yaw = 0;
end

% initializing variables
closestXs = zeros(length(path.X),length(data.traversal));
closestYs = zeros(length(path.X),length(data.traversal));
closestZs = zeros(length(path.X),length(data.traversal));
closestYaws = zeros(length(path.X),length(data.traversal));

% Loop through all traversals other than the 1st one. For each
% traversal, find the closest point on the first traversal, and
% use this to calculate station offset and index offset of the
% ith traversal relative to the first traversal

% calculate the bound we look for similar station, the bound should be no
% more than the difference between longest path and shortest path
safety_factor = 2;
last_station = zeros(length(data.traversal),1); % extract the max station of each path
for i_traversal = 1:length(data.traversal)
    last_station(i_traversal)  = data.traversal{i_traversal}.Station(end);
end

% need to check whether the path is longer than the longest trajectory or
% shorter than the shortest
longest_station = max(max(last_station), path.Station(end));
shortest_station = min(min(last_station), path.Station(end));
station_bound = safety_factor*(longest_station - shortest_station);

for index_reference_path = 1:length(path.X) % loop through all the point on the reference path
    
    % calculating x, y, and z position of one specific point on
    % each traversal given by index_to_match (a common point)
    X_ref = path.X(index_reference_path);
    Y_ref = path.Y(index_reference_path);
    Z_ref = path.Z(index_reference_path);
    Station_ref = path.Station(index_reference_path);
    station_lower = max(Station_ref - station_bound, 0);
    station_upper = min(Station_ref + station_bound, max(last_station)); % need confirm
    
    for i_traversal = 1:length(data.traversal) % loop through all traverasals (already excluded path)
        
        % allowable us to find the indexes within the station range
        index_in_station_range = find(data.traversal{i_traversal}.Station >= station_lower & data.traversal{i_traversal}.Station <= station_upper);
        
        % pull the X,Y,Z data for this range
        X_tra = data.traversal{i_traversal}.X(index_in_station_range);
        Y_tra = data.traversal{i_traversal}.Y(index_in_station_range);
        Z_tra = data.traversal{i_traversal}.Z(index_in_station_range);
        Yaw_tra = data.traversal{i_traversal}.Yaw(index_in_station_range);
        
        % finding the idx of the other traversal that is closet to the current reference point
        id = knnsearch([X_tra Y_tra Z_tra],[X_ref Y_ref Z_ref],'k',2);  %Find k-nearest neighbors using object
        
        id = sort(id); % sort in Ascending Order
        
        % find the projection point on traversal path
        nearest_point1 = [X_tra(id(1)); Y_tra(id(1));Z_tra(id(1))];
        nearest_point2 = [X_tra(id(2)); Y_tra(id(2));Z_tra(id(2))];
        query_point = [X_ref; Y_ref; Z_ref];
        
        % s = S(id(1));  % The station of the nearest point on the map
        % The vector from nearest point1 to nearest point2 on the traversal path
        np1_to_np2 = nearest_point2 - nearest_point1; 
        
        % The vector from nearest point1 to query point
        np1_to_qp = query_point - nearest_point1; 
        
        % the length vector
        np1_to_np2_Length = sum(np1_to_np2.^2).^0.5; 
        
        %pp means proejction point,and np1_to_pp_length is the projection
        %length of np1_to_qp to np1_to_np2, i.e. the distance from nearest
        %point1 to projection point (the value can be negtive)
        np1_to_pp_length = dot(np1_to_qp,np1_to_np2)./ np1_to_np2_Length;  
        
        Proj_point_ref_to_tra = nearest_point1 + np1_to_np2* diag(np1_to_pp_length/np1_to_np2_Length); %projection point  of second station on the reference path (first travel)
        
        closestXs(index_reference_path,i_traversal) = Proj_point_ref_to_tra(1);
        closestYs(index_reference_path,i_traversal) = Proj_point_ref_to_tra(2);
        closestZs(index_reference_path,i_traversal) = Proj_point_ref_to_tra(3);
        
        if flag_yaw == 1
            % find the closet yaw
            nearest_yaw1 = Yaw_tra(id(1));
            nearest_yaw2 = Yaw_tra(id(2));
            Proj_yaw_ref_to_tra = nearest_yaw1 + (nearest_yaw2-nearest_yaw1)*(np1_to_pp_length/np1_to_np2_Length);
            closestYaws(index_reference_path,i_traversal) = Proj_yaw_ref_to_tra;
        end
    end
    
end

end

