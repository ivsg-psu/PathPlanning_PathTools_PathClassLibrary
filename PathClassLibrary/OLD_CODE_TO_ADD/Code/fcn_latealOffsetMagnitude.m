function [lateral_offset,distanceToNearestNeighbor] = fcn_latealOffsetMagnitude(targetENUArray,queryArray)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% calcualte the magnitude of offset by projection 
% issues: 1. change the for loop to matrix form 2. two dimensions to three
% dim


    m = size(queryArray, 1);                 % Length of the queryArray
    n = size(targetENUArray, 1);                        % Length of the path
    
    % Creates a matrix by replicating the array 'm' times
%     ENU_e = repelem(targetENUArray(:,1).', m, 1);  % East coords of points on a path
%     ENU_n = repelem(targetENUArray(:,2).', m, 1);  % North coords of points on a path
%     ENU_u = repelem(targetENUArray(:,3).', m, 1);  % Up coords of points on a path
      ENU_e = targetENUArray(:,1);
      ENU_n = targetENUArray(:,2);
    
    [Id_targetENUArray,distanceToNearestNeighbor] = knnsearch(targetENUArray,queryArray,'K',2);  %'euclidean' distance (default) |

     id = sort(Id_targetENUArray,2);
     lateral_offset=zeros(m,1);
     % for two dimensions temporarily
     for i = 1:m
         
        map1 = [ENU_e(id(i,1)); ENU_n(id(i,1))];
        map2 = [ENU_e(id(i,2)); ENU_n(id(i,2))];
        location = [queryArray(i,1);queryArray(i,2)];  
        
        ab = map2 - map1;   %vector 
        ab_squared = dot(ab,ab);% the square of length 
        ap = location - map1;
        t = dot(ap,ab) ./ ab_squared;  %%the normalized factor 
        point_on_line = map1 + ab * diag(t); %projection point  of second station on the reference path (first travel)
        lateral_offset(i) = sqrt((point_on_line(1)-location(1))^2+(point_on_line(2)-location(2))^2);  %

     end
%             
%             map1 = [X_map(id(1)); Y_map(id(1))];
%             map2 = [X_map(id(2)); Y_map(id(2))];
%             location = [x; y];  
%             
%             s = S(id(1));
%             ab = map2 - map1;
%             ab_squared = dot(ab,ab);%%%==========?????square root
%             ap = location - map1;
%             t = dot(ap,ab) ./ ab_squared;  %%the normalized factor 
%             point_on_line = map1 + ab * diag(t); %projection point  of second station on the reference path (first travel)
%             lateral_offset = sqrt((point_on_line(1)-location(1))^2+(point_on_line(2)-location(2))^2);  %
%             heading_map = psi_map(id(1));% atan2(map2(2)-map1(2),map2(1)-map1(1));
%     
    
    
   
   
    
end

