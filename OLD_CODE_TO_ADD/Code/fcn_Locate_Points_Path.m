% function [location] = fcn_Locate_Points_Path(queryArray, ENU)
% Find the position of points in queryArray wrt path represented by ENU - Left/Right
% queryArray - Set of points for which location needs to be determined
% ENU - ENU coordinates of the path
%
% Author: Srivenkata Satya Prasad Maddipatla
% Date: 07th Oct, 2019
% Edited: 16th Oct, 2019

function [location] = fcn_Locate_Points_Path(queryArray, ENU)
    
    m = size(queryArray, 1);                 % Length of the queryArray
    n = size(ENU, 1);                        % Length of the path
    
    % Creates a matrix by replicating the array 'm' times
    ENU_e = repelem(ENU(:,1).', m, 1);  % East coords of points on a path
    ENU_n = repelem(ENU(:,2).', m, 1);  % North coords of points on a path
    ENU_u = repelem(ENU(:,3).', m, 1);  % Up coords of points on a path
    
    % m x n matrix
    % Each row contains square of the distance from a point in the 
    % queryArray to each point on the path
    distances = ( (queryArray(:,1) - ENU_e).^2 +...
                  (queryArray(:,2) - ENU_n).^2 +...
                  (queryArray(:,3) - ENU_u).^2 );
    
    % Get index of the point on the path that is closest to the point
    % in the queryArray
    [~,I] = min(distances, [], 2);
    
    % Reduce the index by one if the index reaches the limit
    correct = find(I == n); % Find the indices where it reaches the limit
    I_copy = I;             % creates a copy of the indices
    I(correct) = I_copy(correct) - 1;
    
    % Vector in the path
    s_vector = ENU(I+1,:) - ENU(I,:);
    % Vector pointing towards the point, from the nearest point on the path
    p_vector = queryArray - ENU(I,:);
    
    % cross product between both the vectors
    c = cross(s_vector, p_vector);
    % Upward pointing component of the vectors
    up = c(:,3);
    
    % '0' - on the path; '-1' - Left of the path; '1' - Right of the path
    % Intuitively similar to number line
    location = sign(-up);
    
end