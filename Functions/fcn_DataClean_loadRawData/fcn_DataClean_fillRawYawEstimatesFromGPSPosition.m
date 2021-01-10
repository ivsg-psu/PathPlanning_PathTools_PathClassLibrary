function structure = fcn_DataClean_fillRawYawEstimatesFromGPSPosition(structure)

structure.Yaw_deg_from_position = ...
    atan2d(structure.yNorth_increments, structure.xEast_increments );

% Fix yaw angles to be positive numbers (0 to 360);
neg_yaw_indices = find(structure.Yaw_deg_from_position<0);
structure.Yaw_deg_from_position(neg_yaw_indices) = ...
    structure.Yaw_deg_from_position(neg_yaw_indices)+360;
return
