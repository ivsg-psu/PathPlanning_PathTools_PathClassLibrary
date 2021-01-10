function structure = fcn_DataClean_fillRawYawEstimatesFromGPSVelocity(structure)

structure.Yaw_deg_from_velocity = ...
    atan2d(structure.velNorth,structure.velEast);

% Fix yaw angles to be positive numbers (0 to 360);
neg_yaw_indices = find(structure.Yaw_deg_from_velocity<0);
structure.Yaw_deg_from_velocity(neg_yaw_indices) = ...
    structure.Yaw_deg_from_velocity(neg_yaw_indices)+360;

return