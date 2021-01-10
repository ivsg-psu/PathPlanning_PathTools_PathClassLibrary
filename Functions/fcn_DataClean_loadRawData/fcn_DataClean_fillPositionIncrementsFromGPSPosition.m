function structure = fcn_DataClean_fillPositionIncrementsFromGPSPosition(structure)

% Find the position increments first
structure.xEast_increments = [0; diff(structure.xEast)];
structure.yNorth_increments = [0; diff(structure.yNorth)];
structure.zUp_increments = [0; diff(structure.zUp)];

% Fill in the correct standard deviations
% indices_GPS_active = structure.DGPS_is_active==1;
% structure.xEast_increments_Sigma   = structure.sigma_DGPS_inactive;
% structure.yNorth_increments_Sigma  = structure.sigma_DGPS_inactive;
% structure.zUp_increments_Sigma     = structure.sigma_DGPS_inactive;
% structure.xy_increments_Sigma      = structure.sigma_DGPS_inactive;
%
% structure.xEast_increments_Sigma(indices_GPS_active)  = structure.sigma_DGPS_active;
% structure.yNorth_increments_Sigma(indices_GPS_active) = structure.sigma_DGPS_active;
% structure.zUp_increments_Sigma(indices_GPS_active)    = structure.sigma_DGPS_active;
% structure.xy_increments_Sigma(indices_GPS_active)     = structure.sigma_DGPS_active;

% For increments, we use the sigmas from xEast, yNorth, zUp

structure.xEast_increments_Sigma   = structure.xEast_Sigma + [structure.xEast_Sigma(2:end,1); structure.xEast_Sigma(end,1)];
structure.yNorth_increments_Sigma  = structure.yNorth_Sigma + [structure.yNorth_Sigma(2:end,1); structure.yNorth_Sigma(end,1)];
structure.zUp_increments_Sigma     = structure.zUp_Sigma + [structure.zUp_Sigma(2:end,1); structure.zUp_Sigma(end,1)];

% Calculate xy increments
structure.xy_increments = ...
    (structure.xEast_increments.^2+structure.yNorth_increments.^2).^0.5;  %distance in meters
structure.xy_increments_Sigma = max(structure.xEast_increments_Sigma,structure.yNorth_increments_Sigma);

% For debugging
%std(diff(structure.xEast_increments(indices_GPS_active)))
%std(diff(structure.xEast_increments(indices_GPS_inactive)))
%figure; plot(diff(structure.xEast_increments))

return

