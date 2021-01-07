function Merged = fcn_mergeAllYawSources(d)

%% Set up the data
data = [d.GPS_Hemisphere.Yaw_deg(:,1), ...
    d.GPS_Hemisphere.Yaw_deg_from_position(:,1),...
    d.GPS_Hemisphere.Yaw_deg_from_velocity(:,1),...
    d.GPS_Novatel.Yaw_deg(:,1),...
    d.GPS_Novatel.Yaw_deg_from_position(:,1),...
    d.GPS_Novatel.Yaw_deg_from_velocity(:,1)];



%% Find outliers
[~, L,U,C] = isoutlier(data,2);

Merged.Yaw_Median = C;
Merged.Yaw_Upper = U;
Merged.Yaw_Lower = L;
end