function Merged = fcn_mergeByTakingMedianOfSignals(d,sensors_list,field_names)


%% Set up the data
n_data = 0;
first_time = round(d.(sensors_list{1}).GPS_Time_deltaT/0.01);
for i_data = 1:length(sensors_list)
    sensor_name = sensors_list{i_data};
    for i_field = 1:length(field_names)
        field_name = field_names{i_field};
        n_data = n_data+1;
        data(:,n_data) = d.(sensor_name).(field_name); %#ok<AGROW>
        this_time = round(d.(sensor_name).GPS_Time_deltaT/0.01);
        if this_time~=first_time
            error('All dat must be equally time sampled');
        end
    end
end
%     [d.GPS_Hemisphere.Yaw_deg(:,1), ...
%     d.GPS_Hemisphere.Yaw_deg_from_position(:,1),...
%     d.GPS_Hemisphere.Yaw_deg_from_velocity(:,1),...
%     d.GPS_Novatel.Yaw_deg(:,1),...
%     d.GPS_Novatel.Yaw_deg_from_position(:,1),...
%     d.GPS_Novatel.Yaw_deg_from_velocity(:,1)];



%% Find outliers
[~, L,U,C] = isoutlier(data,2);

Merged.Center = C;
Merged.Upper = U;
Merged.Lower = L;
Merged.centiSeconds = first_time;
end