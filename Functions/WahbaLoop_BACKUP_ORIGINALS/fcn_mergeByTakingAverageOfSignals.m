function Merged = fcn_mergeByTakingAverageOfSignals(d,sensors_list,field_names,truth_sensor,truth_field)

flag_do_debug = 0;

%% Set up the data
n_data = 0;  % The number of data that has been processed
first_deltat = round(d.(truth_sensor).ROS_Time_deltaT/0.01);  % The centiseconds for the truth data
first_t_vec = d.(truth_sensor).ROS_Time;  % Time vector for the truth data

% Loop through all the sensors listed, and gather their data. But check to
% see if that data is equally time sampled. If it isn't (for example, if
% one sensor is sampled at 100 Hz and another at 20 Hz), then throw a
% warning and resample it to the truth_data time base.
for i_data = 1:length(sensors_list)
    sensor_name = sensors_list{i_data};
    for i_field = 1:length(field_names)
        field_name = field_names{i_field};
        temp_data = d.(sensor_name).(field_name);
        temp_t    = d.(sensor_name).ROS_Time;
        
        this_time = round(d.(sensor_name).ROS_Time_deltaT/0.01);
        if this_time~=first_deltat
            warning('Data is not equally time sampled. Interpolation will be used...');
            temp_data = interp1(temp_t,temp_data,first_t_vec);
        end
        n_data = n_data+1;
        data(:,n_data) = temp_data; %#ok<AGROW>
    
    end
end

%% Perform cross-correllation
% Some of the data is not perfectly time aligned (the ADIS for example),
% and so we perform a cross-correllation here to re-align the data
% correctly in time. 

truth_data = d.(truth_sensor).(truth_field);

% Loop through each data vector and check cross-correllation with the truth
% data.
for i_data = 1:length(data(1,:))
    test_data = data(:,i_data);
    
    % Check cross-correllation over 200 milliseconds (sample rate is 10 ms);
    [cross_correlation,lags] = xcorr(truth_data,test_data,40);
    
    if 1==flag_do_debug
        figure(363563);
        plot(lags,cross_correlation);
    end
    [~,max_correllation_index] = max(cross_correlation);
    index_offset = lags(max_correllation_index);
    
    % Shift the data to match
    N = length(test_data(:,1));
    if index_offset<0  % Shift foward in time
        index_offset   = -index_offset;
        dest_indices   = 1:(N-index_offset);
        source_indices = (index_offset+1):N;
    else  % Shift backward in time
        source_indices = 1:(N-index_offset);
        dest_indices   = (index_offset+1):N;
    end
    test_data(dest_indices) = test_data(source_indices);
    
    if 1==flag_do_debug
        % Check the result
        [cross_correlation,lags] = xcorr(truth_data,test_data,40);
        
        figure(363563);
        plot(lags,cross_correlation);
        [~,max_correllation_index] = max(cross_correlation);
        index_offset = lags(max_correllation_index); %#ok<NASGU>
    end
    
    data(:,i_data) =  test_data;
    
end


%% Find averages, upper, and lower
% Now that all the data is aligned, we can find the mean, min, and max.
% Note: need to replace the mean calculation here with the bayesian
% averaging.

C = mean(data,2);
U = max(data,[],2);
L = min(data,[],2);

Merged.Center = C;
Merged.Upper = U;
Merged.Lower = L;

Merged.centiSeconds = first_deltat;
end

