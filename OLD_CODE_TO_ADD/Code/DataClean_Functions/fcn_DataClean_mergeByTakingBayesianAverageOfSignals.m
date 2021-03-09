function Merged = fcn_DataClean_mergeByTakingBayesianAverageOfSignals(d,sensors_list,field_names,truth_sensor,truth_field)
% This is a function that merges data by using Baysian averaging (by
% sigmas)
% 
% INPUTS:
%
%     d - time filetered data
%     sensors_list - List of sensors under consideration
%     field_names - Different ways in which data is obtained
%     truth_sensor - Reference sensor/Novatel
%     truth_field - Reference data/Output
% 
% Field name of std-dev of 'field_name' should be 'field_name_Sigma'
%
% Authors: Dr. Sean Brennan
% Contact: sbrennan@psu.edu
%
% Edited by: Satya Prasad Maddipatla, szm888@psu.edu
%
% Revision history:
%
%     2019_10_21 - First write based on Dr.Brennan's version, Replaced arithmetic average with bayesian average
%     2019_11_25 - Cleaned up bugs for when empty fields are passed into this
%     function. Made process more verbose.
%     2019_11_26 - Added more verbose information for when interpolation
%     throws a warning. Added cprintf functionality.
%     2020_11_10 - changed function names in prep for DataClean class (SNB)

flag_do_debug = 1;
flag_show_plots = 0;

if flag_do_debug
    % Grab function name
    st = dbstack;
    namestr = st.name;

    % Show what we are doing
    fprintf(1,'\nWithin function: %s\n',namestr);
    fprintf(1,'\tMerging information from sensors: \n');
    for i_sensor = 1:length(sensors_list)
        fprintf(1,'\t\t Sensor %d: %s',i_sensor,sensors_list{i_sensor});
        if strcmp(sensors_list{i_sensor},truth_sensor)
            fprintf(1,'<-- Assumed to be truth sensor for time base and correlation\n');
        else
            fprintf(1,'\n');
        end
    end
    fprintf(1,'\tMerging information from the following fields in each sensor: \n');
    for i_field = 1:length(field_names)
        fprintf(1,'\t\t Field %d: %s',i_field,field_names{i_field});    
        if strcmp(field_names{i_field},truth_field)
            fprintf(1,'<-- Assumed to be truth field for time base and correlation\n');
        else
            fprintf(1,'\n');
        end
    end   
end

%% Error checking
% Check to see if the number of arguments is correct
if nargin ~= 5
    error('Invalid number of arguments: Expecting five arguments.');
end



%% Set up the time and loop data, and do basic error checking
n_data       = 0;   % The number of data that has been processed

% Make sure truth sensor exists
if ~isfield(d,truth_sensor)
    error('Looking for the truth sensor %s, but cannot find it.',truth_sensor);
end

% Make sure centiSeconds field exists, and grab it from truth sensor
if ~isfield(d.(truth_sensor),'centiSeconds')
    error('Looking for the centiSeconds field in sensor %s, but cannot find it.',truth_sensor);
else
    truth_centiSeconds = d.(truth_sensor).centiSeconds;    % The centiseconds for the truth data
end

% Make sure Clocks field exists, and grab master time from truth sensor's
% value
if ~isfield(d,'Clocks')
    error('Looking for the Clocks sensor field, but cannot find it.');
end
if ~isfield(d.Clocks,'targetTimeVector_GPS')
    error('Looking for the targetTimeVector_GPS field within the Clocks sensor field, but cannot find it.');
end
if length(d.Clocks.targetTimeVector_GPS)<truth_centiSeconds
    error('The number of target time vectors in GPS time is not sufficient for the number of centiSeconds in the truth sensor.');
end

truth_t_vec  = d.Clocks.targetTimeVector_GPS{truth_centiSeconds};   % Time vector for the truth data

%% Start the loop that gathers all the data and sigma values into arrays for each
% Loop through all the sensors listed, and gather their data. But check to
% see if that data is equally time sampled. And if it isn't (for example, if
% one sensor is sampled at 100 Hz and another at 20 Hz), then throw a
% warning and resample it to the truth_data time base. Each sensor/field
% combination contributes one column to the resulting data and sigma
% matrices.

for i_data = 1:length(sensors_list)
    sensor_name = sensors_list{i_data};

    % Grab time values
    if ~isfield(d.(sensor_name),'centiSeconds')
        error('Looking for the centiSeconds field in sensor %s, but cannot find it.',sensor_name);
    else
        this_centiSeconds   = d.(sensor_name).centiSeconds;
    end
    temp_t  = d.Clocks.targetTimeVector_GPS{this_centiSeconds};   % Time vector for the current data

    
    for i_field = 1:length(field_names)
        field_name = field_names{i_field};
        sigma_name = strcat(field_name,'_Sigma');
        
        % Test to see if field exists
        if ~isfield(d.(sensor_name),field_name) || ~isfield(d.(sensor_name),sigma_name) 
            fprintf(1,'\t\t WARNING: Searching within sensor: %s for fields: %s, and %s but one of these does not exist. Skipping...',sensor_name,field_name,sigma_name);            
        else % The field exists, and can proceed
            temp_data  = d.(sensor_name).(field_name);
            temp_sigma = d.(sensor_name).(sigma_name);           

            flag_data_good = 1;
            % Confirm that the time and data vectors are the same length
            if length(temp_t(:,1))~=length(temp_data(:,1))
                flag_data_good = 0; %#ok<NASGU>
                fprintf(1,'Time and data lengths mismatch in field: %s for sensor %s\n',field_name,sensor_name);
                fprintf(1,'Time length: %d\n',length(temp_t(:,1)));
                fprintf(1,'Data length: %d\n',length(temp_data(:,1)));

                error('Time and data lengths mismatch in field: %s for sensor %s',field_name,sensor_name);
            end
            
            % Confirm there are non NaNs within data fields
            if any(isnan(temp_data))
                warning('NaN values found in field: %s for sensor %s',field_name,sensor_name);
                flag_data_good = 0;
            end
            
            % Confirm there are non NaNs within sigma fields
            if any(isnan(temp_sigma))
                warning('NaN values found in data sigma field: %s for sensor %s',strcat(field_name,'_Sigma'),sensor_name);
                flag_data_good = 0;
            end
            
            % Confirm there are no negative or zero values within the sigma
            % field
            if any(temp_sigma<=0)
                warning('Zero or negative values found in data sigma field: %s for sensor %s',strcat(field_name,'_Sigma'),sensor_name);
                flag_data_good = 0;
            end

            % Cannot enter here unless all data is good
            if 1==flag_data_good
                if this_centiSeconds ~= truth_centiSeconds
                    cprintf('Errors','\tWARNING: Data is not equally time sampled. Interpolation will be used.\n');
                    fprintf(1,'\tThe truth sensor: %s, has a time sampling of %d centiSeconds.\n',truth_sensor,truth_centiSeconds);
                    fprintf(1,'\tThis sensor: %s, has a time sampling of %d centiSeconds.\n',sensor_name,this_centiSeconds);
                    
                    temp_data  = interp1(temp_t,temp_data,truth_t_vec);
                    % sigma doesn't vary much
                    if 1==length(temp_sigma(:,1))
                        temp_sigma = temp_sigma * ones(size(temp_t,1),1);
                    end
                    temp_sigma = interp1(temp_t,temp_sigma,truth_t_vec);
                    
                end
                
                % If the sigma value is a scalar, convert it to a vector as
                % all sigma values need to be same length as data vector
                if 1==length(temp_sigma)
                    temp_sigma = temp_sigma * ones(size(temp_t,1),1);
                end
                
                n_data          = n_data+1;
                data_columns(:,n_data)  = temp_data;    %#ok<AGROW>
                sigma_columns(:,n_data) = temp_sigma;   %#ok<AGROW>
                
            else
                warning('NaN detected in %s.%s or its correspoding std-dev, but it works unless every field contains NaN',sensor_name,field_name);
                
            end % ends if statement checking if there is any NaN values
        end % ends if statement checking if field exists
        
    end % ends for loop through field names
    
end % ends for loop through sensors

%% Perform cross-correllation
% Some of the data is not perfectly time aligned (the ADIS for example),
% and so we perform a cross-correllation here to re-align the data
% correctly in time.

truth_data = d.(truth_sensor).(truth_field);

% Loop through each data vector and check cross-correllation with the truth
% data.
for i_data = 1:length(data_columns(1,:))
    test_data  = data_columns(:,i_data);
    test_sigma = sigma_columns(:,i_data);
    
    % Check cross-correllation over 200 milliseconds (sample rate is 10 ms);
    [cross_correlation,lags] = xcorr(truth_data,test_data,40);
    
    if 1 == flag_show_plots
        figure(363563);
        plot(lags,cross_correlation);
    end
    
    [~,max_correllation_index] = max(cross_correlation);
    index_offset               = lags(max_correllation_index);
    
    % Shift the data to match
    N = length(test_data(:,1));
    if index_offset < 0     % Shift foward in time
        index_offset   = -index_offset;
        dest_indices   = 1:(N-index_offset);
        source_indices = (index_offset+1):N;
    else        % Shift backward in time
        source_indices = 1:(N-index_offset);
        dest_indices   = (index_offset+1):N;
    end
    test_data(dest_indices) = test_data(source_indices);
    % shift std-dev in sync with the data
    test_sigma(dest_indices) = test_sigma(source_indices);
    
    if 1 == flag_show_plots
        % Check the result
        [cross_correlation,lags] = xcorr(truth_data,test_data,40);
        figure(363563);
        plot(lags,cross_correlation);
        [~,max_correllation_index] = max(cross_correlation);
        index_offset = lags(max_correllation_index); %#ok<NASGU>
    end
    
    data_columns(:,i_data) =  test_data;
    sigma_columns(:,i_data) = test_sigma;
    
end

%% Find averages, upper, and lower
% Now that all the data is aligned, we can find the bayesian average, min, and max.

% Check that all the sigma values are positive
if max(any(sigma_columns<0))
    fprintf('\t\t WARNING: zero sigma values detected in following columns:\n');
    fprintf('\t\t Column:   ');
    fprintf('%d \t',1:length(sigma_columns(1,:)));
    fprintf('\n');
    fprintf('\t\t Detected? ');
    fprintf('%d \t',any(sigma_columns<0));
    fprintf('\n');
end

if max(any(sigma_columns==0))
    fprintf('\t\t WARNING: zero sigma values detected in following columns:\n');
    fprintf('\t\t Column:   ');
    fprintf('%d \t',1:length(sigma_columns(1,:)));
    fprintf('\n');
    fprintf('\t\t Detected? ');
    fprintf('%d \t',any(sigma_columns==0));
    fprintf('\n');
end



[C, C_sigma] = fcn_DataClean_bayesianAverageMatrixForm(data_columns,sigma_columns);
U = max(data_columns,[],2);
L = min(data_columns,[],2);

Merged.Center = C;
Merged.Upper  = U;
Merged.Lower  = L;
Merged.Sigma  = C_sigma;

Merged.centiSeconds = truth_centiSeconds;

if flag_do_debug
    fprintf(1,'\nExiting function: %s\n',namestr);
end
return