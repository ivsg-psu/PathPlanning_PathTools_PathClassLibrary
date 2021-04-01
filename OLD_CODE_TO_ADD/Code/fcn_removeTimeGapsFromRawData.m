function RawDataWithoutTimeGaps = fcn_removeTimeGapsFromRawData(rawData)
% This function enforces a time model in all the data that
% t(k) = t(k-1) + deltat;
%

flag_do_debug = 1;
flag_verbose_slip = 0;

%% Set the maximum sigma values one would expect for each field.
% This is what the sigma values become if missing
% fieldOrdering = [...
angleUncertainty = 200; % Units are degrees
velUncertainty = 100;  % Units are meters/second
positionUncertainty = 100; % Units are meters

maxSigmas.('Yaw_deg')=angleUncertainty;                 % Yaw variables
maxSigmas.('Yaw_deg_from_position')=angleUncertainty;
maxSigmas.('Yaw_deg_from_velocity')=angleUncertainty;
maxSigmas.('Yaw_deg_merged')=angleUncertainty;
maxSigmas.('ZGyro')=angleUncertainty*pi/180;                   % Yawrate (ZGyro) variables
maxSigmas.('velMagnitude')=velUncertainty;            % velMagnitude variables
maxSigmas.('XAccel')=velUncertainty;                  % XAccel variables
maxSigmas.('xEast_increments')=positionUncertainty;        % Position increment variables
maxSigmas.('yNorth_increments')=positionUncertainty;
maxSigmas.('zUp_increments')=positionUncertainty;
maxSigmas.('XYplot')=positionUncertainty;                  % XY plots
maxSigmas.('xEast')=positionUncertainty;                   % xEast and yNorth plots
maxSigmas.('yNorth')=positionUncertainty;
maxSigmas.('zUp')=positionUncertainty;
maxSigmas.('DGPS_is_active')=2;
maxSigmas.('velNorth')=velUncertainty;                % Remaining are not yet plotted - just kept here for now as  placeholders
maxSigmas.('velEast')=velUncertainty;
maxSigmas.('velUp')=velUncertainty;
maxSigmas.('Roll_deg')=angleUncertainty;
maxSigmas.('Pitch_deg')=angleUncertainty;
maxSigmas.('xy_increments') = positionUncertainty;
maxSigmas.('YAccel')=velUncertainty;
maxSigmas.('ZAccel')=velUncertainty;
maxSigmas.('XGyro')=angleUncertainty*pi/180; 
maxSigmas.('YGyro')=angleUncertainty*pi/180; 
maxSigmas.('VelocityR')=velUncertainty;




%% Start the code
if flag_do_debug
    % Grab function name
    st = dbstack;
    namestr = st.name;
    
    % Show what we are doing
    fprintf(1,'\nWithin function: %s\n',namestr);
    fprintf(1,'Starting iterations through rawData structure to determine if there are time gaps.\n');
end

fields_to_calculate_timegaps_for = [...
    {'ROS_Time'},...
    {'GPS_Time'}...
    ];



names = fieldnames(rawData); % Grab all the fields that are in rawData structure
for i_data = 1:length(names)
    % Grab the data subfield name
    data_name = names{i_data};
    d = rawData.(data_name);  % This is the input (current) data structure
    dout = d; % This is the output data structure
    
    if flag_do_debug
        fprintf(1,'\n Sensor %d of %d: ',i_data,length(names));
        fprintf(1,'Calculating timegap status for sensor: %s\n',data_name);
    end
    
    subfieldNames = fieldnames(d); % Grab all the subfields
    for i_subField = 1:length(subfieldNames)
        % Grab the name of the ith subfield
        subFieldName = subfieldNames{i_subField};
        
        % Check to see if this subField is in the list
        if any(strcmp(subFieldName,fields_to_calculate_timegaps_for))
            
            % This is the "meat" of this function, where we grab and
            % analyze the time data.
            
            % Grab the time vector
            t = d.(subFieldName);
            centiSeconds = d.centiSeconds;
            t_interval = t(end,1)-t(1,1);
            
            % Check if the number of samples make sense
            num_expected = round(t_interval/(centiSeconds*0.01)) + 1;
            
            % Label the result
            flag_length_error_detected = 0;
            if (length(t)>num_expected+2)|| (length(t)<num_expected-2)
                flag_length_error_detected = 1;
            end
            
            % Check to see if there are time jumps out of the ordinary
            diff_t = diff(t);
            mean_dt = mean(diff_t);
            std_dt = std(diff_t);
            max_dt = mean_dt+5*std_dt;
            min_dt = max(0.00001,mean_dt-5*std_dt);
            flag_jump_error_detected = 0;
            if any(diff_t>max_dt) || any(diff_t<min_dt)
                flag_jump_error_detected = 1;
            end
            
            % Update the user with progress
            if flag_do_debug
                fprintf(1,'\tProcessing subfield: %s \n',subFieldName);
                fprintf(1,'\t\tDuration in seconds: %f \n ',t_interval);
                fprintf(1,'\t\tSampled every %d ms\n ',centiSeconds);
                fprintf(1,'\t\tNum samples: %d\n',length(t(:,1)));
                fprintf(1,'\t\tNum expected: %d',num_expected);
                if 1==flag_length_error_detected
                    fprintf(1,' <--FAIL\n');
                else
                    fprintf(1,' <--Pass\n');
                end
                fprintf(1,'\t\tStart time: %f\n',t(1,1));
                fprintf(1,'\t\tEnd time: %f\n',t(end,1));
                fprintf(1,'\t\tEnd time expected: %f \n',t(1,1)+(num_expected-1)*0.01*centiSeconds);
                fprintf(1,'\t\tMean delta_t: %f \n',mean_dt);
                fprintf(1,'\t\tStandard deviation in delta_t: %f ',std_dt);
                if 1==flag_jump_error_detected
                    fprintf(1,' <--FAIL\n');
                else
                    fprintf(1,' <--Pass\n');
                end
                
            end
            
            
            if flag_length_error_detected || flag_jump_error_detected
                % Errors detected - find out why!
                if flag_do_debug
                    fprintf(1,'\t\tBeginning correction process: \n');
                end
                
                if flag_jump_error_detected
                    % Check to see if all the differences make sense, and where
                    % they do not. The way this is done is to, at each time
                    % step after t = 1, calculate the time step and compare
                    % this to the expected time step, plus/minus a tolerance.
                    
                    if flag_do_debug
                        fprintf(1,'\t\tLooking for incorrect changes in delta t: \n ');
                    end
                    flag_detecting_bad_now = 0;
                    flag_delta_errors_detected = 0;
                    bad_indices = [];
                    
                    for i_time=2:length(t)
                        tolerance = centiSeconds*0.5;
                        t_diff = (t(i_time)-t(i_time-1))*100;
                        if t_diff>(centiSeconds+tolerance) || t_diff<(centiSeconds-tolerance)
                            flag_delta_errors_detected = 1;
                            if flag_detecting_bad_now==0
                                i_bad_start = i_time;
                                flag_detecting_bad_now = 1;
                            end
                        else
                            if flag_detecting_bad_now==1
                                flag_detecting_bad_now = 0;
                                i_bad_end = i_time-1;
                                bad_indices = [bad_indices; (i_bad_start:i_bad_end)']; %#ok<AGROW>
                                tvec = t(i_bad_start-1:i_bad_end+1,1)';
                                diff_tvec = diff(tvec);
                                if 1==0 % For severe debugging
                                    fprintf(1,'\t\tBad interval detected from %d to %d\n',i_bad_start,i_bad_end);
                                    fprintf(1,'\t\t\tSample: ');
                                    fprintf(1,'%f\t',tvec);
                                    fprintf(1,'\n');
                                    fprintf(1,'\t\t\tDiff: ');
                                    fprintf(1,'%f\t',diff_tvec);
                                    fprintf(1,'\n\n');
                                end
                            end
                        end
                    end % Ends the for loop for the differences check
                    if flag_delta_errors_detected
                        if flag_do_debug
                            fprintf(1,'\t\tFixing time vector to repair incorrect jumps: \n ');
                        end
                        % Initialize a new time vector and fill it
                        t_new = t;
                        for i_time = 1:length(bad_indices)
                            bad_index = bad_indices(i_time);
                            t_new(bad_index,1) = t_new(bad_index-1,1)+centiSeconds*0.01;
                        end
                        dout.(subFieldName) = t_new;
                        
                        if flag_do_debug
                            
                            if 1==0 % For debugging only
                                figure(5684);
                                clf;
                                grid minor
                                hold on;
                                plot(diff(t),'r');
                                plot(diff(t_new),'b')
                            end
                            
                            fprintf(1,'\t\tThe following should match:\n');
                            fprintf(1,'\t\tOld end time: %f \n ',t(end));
                            fprintf(1,'\t\tNew end time: %f \n ',t_new(end));
                            
                            
                        end
                        
                        if flag_do_debug
                            fprintf(1,'\t\tFixing locations of incorrect jumps: \n ');
                        end
                        
                        
                        for i_subFields = 1:length(subfieldNames)
                            % Grab the name of the ith subfield
                            subFieldName2 = subfieldNames{i_subFields};
                            
                            % Check to see if this subField is in the list that
                            % does NOT include the time information and is
                            % not a sigma field
                            if ~any(strcmp(subFieldName2,fields_to_calculate_timegaps_for)) && ~contains(subFieldName2,'_Sigma')
                                % Grab the data vector
                                data = d.(subFieldName2);
                                
                                % If the data is same length of time, set
                                % values to NaN, then fix
                                if length(data(:,1))==length(t)
                                    data(bad_indices,1) = nan;
                                end

                                data = fillmissing(data,'linear');                                

                                % Push corrected data into the dout data
                                % structure
                                dout.(subFieldName2) = data;

                                % Now fix the sigma values
                                sigma_name = cat(2,subFieldName2,'_Sigma');
                                if isfield(d,sigma_name)
                                    sigma_data = d.(sigma_name);
                                    
                                    % Check to see if this is a scalar. If
                                    % it is, need to convert it to a vector
                                    if length(sigma_data(:,1))==1
                                        sigma_data = sigma_data*ones(length(data(:,1)),1);
                                    elseif length(sigma_data(:,1)) ~= length(data(:,1))
                                        error('Data and corresponding sigma values are unequal in length. Unsure how to continue.');
                                    end
                                    
                                    if isfield(maxSigmas,subFieldName2)
                                        sigma_data(bad_indices,1) = maxSigmas.(subFieldName2);
                                    else
                                        error('Sigma field detected, but corresponding field not found in maxSigma list: %s',subFieldName2);
                                    end

                                    dout.(sigma_name) = sigma_data;
                                end
                                
                            end
                        end
                        t = t_new;
                    end % Ends flag for if on flag_delta_errors_detected
                    
                    if flag_do_debug
                        if ~flag_delta_errors_detected
                            fprintf(1,'\t\tNo delta errors detected.\n');
                        else
                            fprintf(1,'\t\tDone fixing jumps in delta t.\n ');
                        end
                    end
                end % end if for flag_jump_error_detected
                
                if 1==0 && flag_length_error_detected
                    % THE FOLLOWING DOES NOT WORK YET, AS OF 2019_11_22
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Check to see if the land points make sense, and where
                    % they do not.
                    if flag_do_debug
                        fprintf(1,'\t\tLooking for slipping between t and reference: \n ');
                    end
                    
                    slip = 0;
                    error_time = 0*t;
                    slipping = 0*t;
                    predicted_t_land_points = t(1,1):centiSeconds*0.01:(t(1,1)+centiSeconds*0.01*(num_expected-1));
                    flag_detecting_bad_now = 0;
                    tolerance = centiSeconds*0.5*0.01;
                    for i_time=1:length(t)
                        try
                            tpred = predicted_t_land_points(i_time+slip);
                        catch
                            disp('Debug here');
                        end
                        time_error = t(i_time) - tpred;
                        error_time(i_time) = time_error;
                        slipping(i_time) = slip;
                        if time_error>tolerance || time_error<(-tolerance)
                            if flag_detecting_bad_now==0
                                i_bad_start = i_time;
                                flag_detecting_bad_now = 1;
                            end
                            if time_error>0 % Time is ahead of prediction, increase slip
                                slip = slip + 1;
                            else % prediction is ahead of time, decrease slip
                                slip = slip - 1;
                            end
                        else
                            if 1==flag_detecting_bad_now
                                flag_detecting_bad_now = 0;
                                i_bad_end = i_time-1;
                                tvec = t(i_bad_start-1:i_bad_end+1,1)';
                                tvecpred = predicted_t_land_points(i_bad_start-1:i_bad_end+1);
                                if (flag_do_debug && flag_verbose_slip)
                                    fprintf(1,'\t\tBad interval detected from %d to %d\n',i_bad_start,i_bad_end);
                                    fprintf(1,'\t\t\tSlipped by %d\n',slip);
                                    fprintf(1,'\t\t\tSample: ');
                                    fprintf(1,'%f\t',tvec);
                                    fprintf(1,'\n');
                                    fprintf(1,'\t\t\tPred: ');
                                    fprintf(1,'%f\t',tvecpred);
                                    fprintf(1,'\n\n');
                                end
                            end
                        end
                    end % Ends the for loop for time slippage
                    if flag_do_debug
                        figure(46446473);
                        clf;
                        hold on;
                        grid minor;
                        num_common_points = min(length(t),length(predicted_t_land_points));
                        plot(t(1:num_common_points),(t(1:num_common_points)-predicted_t_land_points(1:num_common_points)')/tolerance,'g');
                        plot(t,error_time/tolerance,'r');
                        plot(t,slipping,'b');
                        legend('Error between t and predicted_t, no correction','Error between t and predicted_t, after correction','Slippage');
                    end
                end % ends if statement for flag_length_error_detected
                
            end % Ends if for error being detected
            
            
        else
            % If enter here, then we have a field that does NOT need a
            % time check calculation. So we can do nothing...
            
        end % Ends the if statement to check if subfield is on list
    end % Ends for loop through the subfields
    
    RawDataWithoutTimeGaps.(data_name) = dout; % Save results to main structure
    
    
    
end  % Ends for loop through all sensor names in rawData

if flag_do_debug
    % Show what we are doing
    fprintf(1,'\nFinished processing function: %s\n',namestr);
end

end % Ends the function





%% Subfunctions start here
function real_sigma = fcn_calcSigmaNoOutliers(data)
differences = diff(data);
deviations = differences - mean(differences);
outlier_sigma = std(deviations);
% Reject outliers
deviations_with_no_outliers = deviations(abs(deviations)<(3*outlier_sigma));
real_sigma = std(deviations_with_no_outliers);

%real_sigma_vector = real_sigma*
end
