function medianFilteredRawDataWithSigmas = fcn_DataClean_medianFilterFromRawAndSigmaData(rawDataWithSigmas)

flag_do_debug = 1;

if flag_do_debug
    % Grab function name
    st = dbstack;
    namestr = st.name;

    % Show what we are doing
    fprintf(1,'\nWithin function: %s\n',namestr);
    fprintf(1,'Starting iterations through rawDataWithSigmas structure to calculate sigma values.\n');    
end

fields_to_calculate_median_filters_for = [...
    {'velNorth'},...
    {'velEast'},...
    {'velUp'},...
    {'velMagnitude'},...
    {'Roll_deg'},...
    {'Pitch_deg'},...
    {'Yaw_deg'}...  % Confirmed
    {'Yaw_deg_from_position'}... % Confirmed
    {'Yaw_deg_from_velocity'}... % Confirmed
    {'xy_increments'}... % Confirmed
    {'XAccel'}...
    {'YAccel'}...
    {'ZAccel'}...
    {'XGyro'}...
    {'YGyro'}...
    {'ZGyro'}...
    {'xEast_increments'}...
    {'yNorth_increments'}...
    {'xEast'}...
    {'yNorth'}...
    ];

fields_to_unwrap_angles = [...
    {'Yaw_deg'}...  % Confirmed
    {'placeholder'}...
    ];

fields_to_check_diff_not_zero = [...
    {'xEast_increments'}... 
    {'yNorth_increments'}...
    ];

names = fieldnames(rawDataWithSigmas); % Grab all the fields that are in rawData structure
for i_data = 1:length(names)
    % Grab the data subfield name
    data_name = names{i_data};
    d = eval(cat(2,'rawDataWithSigmas.',data_name));
    
    if flag_do_debug
        fprintf(1,'\n Sensor %d of %d: ',i_data,length(names));
        fprintf(1,'Performing median filter calculations for sensor: %s\n',data_name);
    end
    
    subfieldNames = fieldnames(d); % Grab all the subfields
    clear dout; % Initialize this structure
    for i_subField = 1:length(subfieldNames)
        % Grab the name of the ith subfield
        subFieldName = subfieldNames{i_subField};

        if flag_do_debug
            fprintf(1,'\tProcessing subfield: %s ',subFieldName);
        end
        
        % Copy over the field itself first
        dout.(subFieldName) = d.(subFieldName);
        
        % Check to see if this subField is in the list
        if any(strcmp(subFieldName,fields_to_calculate_median_filters_for))
            % Grab data
            data = d.(subFieldName);
                       
            % Grab sigma of data (MUST exist)
            subFieldNameSigma = cat(2,subFieldName,'_Sigma');
            if ~isfield(d,subFieldNameSigma)
                % The sigma field does not exist - this would break the code
                fprintf(1,'\n\n SUBFIELD DOES NOT EXIST: %s ',subFieldName);
                error('Unable to continue');
            end
            
            % Grab sigma
            sigmas = d.(subFieldNameSigma);
                        
            % Check if the data needs to be unwrapped first
            if any(strcmp(subFieldName,fields_to_unwrap_angles))
                data = fcn_DataClean_unwrapAngles(data);
            end  
            
            % Check if the data has absolutely zero change
            if any(strcmp(subFieldName,fields_to_check_diff_not_zero))
                good_indices = find(abs(data)>0.000000001); % Set a very low tolerance...
                data = fcn_DataClean_replaceBadIndicesWithNearestGood(data,good_indices);
                
            end
            
            % Perform median filter (Edited this section on 2019_11_27 to
            % include Sigma calculations)
            [median_cleaned_data, sigma_cleaned_data] = fcn_DataClean_medianFilterData(data,2*sigmas);
            dout.(subFieldName) = median_cleaned_data;
            dout.(subFieldNameSigma) = sigma_cleaned_data;
            sigmas = sigma_cleaned_data;
            % End of edits on 2019_11_27
            
            % Check if Sigma field is a vector. If so, NOT overwrite it.
            % Just pass it through.
            if length(sigmas)>1
                % This is a vector of sigma values - should keep it as this
                % means it was entered directly earlier
                
                % Check if sigma vector has nan's within
                if any(isnan(sigmas))
                    sigma_via_next = fillmissing(sigmas,'next');
                    sigma_via_previous = fillmissing(sigmas,'previous');
                    sigmas_new = max([sigma_via_next,sigma_via_previous],2);
                    if flag_do_debug
                        fprintf(1,' <-- fixed NaNs in this sigma, and did median filter\n');
                    end
                else   % There are non NaN's in the sigma vector              
                    if flag_do_debug
                        fprintf(1,' <-- skipped this sigma, but did median filter\n');
                    end
                    sigmas_new = sigmas;
                end
                dout.(subFieldNameSigma) = sigmas_new;
            else
                real_sigma = fcn_DataClean_calcSigmaNoOutliers(median_cleaned_data);
                dout.(subFieldNameSigma) = real_sigma;
                
                if flag_do_debug
                    % fprintf(1,' <-- calculated a sigma, has length: %d\n',length(real_sigma(:,1)));
                    fprintf(1,' <-- calculated a median filter and new sigma\n'); 
                end
            end % Ends the if on whether the subfield already exists
        else 
            % If enter here, then we have a field that does NOT need a
            % sigma calculation. So we can do nothing...
            
            if flag_do_debug
                fprintf(1,'<-- skipped this. Not on list.\n');
            end
        end % Ends the if statement to check if subfield is on list
    end % Ends for loop through the subfields

    medianFilteredRawDataWithSigmas.(data_name) = dout; % Save results to main structure
    
end  % Ends for loop through all sensor names in rawData


if flag_do_debug
    % Show what we are doing
    fprintf(1,'\nFinished processing function: %s\n',namestr);
end

return % Ends the function




%% Subfunctions start here
%% Sigma calculation function
function real_sigma = fcn_DataClean_calcSigmaNoOutliers(data)
differences = diff(data);
deviations = differences - mean(differences);
outlier_sigma = std(deviations);
% Reject outliers
deviations_with_no_outliers = deviations(abs(deviations)<(3*outlier_sigma));
real_sigma = std(deviations_with_no_outliers);
return

%% Median filtering function  
function [data_median,sigma_median] = fcn_DataClean_medianFilterData(data,sigmas)
% Revision history: 2019_11_27 - edited this section to include
% sigma_median as one of the function outputs.

data_median = medfilt1(data,7,'truncate');
sigma_median = medfilt1(sigmas,7,'truncate');

% For debugging:
% figure; plot(data,'b'); hold on; plot(data_median,'c'); plot(highest_expected_data,'r'); plot(lowest_expected_data,'r'); 
 
if 1==0  % This keeps the data, but removes only outliers that are beyond 2-sigma range
    % Calculate bounds
    highest_expected_data = data_median + sigma_median;
    lowest_expected_data = data_median - sigma_median;

    % Find outliers and remove them via the median filter
    out_of_bounds = [...
        find(data>highest_expected_data);
        find(data<lowest_expected_data)];
    
    cleaned = data;
    cleaned(out_of_bounds) = data_median(out_of_bounds);
    data_median = cleaned;
end

% For debugging:
% figure;  plot(cleaned,'b'); hold on; plot(highest_expected_data,'r'); plot(lowest_expected_data,'r');
return

%% Unwrapper of angles
function unwrapped_angle = fcn_DataClean_unwrapAngles(wrapped)
initial_angle = wrapped(1,1);
change_in_angle = [0; diff(wrapped)];
index_jumps = find(change_in_angle>180);
change_in_angle(index_jumps) = change_in_angle(index_jumps)-360;
index_jumps = find(change_in_angle<-180);
change_in_angle(index_jumps) = change_in_angle(index_jumps)+360;
unwrapped_angle = cumsum(change_in_angle) + initial_angle;

% Shift all data up or down 
mean_angle = mean(unwrapped_angle);
good_mean = mod(mean_angle,360);
shift = mean_angle - good_mean;
unwrapped_angle = unwrapped_angle - shift;
return

%% Cleans data with bad indices by nearest neighbor replacement
function new_data = fcn_DataClean_replaceBadIndicesWithNearestGood(data,good_indices)
new_data = data; % Pre-fill new data
all_indices = (1:length(data))'; % Grab all indices
moved_enough = ismember(all_indices,good_indices); % Tag bad ones
for i = 1:length(moved_enough) % Loop through indices
    if 0==moved_enough(i) % Check if bad one
        % Find closest index to stationary one that can work
        distances = (i-good_indices).^2;
        [~,index] = min(distances);
        best_index = good_indices(index);
        new_data(i,1) = data(best_index,1);
    end
end
return



