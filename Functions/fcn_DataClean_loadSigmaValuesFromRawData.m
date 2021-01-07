function RawDataWithSigmas = fcn_DataClean_loadSigmaValuesFromRawData(rawData)
% fcn_DataClean_loadSigmaValuesFromRawData - ensures that key variables have standard
% deviations defined, if they are left empty during the data loading
% process.

% Revision history:
% 2019_10_10 - first write of function by sbrennan@psu.edu
% 2019_11_27 - added more comments, function header.
% 2020_11_10 - changed file name in prep for DataClean class

flag_do_debug = 1;

if flag_do_debug
    % Grab function name
    st = dbstack;
    namestr = st.name;

    % Show what we are doing
    fprintf(1,'\nWithin function: %s\n',namestr);
    fprintf(1,'Starting iterations through rawData structure to calculate sigma values.\n');    
end

%% An old code can be found here that calculates sigmas based on standard deviations
% rawDataWithSigmas = fcn_DataClean_estimateStatesFromIncrementedStatesViaCumsum(rawDataWithSigmasAndMedianFiltered);


%% Define which fields to calculate sigmas for
fields_to_calculate_sigmas_for = [...
    {'velNorth'},...
    {'velEast'},...
    {'velUp'},...
    {'velMagnitude'},...
    {'Roll_deg'},...
    {'Pitch_deg'},...
    {'Yaw_deg'}...
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


names = fieldnames(rawData); % Grab all the fields that are in rawData structure
for i_data = 1:length(names)
    % Grab the data subfield name
    data_name = names{i_data};
    d = eval(cat(2,'rawData.',data_name));
    
    if flag_do_debug
        fprintf(1,'\n Sensor %d of %d: ',i_data,length(names));
        fprintf(1,'Calculating rawData sigmas for sensor: %s\n',data_name);
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
        if any(strcmp(subFieldName,fields_to_calculate_sigmas_for))

            % Check if Sigma field exists - if it does, do NOT overwrite it
            subFieldNameSigma = cat(2,subFieldName,'_Sigma');
            if any(strcmp(subFieldNameSigma,subfieldNames))
                % The Sigma field already exists, just copy it over then
                dout.(subFieldNameSigma) = d.(subFieldNameSigma);
                if flag_do_debug
                    fprintf(1,' <-- skipped this sigma, already defined\n');
                end
            else
                % The Sigma field does not exist - need to calculate it.
                % Some of the dat may have NaN values, so we need to
                % consider this.
                data = d.(subFieldName);
                real_sigma = fcn_DataClean_calcSigmaNoOutliers(data);
                if isnan(real_sigma)
                    error('Sigma calculation produced an NaN value on field %s of sensor %s. Exiting.',subFieldName,data_name);
                end
                dout.(subFieldNameSigma) = real_sigma;
                
                if flag_do_debug
                    % fprintf(1,' <-- calculated a sigma, has length: %d\n',length(real_sigma(:,1)));
                    fprintf(1,' <-- calculated a sigma\n'); 
                end
            end % Ends the if on whether the subfield already exists
        else 
            % If enter here, then we have a field that does NOT need a
            % sigma calculation. So we can do nothing...
            
            if flag_do_debug
                fprintf(1,'\n');
            end
        end % Ends the if statement to check if subfield is on list
    end % Ends for loop through the subfields

    RawDataWithSigmas.(data_name) = dout; % Save results to main structure
    
     
    
end  % Ends for loop through all sensor names in rawData

if flag_do_debug
    % Show what we are doing
    fprintf(1,'\nFinished processing function: %s\n',namestr);
end

return % Ends the function





%% Subfunctions start here
function real_sigma = fcn_DataClean_calcSigmaNoOutliers(data)
% Some of the data may contain NaN values, hence the use of nanmean and
% nanstd below.

differences = diff(data);
deviations = differences - nanmean(differences);
outlier_sigma = nanstd(deviations);

% Reject outliers
deviations_with_no_outliers = deviations(abs(deviations)<(3*outlier_sigma));
real_sigma = nanstd(deviations_with_no_outliers);

%real_sigma_vector = real_sigma*
return
