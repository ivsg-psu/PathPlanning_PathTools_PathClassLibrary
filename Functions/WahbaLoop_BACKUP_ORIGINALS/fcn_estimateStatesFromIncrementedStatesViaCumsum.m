function RawDataWithCumsum = fcn_estimateStatesFromIncrementedStatesViaCumsum(rawData)
% This function estimates the full states via cumsum operations on
% incremental states.

flag_do_debug = 1;
flag_verbose_slip = 0;

if flag_do_debug
    % Grab function name
    st = dbstack;
    namestr = st.name;
    
    % Show what we are doing
    fprintf(1,'\nWithin function: %s\n',namestr);
    fprintf(1,'Estimating states via cumsum.\n');
end

fields_to_calculate_cumsum_for = [...
    {'xEast'},...
    {'Dummy placeholder field'}...
    ];


names = fieldnames(rawData); % Grab all the fields that are in rawData structure
for i_data = 1:length(names)
    % Grab the data subfield name
    data_name = names{i_data};
    d = rawData.(data_name);  % This is the input (current) data structure
    dout = d; % This is the output data structure
    
    if flag_do_debug
        fprintf(1,'\n Sensor %d of %d: ',i_data,length(names));
        fprintf(1,'Calculating cumsum values for sensor: %s',data_name);
    end
    
    subfieldNames = fieldnames(d); % Grab all the subfields
    for i_subField = 1:length(subfieldNames)
        % Grab the name of the ith subfield
        subFieldName = subfieldNames{i_subField};

        if flag_do_debug
            fprintf(1,'\tProcessing subfield: %s ',subFieldName);
        end
        
        % Check to see if this subField is in the list
        if any(strcmp(subFieldName,fields_to_calculate_cumsum_for))
            
            if flag_do_debug
                fprintf(1,' <--REPLACING WITH CUMSUM\n');
            end

            % This is the "meat" of this function, where we grab and
            % sum the increments           
            incrementName = cat(2,subFieldName,'_increments');

            % Grab the incremented data
            x = d.(subFieldName);
            x_increments = d.(incrementName);
            clean_x = x(1,1) + cumsum(x_increments);
            dout.(subFieldName) = clean_x;
            
        else
            % If enter here, then we have a field that does NOT need a
            % time check calculation. So we can do nothing...
            if flag_do_debug
                fprintf(1,'\n');
            end
            
        end % Ends the if statement to check if subfield is on list
    end % Ends for loop through the subfields
    
    RawDataWithCumsum.(data_name) = dout; % Save results to main structure
    
    
    
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
