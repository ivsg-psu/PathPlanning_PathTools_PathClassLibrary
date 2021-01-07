function fcn_searchAllFieldsForNaN(cleanAndTimeAlignedData)

flag_do_debug = 1;

if flag_do_debug
    % Grab function name
    st = dbstack;
    namestr = st.name;

    % Show what we are doing
    fprintf(1,'\nWithin function: %s\n',namestr);
    fprintf(1,'Starting iterations through data structure to ensure there are no NaN.\n');    
end

names = fieldnames(cleanAndTimeAlignedData); % Grab all the fields that are in rawData structure
for i_data = 1:length(names)
    % Grab the data subfield name
    data_name = names{i_data};
    d = cleanAndTimeAlignedData.(data_name);
    
    if flag_do_debug
        fprintf(1,'\n Sensor %d of %d: ',i_data,length(names));
        fprintf(1,'Searching NaN within fields for sensor: %s\n',data_name);
    end
    
    subfieldNames = fieldnames(d); % Grab all the subfields
    for i_subField = 1:length(subfieldNames)
        % Grab the name of the ith subfield
        subFieldName = subfieldNames{i_subField};

        if flag_do_debug
            fprintf(1,'\tProcessing subfield: %s ',subFieldName);
        end
        
        % Check to see if this subField has any NaN
        if ~iscell(d.(subFieldName))
            if any(isnan(d.(subFieldName)))
                if flag_do_debug
                    fprintf(1,' <-- contains an NaN value\n');
                end
            else % No NaNs found
                if flag_do_debug
                    fprintf(1,'\n');
                end
                
            end % Ends the if statement to check if subfield is on list
        end  % Ends if to check if the fiel is a call
    end % Ends for loop through the subfields
    
end  % Ends for loop through all sensor names in rawData
return % Ends the function



