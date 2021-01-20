function cleanData = fcn_DataClean_cleanRawDataBeforeTimeAlignment(rawDataWithSigmasAndMedianFiltered)
%fcn_cleanRawDataBeforeTimeAlignment - cleans artifacts within the data
%that may be affected by time alignment. For example, yaw angles exhibit
%roll-over effects that need to be removed prior to interpolation, to avoid
%the discontinuities from producing intermediate values that do not
%actually exist.

% Revision history:
% 2019_11_27 - first write of function, moving material out of main code
% area.

flag_do_debug = 1;

%% Let the user know what we are doing
if flag_do_debug
    % Grab function name
    st = dbstack;
    namestr = st.name;

    % Show what we are doing
    fprintf(1,'\nWithin function: %s\n',namestr);
    fprintf(1,'Cleaning indix-based artifacts within data fields, from input structure: %s\n',inputname(1));    
end


cleanData.GPS_Hemisphere     = fcn_DataClean_cleanGPSData(rawDataWithSigmasAndMedianFiltered.GPS_Hemisphere);
cleanData.GPS_Novatel        = fcn_DataClean_cleanGPSData(rawDataWithSigmasAndMedianFiltered.GPS_Novatel);
cleanData.GPS_Garmin         = fcn_DataClean_cleanGPSData(rawDataWithSigmasAndMedianFiltered.GPS_Garmin);
cleanData.IMU_Novatel        = fcn_DataClean_cleanIMUData(rawDataWithSigmasAndMedianFiltered.IMU_Novatel);
cleanData.IMU_ADIS           = fcn_DataClean_cleanIMUData(rawDataWithSigmasAndMedianFiltered.IMU_ADIS);
cleanData.Encoder_RearWheels = fcn_DataClean_cleanEncoderData(rawDataWithSigmasAndMedianFiltered.Encoder_RearWheels);



%% Tell user we are leaving
if flag_do_debug
    % Show what we are doing
    fprintf(1,'Exiting function: %s\n',namestr);   
end

return

