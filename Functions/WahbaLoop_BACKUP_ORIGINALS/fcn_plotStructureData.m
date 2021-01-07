function fcn_plotStructureData(rawDataWithSigmas,plottingFlags)
%fcn_plotStructureData - this function plots the data located in the first
%argument, with plotting flags listed in the second argument.
%   
% Authors: Sean Brennan
% Contact: sbrennan@psu.edu
%
% Revision history:
%   2019_10_20 - First write of this function.
%   2019_10_21 - Minor bug fixes (typos)
%   2019_11_21 - Added plotting capabilty for merged data types
%   2019_11_26 - Added fields to support mergedGPS and mergedINS sensors

%% Grab structure information from input
structureString = inputname(1);  % The name of the structure

%% Let the user know what we are doing
flag_do_debug = 1;

if flag_do_debug
    % Grab function name
    st = dbstack;
    namestr = st.name;

    % Show what we are doing
    fprintf(1,'\nWithin function: %s\n',namestr);
    fprintf(1,'Plotting data from input structure: %s\n',structureString);    
end

%% Define the field structure that is being passed

% Here's the structures that have the above fields. The following defines
% what will be plotted, from what sensors, and what labels to use.
structure_from_field.('Yaw_deg') = [{'GPS_Hemisphere'},{'GPS_Novatel'},{'MergedGPS'}];
fields_from_field.('Yaw_deg') = [{'Yaw_deg'}];
ylabel_from_field.('Yaw_deg') = 'Yaw [deg]'; 

structure_from_field.('Yaw_deg_from_position') = [{'GPS_Hemisphere'},{'GPS_Novatel'}]; 
fields_from_field.('Yaw_deg_from_position') = [{'Yaw_deg_from_position'}];
ylabel_from_field.('Yaw_deg_from_position') = 'Yaw [deg]'; 

structure_from_field.('Yaw_deg_from_velocity') = [{'GPS_Hemisphere'},{'GPS_Novatel'}]; 
fields_from_field.('Yaw_deg_from_velocity') = [{'Yaw_deg_from_velocity'}];
ylabel_from_field.('Yaw_deg_from_velocity') = 'Yaw [deg]'; 

structure_from_field.('All_SingleSensor_Yaw_deg') = [{'GPS_Hemisphere'},{'GPS_Novatel'},{'MergedGPS'}]; 
fields_from_field.('All_SingleSensor_Yaw_deg') = [{'Yaw_deg'},{'Yaw_deg_from_position'},{'Yaw_deg_from_velocity'}];
ylabel_from_field.('All_SingleSensor_Yaw_deg') = 'Yaw [deg]'; 

structure_from_field.('All_AllSensors_Yaw_deg') = [{'GPS_Hemisphere'},{'GPS_Novatel'},{'MergedGPS'}]; 
fields_from_field.('All_AllSensors_Yaw_deg') = [{'Yaw_deg'},{'Yaw_deg_from_position'},{'Yaw_deg_from_velocity'}];
ylabel_from_field.('All_AllSensors_Yaw_deg') = 'Yaw [deg]'; 

structure_from_field.('Yaw_deg_merged') = [{'MergedGPS'}];
fields_from_field.('Yaw_deg_merged') = [{'Yaw_deg'}];
ylabel_from_field.('Yaw_deg_merged') = 'Yaw [deg]'; 

structure_from_field.('All_AllSensors_Yaw_deg_merged') = [{'GPS_Hemisphere'},{'GPS_Novatel'},{'MergedGPS'}];
fields_from_field.('All_AllSensors_Yaw_deg_merged') = [{'Yaw_deg'},{'Yaw_deg_from_position'},{'Yaw_deg_from_velocity'}];
ylabel_from_field.('All_AllSensors_Yaw_deg_merged') = 'Yaw [deg]'; 

% Yawrate plots
structure_from_field.('ZGyro') = [{'IMU_Novatel'},{'IMU_ADIS'},{'MergedIMU'}]; 
fields_from_field.('ZGyro') = [{'ZGyro'}];
ylabel_from_field.('ZGyro') = 'YawRate [rad/sec]'; 

structure_from_field.('All_AllSensors_ZGyro') = [{'IMU_Novatel'},{'IMU_ADIS'},{'MergedIMU'}]; 
fields_from_field.('All_AllSensors_ZGyro') = [{'ZGyro'}];
ylabel_from_field.('All_AllSensors_ZGyro') = 'YawRate [rad/sec]'; 
  
structure_from_field.('ZGyro_merged') = [{'MergedIMU'}];
fields_from_field.('ZGyro_merged') = [{'ZGyro'}];
ylabel_from_field.('ZGyro_merged')  = 'YawRate [rad/sec]'; 

structure_from_field.('All_AllSensors_ZGyro_merged') = [{'IMU_Novatel'},{'IMU_ADIS'},{'MergedIMU'}];
fields_from_field.('All_AllSensors_ZGyro_merged') = [{'ZGyro'}];
ylabel_from_field.('All_AllSensors_ZGyro_merged')  = 'YawRate [rad/sec]'; 


% Velocity plots
structure_from_field.('velMagnitude') = [{'GPS_Hemisphere'},{'GPS_Novatel'},{'Encoder_RearWheels'},{'MergedGPS'}]; 
fields_from_field.('velMagnitude') = [{'velMagnitude'}];
ylabel_from_field.('velMagnitude') = 'Velocity [m/s]'; 

structure_from_field.('All_AllSensors_velMagnitude') = [{'GPS_Hemisphere'},{'GPS_Novatel'},{'Encoder_RearWheels'},{'MergedGPS'}]; 
fields_from_field.('All_AllSensors_velMagnitude') = [{'velMagnitude'}];
ylabel_from_field.('All_AllSensors_velMagnitude') = 'Velocity [m/s]'; 

% XAccel plots
structure_from_field.('XAccel') = [{'IMU_Novatel'},{'IMU_ADIS'}]; 
fields_from_field.('XAccel') = [{'XAccel'}];
ylabel_from_field.('XAccel') = 'X-acceleration [m/s^2]'; 

structure_from_field.('All_AllSensors_XAccel') = [{'IMU_Novatel'},{'IMU_ADIS'}]; 
fields_from_field.('All_AllSensors_XAccel') = [{'XAccel'}];
ylabel_from_field.('All_AllSensors_XAccel') = 'X-acceleration [m/s^2]'; 

% Position increments
structure_from_field.('xEast_increments') = [{'GPS_Novatel'},{'GPS_Hemisphere'},{'MergedGPS'},{'VelocityProjectedByYaw'}]; 
fields_from_field.('xEast_increments') = [{'xEast_increments'}];
ylabel_from_field.('xEast_increments') = 'xEast_increments [m]'; 

structure_from_field.('All_AllSensors_xEast_increments') = [{'GPS_Novatel'},{'GPS_Hemisphere'},{'MergedGPS'},{'VelocityProjectedByYaw'}];
fields_from_field.('All_AllSensors_xEast_increments') = [{'xEast_increments'}];
ylabel_from_field.('All_AllSensors_xEast_increments') = 'xEast_increments [m]'; 

structure_from_field.('yNorth_increments') = [{'GPS_Novatel'},{'GPS_Hemisphere'},{'MergedGPS'},{'VelocityProjectedByYaw'}]; 
fields_from_field.('yNorth_increments') = [{'yNorth_increments'}];
ylabel_from_field.('yNorth_increments') = 'yNorth_increments [m]'; 

structure_from_field.('All_AllSensors_yNorth_increments') = [{'GPS_Novatel'},{'GPS_Hemisphere'},{'MergedGPS'},{'VelocityProjectedByYaw'}];
fields_from_field.('All_AllSensors_yNorth_increments') = [{'yNorth_increments'}];
ylabel_from_field.('All_AllSensors_yNorth_increments') = 'yNorth_increments [m]'; 

% XY plots
structure_from_field.('XYplot') = [{'GPS_Novatel'},{'GPS_Hemisphere'},{'MergedGPS'}]; 
fields_from_field.('XYplot') = [{'XYplot'}];
ylabel_from_field.('XYplot') = 'XYplot [m]'; 

structure_from_field.('All_AllSensors_XYplot') = [{'GPS_Novatel'},{'GPS_Hemisphere'},{'MergedGPS'}];
fields_from_field.('All_AllSensors_XYplot') = [{'XYplot'}];
ylabel_from_field.('All_AllSensors_XYplot') = 'XYplot [m]'; 

% xEast and yNorth plots
structure_from_field.('xEast') = [{'GPS_Novatel'},{'GPS_Hemisphere'},{'MergedGPS'}]; 
fields_from_field.('xEast') = [{'xEast'}];
ylabel_from_field.('xEast') = 'xEast [m]'; 

structure_from_field.('All_AllSensors_xEast') = [{'GPS_Novatel'},{'GPS_Hemisphere'},{'MergedGPS'}];
fields_from_field.('All_AllSensors_xEast') = [{'xEast'}];
ylabel_from_field.('All_AllSensors_xEast') = 'xEast [m]'; 

structure_from_field.('yNorth') = [{'GPS_Novatel'},{'GPS_Hemisphere'},{'MergedGPS'}]; 
fields_from_field.('yNorth') = [{'yNorth'}];
ylabel_from_field.('yNorth') = 'yNorth [m]'; 

structure_from_field.('All_AllSensors_yNorth') = [{'GPS_Novatel'},{'GPS_Hemisphere'},{'MergedGPS'}];
fields_from_field.('All_AllSensors_yNorth') = [{'yNorth'}];
ylabel_from_field.('All_AllSensors_yNorth') = 'yNorth [m]'; 

% DGPS mode plots - DGPS_is_active
structure_from_field.('DGPS_is_active') = [{'GPS_Hemisphere'}]; 
fields_from_field.('DGPS_is_active') = [{'DGPS_is_active'}];
ylabel_from_field.('DGPS_is_active') = 'DGPS_is_active [flag]'; 

structure_from_field.('All_AllSensors_DGPS_is_active') = [{'GPS_Hemisphere'}];
fields_from_field.('All_AllSensors_DGPS_is_active') = [{'DGPS_is_active'}];
ylabel_from_field.('All_AllSensors_DGPS_is_active') = 'DGPS_is_active [flag]'; 


%% Make the plots
for i_field = 1:length(plottingFlags.fields_to_plot)
    fieldString = plottingFlags.fields_to_plot{i_field};  % This is the field that is being plotted
    sensors_that_have_this_field = structure_from_field.(fieldString);        
    fields_to_use = fields_from_field.(fieldString);
    ylabel_string = ylabel_from_field.(fieldString);
   
    % Need to figure out if we are plotting all sensors together, which is
    % one function, or sensors by themselves, which is another function.
    if strncmpi(fieldString,'All',3)
        % Plotting all sensors
        % Plot grouped by sensor?
        if strncmpi(fieldString,'All_SingleSensor_',16)
            % Plotting all fields within a single sensor
            for i_plot = 1:length(sensors_that_have_this_field)
                sensorString = sensors_that_have_this_field{i_plot};
                [figNum] = fcn_plotGenerateFigureNumber(fieldString,sensorString,structureString);
                identifierString = cat(2,sensorString,' Allsensors ',fieldString,' for ',structureString);
                fcn_plotAllSources(rawDataWithSigmas,fields_to_use,[{sensorString}],figNum,identifierString,ylabel_string); %#ok<*NBRAK>
                fcn_autoZoom(fieldString, plottingFlags)
            end
        elseif strncmpi(fieldString,'All_All',7)
            identifierString = cat(2,fieldString,' for ',structureString);
            sensorString = fieldString(16:end);
            [figNum] = fcn_plotGenerateFigureNumber(fieldString,sensorString,structureString);
            fcn_plotAllSources(rawDataWithSigmas,fields_to_use,sensors_that_have_this_field,figNum,identifierString,ylabel_string);
            fcn_autoZoom(fieldString, plottingFlags)
        else
            error('Unknown field detected');
        end
        
    else
        % Plotting just single sensors
        for i_plot = 1:length(sensors_that_have_this_field)
            sensorString = sensors_that_have_this_field{i_plot};
            [figNum] = fcn_plotGenerateFigureNumber(fieldString,sensorString,structureString);
            if any(strcmp(sensorString,plottingFlags.SensorsToPlotIndividually))
                % Check that the field exists before trying to plot
                if ~isfield(rawDataWithSigmas,sensorString)
                    fprintf(1,'WARNING: Sensor %s does not exist as a field within the input data. Skipping...\n',sensorString);
                else
                    % Proceed to plot
                    identifierString = cat(2,sensorString,' ',fieldString,' for ',structureString);
                    fcn_plotVariableWithBounds(rawDataWithSigmas.(sensorString),fields_to_use{1},figNum,identifierString, ylabel_string,identifierString);
                    fcn_autoZoom(fieldString, plottingFlags)
                end
            end
        end
    end 
    

end

%% End the function
if flag_do_debug
    fprintf(1,'\nCompleted function: %s\n',namestr);
end

end

function fcn_autoZoom(fieldString, plottingFlags)
    % Add automatic zoom points?
    if any(strcmp(fieldString,[{'XYplot'},{'All_AllSensors_XYplot'}]))
        if isfield(plottingFlags,'XYZoomPoint')
            axis(plottingFlags.XYZoomPoint);
        end
    else
        if isfield(plottingFlags,'TimeZoomPoint')
            xlim(plottingFlags.TimeZoomPoint);
        end
    end
    
    % Set the y-axis zoom point?
    if isfield(plottingFlags,'ylim')
        if isfield(plottingFlags.ylim,fieldString)
            ylim(plottingFlags.ylim.(fieldString));
        end
    end

    % Set the point style?
    if isfield(plottingFlags,'PlotDataDots')
        if 1 == plottingFlags.PlotDataDots
            h_axis = gca;
            h_lines = get(h_axis,'Children');
            for i=1:length(h_lines)
                set(h_lines(i),'Marker','*')
            end
        end
    end
end
