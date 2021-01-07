function fcn_plotAllSources(d,field_names,sensors_list, fig_number,fig_name,ylabel_string)
% Plots all the data together in the field name, from each sensor listed in
% the sensor list.

% Revision history:
% 2019_10_05 - first write of code by S. Brennan, sbrennan@psu.edu
% 2019_10_20 - revised code to use dynamic fieldnames, rather than eval
% command. Also fixed yaw references throughout, allowed multiple field
% names instead of just one, added ylable, and cleaned up unused code. And
% fixed typo in function name (!).
% 2019_11_26 - fixed plotting so that lines are on top of variance bands.
% Fixed if ~exist bug.

flag_do_debug = 1;

%% Let the user know what we are doing
if flag_do_debug
    % Grab function name
    st = dbstack;
    namestr = st.name;
    
    % Show what we are doing
    fprintf(1,'\nWithin function: %s\n',namestr);
    fprintf(1,'\tPlotting data fields: \n\t\t');
    for i = 1:length(field_names)
        fprintf(1,'%s \t',field_names{i});
    end
    fprintf(1,'\n\tFrom sensors: \n\t\t');
    for i = 1:length(sensors_list)
        fprintf(1,'%s \t',sensors_list{i});
    end
    fprintf(1,'\n');
    
    fprintf(1,'\tPlotting to figure number: %d \n',fig_number);
    fprintf(1,'\tWith figure name: %s \n',fig_name);
    fprintf(1,'\tWith ylabel: %s \n',ylabel_string);
end

%% Check that the variables exist - they might not
% Set plotting flag

% Check sensors first. Loop through entire sensor list, checking flags to
% zero if the sensor is not there.
flags_sensors_exist = ones(length(sensors_list),1);
for i_sensor=1:length(sensors_list)
    if ~isfield(d,sensors_list{i_sensor})
        flags_sensors_exist(i_sensor) = 0;
    end
end
% If the largest number in the sensor flags is zero, there are no sensors
if 0==max(flags_sensors_exist)
    error('None of the sensors exist in the data. Unable to plot.');
end
% Check to see if some of the sensors are missing
missing_sensors = sensors_list(flags_sensors_exist==0);
for i_sensor = 1:length(missing_sensors)
    fprintf(1,'\tWARNING: The following sensors are missing from the data and will be ignored:\n');
    fprintf(1,'\t');
    fprintf(1,'%s\t',missing_sensors{i_sensor});
    fprintf(1,'\n');
end
sensors_list = sensors_list(flags_sensors_exist==1);

%% Grab x-data for time plotting
if ~strcmp(field_names,'XYplot')
    % This is a variable versus time plot, time on x-axis
    
    %% Decide if GPS or ROS time will be used - depends on sensors listed
    % Check to see if all data is from GPS time, or ROS time
    % Prefer GPS time as it's ahard standard...
    flag_data_from_GPS_Time_is_Good = 1;
    for i_sensor=1:length(sensors_list)
        sensorNameString = sensors_list{i_sensor};
        sensor = d.(sensorNameString);
        if ~isfield(sensor,'GPS_Time')
            flag_data_from_GPS_Time_is_Good = 0;
        else
            temp_time = sensor.GPS_Time;
            if any(isnan(temp_time))
                flag_data_from_GPS_Time_is_Good = 0;
            end
        end
    end
    
    % Do we need to check ROS time?
    if 0 == flag_data_from_GPS_Time_is_Good
        % Must use ROS time - recheck here that field exists and there is
        % no NaN. But at this point, throw errors if not.
        for i_sensor=1:length(sensors_list)
            sensorNameString = sensors_list{i_sensor};
            sensor = d.(sensorNameString);
            if ~isfield(sensor,'ROS_Time')
                error('Time plotting requested for variable that appears to have no valid time fields. Sensor is: %s:',sensorNameString);
            else
                temp_time = sensor.ROS_Time;
                if any(isnan(temp_time))
                    error('Time plotting requested for variable that appears to have NaN values. Sensor is: %s:',sensorNameString);
                end
            end
        end
    end
    
    %% Find the time offset
    offset = 0;
    for i_sensor=1:length(sensors_list)
        sensor = d.(sensors_list{i_sensor});
        if flag_data_from_GPS_Time_is_Good
            offset_here = sensor.GPS_Time(1,1);
        else
            offset_here = sensor.ROS_Time(1,1);
        end
        offset = max(offset,offset_here);
    end
    
end  % Ends the if statement for XY plotting


%% Grab all the data
% Here we loop through each sensor, and each field, and save the time,
% data, and names of each data we encounter into a cell array for each.

flag_do_plotting = 1;

if strcmp(field_names,'XYplot') % in this case, we are doing an XY plot
    num_XYplots = 0;
    for i_sensor=1:length(sensors_list)
        sensorNameString = sensors_list{i_sensor};
        sensor = d.(sensorNameString);
        if isfield(sensor,'xEast') && isfield(sensor,'yNorth')
            num_XYplots = num_XYplots + 1;
            times{num_XYplots} = sensor.xEast;
            datas{num_XYplots} = sensor.yNorth;
            names{num_XYplots} = cat(2,sensorNameString,', XY plot');
        else
            fprintf(1,'WARNING: sensor %s is missing either xEast or yNorth field, so must skip this in plotting.\n',sensorNameString);
        end
    end
    if 0 == num_XYplots
        flag_do_plotting = 0;
        fprintf(1,'WARNING: no sensors found with appropriate XY fields, so nothing to plot\n');
    end
else
    % Only enter here if creating a time plot, so we create data structures
    i_data = 0;
    for i_sensor=1:length(sensors_list)
        sensorNameString = sensors_list{i_sensor};
        sensor = d.(sensorNameString);
        for j_field = 1:length(field_names)
            field_name = field_names{j_field};
            field_sigma_name = cat(2,field_name,'_Sigma');
            
            if ~isfield(sensor,field_name)||~isfield(sensor,field_sigma_name)
                fprintf(1,'Sensor: %s is missing fields: %s or %s - skipping this sensor.\n',sensorNameString, field_name,field_sigma_name);
            else % Data should be good
                i_data = i_data+1;
                
                % Grab time
                if flag_data_from_GPS_Time_is_Good
                    times{i_data} = sensor.GPS_Time - offset;
                else
                    times{i_data} = sensor.ROS_Time - offset;
                end
                
                % Grab data
                data = sensor.(field_name);
                datas{i_data} = data;
                
                % Grab median value, which only exists for singles/doubles
                if isa(data,'single')||isa(data,'double')
                    sigma = sensor.(field_sigma_name);
                    data_median = medfilt1(data,7,'truncate');
                    highs{i_data} = data_median + 2*sigma; %#ok<*AGROW>
                    lows{i_data}  = data_median - 2*sigma;
                else
                    highs{i_data} = data;
                    lows{i_data}  = data;
                end
                
                
                names{i_data} = cat(2,sensors_list{i_sensor},',',field_names{j_field});
            end % Ends if check on isfields for field name and sigma name
        end
    end
    
    if 0 == i_data
        flag_do_plotting = 0;
        fprintf(1,'WARNING: no sensors found with appropriate time or data fields, so nothing to plot\n');
    end
    
end

%% Only enter here if plot data looks good
if 1==flag_do_plotting
    %% Set up the figure
    h_fig = figure(fig_number);
    clf;
    set(h_fig,'Name',fig_name);
    hold on;
    grid minor;

    
    %% Create plots
    % Set the colors
    c_temp = colormap('lines');
    plotColor = c_temp(1:length(datas),:); % Choose N colors 
    
    % Plot the data
    for i_data=1:length(datas)
        h_plot{i_data} = plot(times{i_data},datas{i_data},'Linewidth',1);
        set(h_plot{i_data},'Color',plotColor(i_data,:));
    end
    
    % Plot the variance bands
    if ~strcmp(field_names,'XYplot')
        for i_data=1:length(datas)
            newColor = plotColor(i_data,:) * 0.5 + [1 1 1]*0.5;
            fcn_plotVarianceBand(times{i_data},lows{i_data},highs{i_data},newColor);
        end
    end

    % Replot the data, so it is on top of the bands, but so legend is right
    for i_data=1:length(datas)
        h_plot2{i_data} = plot(times{i_data},datas{i_data},'Linewidth',1);
        set(h_plot2{i_data},'Color',plotColor(i_data,:));
    end

    
    % Add the plot legend
    h_legend = legend(names);
    set(h_legend, 'Interpreter', 'none');           
    set(h_legend,'String',names);
    
    % Add the plot's X and Y axis labels
    if strcmp(field_names,'XYplot') % in this case, we are doing an XY plot
        xlabel('xEast [m]');
        ylabel('yNorth [m]');
        axis square;
    else
        if flag_data_from_GPS_Time_is_Good==0
            xlabel('ROS Time [sec]');   % set  x label
        else
            xlabel('GPS Time [sec]');   % set  x label
        end
        h_ylabel = ylabel(ylabel_string); % set y label
        set(h_ylabel, 'Interpreter', 'none');
    end
    
    % Set the title
    h_title = title(cat(2,'Plot of ',fig_name,' from all sources together'));
    set(h_title, 'Interpreter', 'none');
    
    
end % Ends if statement on flag_do_plotting

if 1==flag_do_debug
    fprintf(1,'Exiting function: %s\n',namestr);
end
return;

function fcn_plotVarianceBand(x,low_y,high_y,color)
% See: https://www.mathworks.com/matlabcentral/fileexchange/58262-shaded-area-error-bar-plot
% options.color_area = [128 193 219]./255;    % Blue theme

if 1==1  % This one looks best, but is memory intensive
    % Plotting the result with a patch object
    
    % We want to make a patch object, but it doesn't work with NaN values.
    % So we parse out the groups of data that contain NaN values
    [~,indices] = fcn_parseVectorByNaN(low_y);
    
    for i=1:length(indices)
        x_i = x(indices{i});
        high_y_i = high_y(indices{i});
        low_y_i = low_y(indices{i});
        
        % Now make the patch
        x_vector = [x_i', fliplr(x_i')];
        y_vector = [high_y_i',fliplr(low_y_i')];
        patch = fill(x_vector, y_vector,color);
        set(patch, 'edgecolor', 'none');
        set(patch, 'FaceAlpha', 0.3); % Lower values are more transparent
    end
else % Less memory intensive way is here - it just plots lines
    plot(x, low_y, 'Color',color, 'LineWidth', 1);
    hold on;
    plot(x, high_y, 'Color',color, 'LineWidth', 1);
    %legend('Data','median filtered','95% high','95% low');
end
return