function fcn_plotAllSources(d,field_names,sensors_list, fig_number,fig_name,ylabel_string)
% Plots all the data together in the field name, from each sensor listed in
% the sensor list.

% Revision history:
% 2019_10_05 - first write of code by S. Brennan, sbrennan@psu.edu
% 2019_10_20 - revised code to use dynamic fieldnames, rather than eval
% command. Also fixed yaw references throughout, allowed multiple field
% names instead of just one, added ylable, and cleaned up unused code. And
% fixed typo in function name (!).


%% Set up the figure
h_fig = figure(fig_number);
clf;
set(h_fig,'Name',fig_name);
hold on;


%% Decide if GPS or ROS time will be used - depends on sensors listed
% Check to see if all data is from GPS time, or ROS time
% Prefer GPS time as it's ahard standard...
flag_dataGood = 1;
for i_sensor=1:length(sensors_list)
    sensor = d.(sensors_list{i_sensor});
    GPS_Time = sensor.GPS_Time;
    if any(isnan(GPS_Time))
        flag_dataGood = 0;
    end    
end

%% Find the time offset
offset = 0;
for i_sensor=1:length(sensors_list)
    sensor = d.(sensors_list{i_sensor});
    if flag_dataGood
        offset_here = sensor.GPS_Time(1,1);
    else
        offset_here = sensor.ROS_Time(1,1);        
    end
    offset = max(offset,offset_here);
end

%% Grab all the data
% Here we loop through each sensor, and each field, and save the time,
% data, and names of each data we encounter into a cell array for each.
    

if strcmp(field_names,'XYplot') % in this case, we are doing an XY plot
    for i_sensor=1:length(sensors_list)
        sensor = d.(sensors_list{i_sensor});
        times{i_sensor} = sensor.xEast;
        datas{i_sensor} = sensor.yNorth;
        names{i_sensor} = cat(2,sensors_list{i_sensor},', XY plot'); 
    end
else
    % Not an XY plot, so we create data structures
    i_data = 0;
    for i_sensor=1:length(sensors_list)
        sensor = d.(sensors_list{i_sensor});
        for j_field = 1:length(field_names)
            i_data = i_data+1;
            if flag_dataGood
                times{i_data} = sensor.GPS_Time - offset; %#ok<AGROW>
            else
                times{i_data} = sensor.ROS_Time - offset;     %#ok<AGROW>
            end
            
            try
                datas{i_data} = sensor.(field_names{j_field});     %#ok<AGROW>
            catch
                error('Unknown field detected');
            end
            names{i_data} = cat(2,sensors_list{i_sensor},',',field_names{j_field}); %#ok<AGROW>
        end
    end
    
end

%% Create plots
for i_data=1:length(datas)
    plot(times{i_data},datas{i_data},'Linewidth',1);
end
h_legend = legend(names);
set(h_legend, 'Interpreter', 'none');

if flag_dataGood==0
    xlabel('ROS Time [sec]');   % set  x label
else
    xlabel('GPS Time [sec]');   % set  x label
end
h_ylabel = ylabel(ylabel_string); % set y label
set(h_ylabel, 'Interpreter', 'none');

if strcmp(field_names,'XYplot') % in this case, we are doing an XY plot
    xlabel('xEast [m]');
    ylabel('yNorth [m]');
    xlim([-4700 500]);
    ylim([0 3000]);
    axis square;
    grid minor;
end

grid minor;
h_title = title(cat(2,'Plot of ',fig_name,' from all sources together'));
set(h_title, 'Interpreter', 'none');

end