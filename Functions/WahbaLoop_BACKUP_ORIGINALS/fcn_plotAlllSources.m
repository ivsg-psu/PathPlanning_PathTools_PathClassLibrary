function fcn_plotAlllSources(d,field_name,sensors_list, fig_number,fig_name,varargin)

if nargin>5
    MergedXAccel= varargin{1};
end

%% Set up the figure
h_fig = figure(fig_number);
clf;
set(h_fig,'Name',fig_name);
hold on;

%% Create time vectors
% Check to see if all dat is from GPS time, or ROS time
% Prefer GPS time as it's ahard standard...
flag_dataGood = 1;
for i_sensor=1:length(sensors_list)
    sensor = d.(sensors_list{i_sensor});
    GPS_Time = sensor.GPS_Time;
    if any(isnan(GPS_Time))
        flag_dataGood = 0;
    end    
end

offset = 0;
% Find the time offset
for i_sensor=1:length(sensors_list)
    sensor = d.(sensors_list{i_sensor});
    if flag_dataGood
        offset_here = sensor.GPS_Time(1,1);
    else
        offset_here = sensor.ROS_Time(1,1);        
    end
    offset = max(offset,offset_here);
end


% Fill in the time and data vectors
for i_sensor=1:length(sensors_list)
    sensor = d.(sensors_list{i_sensor});
    if flag_dataGood
        times{i_sensor} = sensor.GPS_Time - offset; %#ok<AGROW>
    else
        times{i_sensor} = sensor.ROS_Time - offset;     %#ok<AGROW>
    end
    datas{i_sensor} = sensor.(field_name);     %#ok<AGROW>
end

if strcmp(field_name,'yNorth') % in this case, we are doing an XY plot
    for i_sensor=1:length(sensors_list)
        sensor = d.(sensors_list{i_sensor});
        times{i_sensor} = sensor.xEast;
        datas{i_sensor} = sensor.yNorth;
    end
end

%% Create plots
for i_sensor=1:length(sensors_list)
    plot(times{i_sensor},datas{i_sensor},'Linewidth',1);
end

if nargin == 5
    legend(sensors_list);
elseif nargin > 5
    plot(t_IMU_Novatel,MergedXAccel.XAccel_Average,'g','Linewidth',1.5);   
    
    legend(...
    'Novatel IMU',...
    'ADIS IMU',...
    'Average');

end

if flag_dataGood==0
    xlabel('ROS Time [sec]')   % set  x label
else
    xlabel('GPS Time [sec]')   % set  x label
end
ylabel(field_name) % set y label

if strcmp(field_name,'yNorth') % in this case, we are doing an XY plot
    xlabel('xEast [m]');
    ylabel('yNorth [m]');
    xlim([-4700 500]);
    ylim([0 3000]);
    axis square
    grid minor
end

grid minor;
title(cat(2,'Plot of ',field_name,' from all sources together'));

end