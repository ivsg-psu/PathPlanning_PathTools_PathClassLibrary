function fcn_plotAllVelocitySources(d,fig_number,fig_name,varargin)

if nargin>3
    mergedVelocity = varargin{1};
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
if any(isnan(d.IMU_ADIS.GPS_Time))
    flag_dataGood = 0;
end
if any(isnan(d.IMU_Novatel.GPS_Time))
    flag_dataGood = 0;
end

% Fill in the time vectors
if flag_dataGood
    offset                  = d.GPS_Novatel.GPS_Time(1,1);
    t_Encoder_RearWheels    = d.Encoder_RearWheels.GPS_Time - offset;
    t_GPS_Novatel           = d.GPS_Novatel.GPS_Time - offset;
else
    offset = d.GPS_Hemisphere.ROS_Time(1,1);
    t_Encoder_RearWheels    = d.Encoder_RearWheels.ROS_Time - offset;
    t_GPS_Novatel           = d.GPS_Novatel.ROS_Time - offset;
end

%% Insert Plots
plot(t_Encoder_RearWheels,  d.Encoder_RearWheels.velMagnitude,'b','Linewidth',1);
plot(t_GPS_Novatel,         d.GPS_Novatel.velMagnitude,       'r','Linewidth',1);



if nargin == 3
    plot(t_GPS_Novatel,      d.GPS_Hemisphere.velMagnitude,    'c','Linewidth',1);
    legend(...
        'Novatel Velocity',...
        'Encoder Velocity',...
        'Hemisphere Velocity');
elseif nargin == 4
    plot(t_Encoder_RearWheels , mergedVelocity.Velocity_Average,'g','Linewidth',1.5);

    legend(...
    'Novatel Velocity',...
    'Encoder Velocity',...
    'Average Velocity');

end

% Label depends on data source
if flag_dataGood==0
    xlabel('ROS Time [sec]')   % set  x label
else
    xlabel('GPS Time [sec]')   % set  x label
end

% Standard labels
grid on;
ylabel('Velocity [m/s]') % set y label
title('Plot of velocities from all sources together');

end