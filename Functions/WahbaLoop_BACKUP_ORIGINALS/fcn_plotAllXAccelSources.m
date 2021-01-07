function fcn_plotAllXAccelSources(d,fig_number,fig_name,varargin)

if nargin>3
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
if any(isnan(d.IMU_ADIS.GPS_Time))
    flag_dataGood = 0;
end
if any(isnan(d.IMU_Novatel.GPS_Time))
    flag_dataGood = 0;
end

% Fill in the time vectors
if flag_dataGood
   offset        = d.IMU_Novatel.GPS_Time(1,1);
   t_IMU_ADIS    = d.IMU_ADIS.GPS_Time - offset;
   t_IMU_Novatel = d.IMU_Novatel.GPS_Time - offset;
else
   offset        = d.IMU_Novatel.ROS_Time(1,1);
   t_IMU_ADIS    = d.IMU_ADIS.ROS_Time;
   t_IMU_Novatel = d.IMU_Novatel.ROS_Time;      
end


%% Create plots

plot(t_IMU_ADIS,   d.IMU_ADIS.XAccel,   'r','Linewidth',1);
plot(t_IMU_Novatel,d.IMU_Novatel.XAccel,'b','Linewidth',1);
    
if nargin == 3
    legend(...
    'Novatel IMU',...
    'ADIS IMU');

elseif nargin > 3
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


grid on;
ylabel('X-accelreation [m/s^2]') % set y label
title('Plot of XAccel from all sources together');

end