function fcn_DataClean_plotAllYawSources(d,fig_number,fig_name,varargin)

if nargin>3
    mergedYaw = varargin{1};
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
if any(isnan(d.GPS_Hemisphere.GPS_Time))
    flag_dataGood = 0;
end
if any(isnan(d.GPS_Novatel.GPS_Time))
    flag_dataGood = 0;
end

% Fill in the time vectors
if flag_dataGood
   offset = d.GPS_Novatel.GPS_Time(1,1);
   t_Hemisphere = d.GPS_Hemisphere.GPS_Time - offset;
   t_Novatel    = d.GPS_Novatel.GPS_Time - offset;
else
    offset = d.GPS_Novatel.ROS_Time(1,1);
   t_Hemisphere = d.GPS_Hemisphere.ROS_Time;
   t_Novatel    = d.GPS_Novatel.ROS_Time;      
end


%% Insert Plots
plot(t_Hemisphere,d.GPS_Hemisphere.Yaw_deg,'k','Linewidth',1);
plot(t_Hemisphere,d.GPS_Hemisphere.Yaw_deg_from_position,'r','Linewidth',1);
plot(t_Hemisphere,d.GPS_Hemisphere.Yaw_deg_from_velocity,'b','Linewidth',1);

plot(t_Novatel,d.GPS_Novatel.Yaw_deg,'k--','Linewidth',1);
plot(t_Novatel,d.GPS_Novatel.Yaw_deg_from_position,'r--','Linewidth',1);
plot(t_Novatel,d.GPS_Novatel.Yaw_deg_from_velocity,'b--','Linewidth',1);

if nargin == 3
        legend(...
        'Hemisphere Yaw from GPS (none)',...
        'Hemisphere Yaw from position',...
        'Hemisphere Yaw from velocity',...
        'Novatel Yaw from GPS',...
        'Novatel Yaw from position',...
        'Novatel Yaw from velocity');
elseif nargin>3
    plot(t_Novatel, mergedYaw.Center,'g','Linewidth',1.5);
    plot(t_Novatel, mergedYaw.Upper,'c','Linewidth',1);
    plot(t_Novatel, mergedYaw.Lower,'c','Linewidth',1);

    legend(...
        'Hemisphere Yaw from GPS (none)',...
        'Hemisphere Yaw from position',...
        'Hemisphere Yaw from velocity',...
        'Novatel Yaw from GPS',...
        'Novatel Yaw from position',...
        'Novatel Yaw from velocity',...
        'Median Yaw',...
        'Upper Bound on Median',...
        'Lower Bound on Median');
end
grid on;
if flag_dataGood
    xlabel('TIme from GPS [sec]') %set  x label
else
    xlabel('Time from ROS [sec]') %set  x label
end
ylabel('Yaw angle [deg]') % set y label
title('Plot of yaw angles from all sources together');

end