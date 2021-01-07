function dout = fcn_alignToGPSTime(d,offset,offset_time_in_ROS)

offset_in_d_time = interp1(offset_time_in_ROS,offset,d.ROS_Time);

%% Creates a pseudo GPS time
dout.GPS_Time               = d.ROS_Time - offset_in_d_time;
dout.GPS_Time_is_emulated   = 1;
% For debugging:
% plot(dout.GPS_Time - min(dout.GPS_Time),'k.'); grid minor;

