%%
% This program is used to plot the mapping van DGPS data collected for the
% Wahba route on 2019_09_17 with the Penn State Mapping Van.
%
% Author: Liming Gao
% Create Date: 2019_09_24
% Modify Date: 2019_10_05

%% Clear the workspace
clear all
%close all

%% Load the data?
try
    temp = tripTime(1);
catch
    load('Route_Wahba.mat');
    Time_start_ROS = Route_WahbaLoop.Hemisphere_DGPS.ROSTime(1);
    Time_start_Hemisphere_GPS = Route_WahbaLoop.Hemisphere_DGPS.GPSTimeOfWeek(1);
    Time_start_Novatel_GPS = Route_WahbaLoop.GPS_Novatel.Seconds(1);
    
    latitude = Route_WahbaLoop.Hemisphere_DGPS.Latitude;
    longitude = Route_WahbaLoop.Hemisphere_DGPS.Longitude;
    altitude = Route_WahbaLoop.Hemisphere_DGPS.Height;
    navMode = Route_WahbaLoop.Hemisphere_DGPS.NavMode;
    rawTime = Route_WahbaLoop.Hemisphere_DGPS.ROSTime-Time_start_ROS;
    velNorth = Route_WahbaLoop.Hemisphere_DGPS.VNorth;    
    velEast = Route_WahbaLoop.Hemisphere_DGPS.VEast;
    velUp = Route_WahbaLoop.Hemisphere_DGPS.VUp;
    velMagnitude = 1*sqrt(velNorth.^2+ velEast.^2); %1.025
    numSatellites = Route_WahbaLoop.Hemisphere_DGPS.NumOfSats;
    xEast = Route_WahbaLoop.Hemisphere_DGPS.xEast;
    yNorth = Route_WahbaLoop.Hemisphere_DGPS.yNorth;
    zUp = Route_WahbaLoop.Hemisphere_DGPS.zUp;
%     GPSWeek = Route_WahbaLoop.Hemisphere_DGPS.GPSWeek;
    GPSTimeOfWeek = Route_WahbaLoop.Hemisphere_DGPS.GPSTimeOfWeek;  %GPS tow (sec) associated with this message
    StdDevResid=Route_WahbaLoop.Hemisphere_DGPS.StdDevResid;
%     AgeOfDiff = Route_WahbaLoop.Hemisphere_DGPS.AgeOfDiff;
%     ExtendedAgeOfDiff=Route_WahbaLoop.Hemisphere_DGPS.ExtendedAgeOfDiff;
    GPSTime_Novatel =Route_WahbaLoop.GPS_Novatel.Seconds;

    %% Define plot attributes
    skinny_width_Line = 1.2;
    thick_width_Line = 2.4;
    very_thick_width_line = 12;
    very_very_thick_width_line = 24;
    % Fill in empty vector (this is useful later)
    empty_large_data_vector = NaN*latitude;

    %% compare the ROS time and hemisphere GPS time 
     h_fig = figure(18568464);
    set(h_fig,'Name','ROS time and GPS time ');
    p1 = subplot(3,1,1);
    plot(rawTime,'b.');
    hold on;
    plot(GPSTimeOfWeek - Time_start_Hemisphere_GPS,'r.');
    plot(GPSTime_Novatel - Time_start_Novatel_GPS,'g--');
    
    grid on;
    xlabel('Index [unitless]') %set  x label
    ylabel('time [second]') % set y label
    
    legend('ROS time', "GPS time")
    
    p2 = subplot(3,1,2);
    plot(diff(rawTime),'b.');
    hold on 
    plot(diff(GPSTimeOfWeek),'r.');
    plot(diff(GPSTime_Novatel),'g--')
    grid on;
    
    p3 = subplot(3,1,3);
    plot(rawTime-GPSTimeOfWeek,'b.');
    hold on 
       
    grid on;
    
    linkaxes([p1,p2,p3],'x')
    
    
    figure(8383838);
    clf;
    p1 = subplot(2,1,1);
    plot(GPSTime_Novatel(1:end-1)-GPSTimeOfWeek,'b.');
    hold on;
    grid on;
    xlabel('Index [unitless]') %set  x label
    ylabel('time [second]') % set y label
    
    
    p2 = subplot(2,1,2);
    hold on;
    grid on;
    plot(GPSTimeOfWeek - Time_start_Hemisphere_GPS,'r.');
    plot(GPSTime_Novatel - Time_start_Novatel_GPS,'g--');
    xlabel('Index [unitless]') %set  x label
    ylabel('time [second]') % set y label
    legend('Hemisphere', 'Novatel');
    linkaxes([p1 p2],'x');
        
end






