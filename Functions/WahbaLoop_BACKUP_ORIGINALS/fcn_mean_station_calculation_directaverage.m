function [east_gps_mean,north_gps_mean,Num_laps_mean,station_equidistance] = fcn_mean_station_calculation_directaverage(lapData,numLaps)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%%calculate station by laps
station_equidistance=0:1:min([lapData{1}.station(end) lapData{2}.station(end) lapData{3}.station(end) lapData{4}.station(end) ]);
Num_laps_mean= numLaps; % whcih laps you choice to calculate the mean 
north_gps_sum=0;
east_gps_sum=0;
for i_Laps=1:Num_laps_mean
    
        east_gps_interp=interp1(lapData{i_Laps}.station,lapData{i_Laps}.clean_xEast,station_equidistance);
        north_gps_interp=interp1(lapData{i_Laps}.station,lapData{i_Laps}.clean_yNorth,station_equidistance);
%         east_gps_interp=interp1(lapData{i_Laps}.station,lapData{i_Laps}.xEast,station_equidistance);
%         north_gps_interp=interp1(lapData{i_Laps}.station,lapData{i_Laps}.yNorth,station_equidistance);
        
        north_gps_sum=north_gps_interp+north_gps_sum;
        east_gps_sum=east_gps_interp+east_gps_sum;
    
end

% calculate mean ENU data 

east_gps_mean=east_gps_sum/Num_laps_mean;
north_gps_mean=north_gps_sum/Num_laps_mean;

h_fig = figure(167);
set(h_fig,'Name','mean_ENU_in_Laps');
legend_string = ''; % Initialize an empty string
for i_Laps = 1:Num_laps_mean
    plot(lapData{i_Laps}.clean_xEast,lapData{i_Laps}.clean_yNorth,'LineWidth', 1.2);
    hold on;
    legend_string{i_Laps} = sprintf('Lap %d',i_Laps);
end
plot(east_gps_mean,north_gps_mean,'r','LineWidth', 2.4)
grid on; 
xlabel('time [s]') %set  x label 
ylabel('station [m]') % set y label 
title('Plot of mean ENU for Wahba loop'); 
legend(legend_string, 'mean');


end

