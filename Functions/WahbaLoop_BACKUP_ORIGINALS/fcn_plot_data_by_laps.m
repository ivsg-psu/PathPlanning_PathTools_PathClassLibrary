function fcn_plot_data_by_laps(lapData,numLaps,start_point,fig_number) 
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Plot the raw lla data by laps
%set figure
h_fig = figure(fig_number);
set(h_fig,'Name','Lat_vs_Long_in_Laps');
legend_string = ''; % Initialize an empty string

for i_Laps = 1:numLaps
    plot(lapData{i_Laps}.longitude,lapData{i_Laps}.latitude,'LineWidth', 1);
    hold on;
    legend_string{i_Laps} = sprintf('Lap %d',i_Laps);
end
plot (start_point(2),start_point(1),'ro','MarkerSize',10)
grid on; 
xlabel('Latitude [deg]') %set  x label 
ylabel('Longitude [deg]') % set y label 
title('Plot of raw LLA data by laps'); 
legend(legend_string,'start\_point');

%--------------------------------------------------------------------------------------------------------
% Plot the raw ENU data by laps
h_fig = figure(fig_number+2);
set(h_fig,'Name','ENU_in_Laps');

legend_string = ''; % Initialize an empty string

for i_Laps = 1:numLaps    
    plot(lapData{i_Laps}.xEast,lapData{i_Laps}.yNorth,'LineWidth', 1);
    hold on;
    legend_string{i_Laps} = sprintf('Lap %d',i_Laps);
end
grid on; 
plot (start_point(3),start_point(4),'ro','MarkerSize',10)
xlabel('xEast [m]') %set  x label 
ylabel('yNorth [m]') % set y label 
title('Plot of raw station data by laps'); 
legend(legend_string,'start\_point');

%
end

