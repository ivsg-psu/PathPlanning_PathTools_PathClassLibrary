function  [start_longitude,start_latitude, start_xEast,start_yNorth] =fcn_plotraw_lla(d,fig_number,fig_name,start_index)

%% Set up the data
delta = 1000;
%% Set up the figure
h_fig = figure(fig_number);
clf;
set(h_fig,'Name',fig_name);
%p1 = subplot(2,1,1);
hold on;

%% Insert plots
plot(d.Longitude,d.Latitude,'b','Linewidth',1);
plot(d.Longitude(start_index),d.Latitude(start_index),'ro','MarkerSize',10);
plot(d.Longitude(1:start_index),d.Latitude(1:start_index),'r','Linewidth',1);
grid on;
xlabel('Longitude[deg]') %set  x label
ylabel('Latitude [deg]') % set y label
title('Plot of raw LL');
start_longitude = d.Longitude(start_index);
start_latitude= d.Latitude(start_index);


%% Set up the figure
h_fig = figure(fig_number+1);
clf;
set(h_fig,'Name',fig_name);
%p1 = subplot(2,1,1);
hold on;

%% Insert plots
plot(d.xEast,d.yNorth,'b','Linewidth',1);
plot(d.xEast(start_index),d.yNorth(start_index),'ro','MarkerSize',10);
plot(d.xEast(start_index:start_index+delta),d.yNorth(start_index:start_index+delta),'r','Linewidth',1);
grid on;
xlabel('xEast [m]') %set  x label
ylabel('yNorth [m]') % set y label
title('Plot of raw EN');

start_xEast =  d.xEast(start_index);
start_yNorth = d.yNorth(start_index);


end
