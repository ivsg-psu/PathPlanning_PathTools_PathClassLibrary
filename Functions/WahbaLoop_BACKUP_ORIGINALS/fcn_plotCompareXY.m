function fcn_plotCompareXY(raw, clean, fig_number,fig_name)

%% Set up the figure
h_fig = figure(fig_number);
clf;
set(h_fig,'Name',fig_name);
hold on;

%% Fill in data
x = clean.xEast.Center;
y = clean.yNorth.Center;
raw_x = raw.GPS_Hemisphere.xEast;
raw_y =  raw.GPS_Hemisphere.yNorth;

raw_Novatel_x = raw.GPS_Novatel.xEast;
raw_Novatel_y =  raw.GPS_Novatel.yNorth;

raw_x_d = raw_x(raw.GPS_Hemisphere.DGPS_is_active==0);
raw_y_d = raw_y(raw.GPS_Hemisphere.DGPS_is_active==0);

%% Create plots
figure(fig_number);
clf;
hold on;
plot(x,y,'b','Linewidth',1);

plot(raw_x,raw_y,'r.','Linewidth',1);

plot(raw_Novatel_x,raw_Novatel_y,'g','Linewidth',1);

plot(raw_x_d,raw_y_d,'ko','Linewidth',1);


xlabel('xEast [m]') 
ylabel('yNorth [m]'); 
grid minor;

title(cat(2,'Plot of XY from merged data'));
legend( 'clean XY','raw XY','bad gps')

%legend('Central estimate','Upper estimate','Lower estimate');
end