function fcn_plotMergedXY(d,fig_number,fig_name)

%% Set up the figure
h_fig = figure(fig_number);
clf;
set(h_fig,'Name',fig_name);
hold on;

%% Fill in data
x = d.xEast.Center;
y = d.yNorth.Center;


%% Create plots
figure(fig_number);
clf;
hold on;
plot(x,y,'Linewidth',1);

xlabel('xEast [m]') 
ylabel('yNorth [m]'); 
grid minor;

title(cat(2,'Plot of XY from merged data'));
%legend('Central estimate','Upper estimate','Lower estimate');
end