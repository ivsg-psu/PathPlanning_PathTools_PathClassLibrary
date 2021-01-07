%% Perform a series of plots to see better ways to smooth data 
% Try spline fits on state-space trajectory

%data_to_fit = rawData.GPS_Hemisphere;
data_to_fit = cleanAndTimeAlignedData.GPS_Hemisphere;

% First, grab the data
t_min = min(plottingFlags.TimeZoomPoint)+data_to_fit.GPS_Time(1,1);
t_max = max(plottingFlags.TimeZoomPoint)+data_to_fit.GPS_Time(1,1);

indices_of_interest = find(data_to_fit.GPS_Time>t_min & data_to_fit.GPS_Time<t_max);
xEast = data_to_fit.xEast(indices_of_interest);
xEast_increment = data_to_fit.xEast_increments(indices_of_interest);
fcn_medianFilterViaSplineFitting(xEast,xEast_increment);

%% Testing of generic splines to learn how they work
% Splines use an independent variable, x, to determine how the dependent
% variables should evolve. There can be more than one dependent variable,
% represented by rows in the y matrix. If 2 rows, this allows each row of y
% to represent the respective cartesian coordinates of the spline's output,
% thus creating a smooth xy curve.

%x = pi*[0.5:.5:2]; 
x = [0.5 1 1.5 2];
x = 0.5*[0 1 2 3];

y = [-1 0 -1  0  1  0; 
     0  1  0 -1  0  1];
pp = spline(x,y);
%yy = ppval(pp, linspace(0.5*pi,2*pi,101));
yy = ppval(pp, linspace(min(x),max(x),101));
figure(3333);
plot(yy(1,:),yy(2,:),'-b',y(1,2:5),y(2,2:5),'or')
axis equal

% see: https://www.mathworks.com/matlabcentral/fileexchange/38862-contour-line-smoothing?focused=5249851&tab=function


%%
% Check recursive spline
time_duration = 30; % 30 second window
time_step = 1; % Units are seconds

data_to_fit = cleanAndTimeAlignedData.GPS_Hemisphere;

for t_start = 0:time_step:(data_to_fit.GPS_Time(end)-data_to_fit.GPS_Time(1,1) - time_duration)
    t_min = t_start + data_to_fit.GPS_Time(1,1);
    t_max = t_min + time_duration;
    indices_of_interest = find(data_to_fit.GPS_Time>t_min & data_to_fit.GPS_Time<t_max);
    xEast = data_to_fit.xEast(indices_of_interest);
    xEast_increment = data_to_fit.xEast_increments(indices_of_interest);
    
    fcn_medianFilterViaSplineFitting(xEast,xEast_increment);
    mean_vel = mean(data_to_fit.velMagnitude(indices_of_interest));
    current_time = data_to_fit.GPS_Time(indices_of_interest(1,1),1) - data_to_fit.GPS_Time(1,1);
    title(sprintf('Velocity is: %.2f m/s,Time is: %f.2',mean_vel,current_time));

    %     % Calculate a cubic spline fit to the data
    %     % Format is: csaps(X,Y,p_value,eval_points,weights)
    %     % Pvalue = 1 gives cubic spline, p of 0 produces line fit;
    %     spline_fit = csaps(xEast,xEast_increment,0.4,[]);
    %
    %     figure(38383);
    %     clf;
    %     plot(xEast,xEast_increment,'b');
    %     fnplt(spline_fit,'r');
    %
    %     mean_xEast = mean(xEast);
    %     xlim([mean_xEast-time_duration*30 mean_xEast+time_duration*30]);
    %     ylim([-1.5 1.5]);
    pause(0.01);
end