function fcn_plotYawRate(d,fig_number,fig_name)

%% Set up the data
x = (1:length(d.ZGyro(:,1)))';

yawRate = d.ZGyro;

% Calculate median filtered values for yawRate
yawRate_median = medfilt1(yawRate,7,'truncate');

% Add bounds to the position-based yawRate
highest_expected_yawRate = yawRate_median + 2*d.ZGyro_Sigma;
lowest_expected_yawRate  = yawRate_median - 2*d.ZGyro_Sigma;

%% Set up the figure
h_fig = figure(fig_number);
clf;
set(h_fig,'Name',fig_name);
%p1 = subplot(2,1,1);
hold on;

%% Insert plots
plot(yawRate,'b','Linewidth',1);
plot(yawRate_median,'c');
fcn_plotVarianceBand(x,lowest_expected_yawRate,highest_expected_yawRate)

grid on;
xlabel('Index [unitless]') %set  x label
ylabel('Yaw Rate [deg/sec]') % set y label
title('Plot of yaw rate from IMU');

end

%% Function to plot the band of variance

function fcn_plotVarianceBand(x,low_y,high_y)
% See: https://www.mathworks.com/matlabcentral/fileexchange/58262-shaded-area-error-bar-plot
% options.color_area = [128 193 219]./255;    % Blue theme

% Plotting the result
x_vector = [x', fliplr(x')];
y_vector = [high_y',fliplr(low_y')];
patch = fill(x_vector, y_vector,[128 193 219]./255);
set(patch, 'edgecolor', 'none');
set(patch, 'FaceAlpha', 0.5);
legend('YawRate','median filtered','95% range');   

%%% Less memory intensive way is here - it just plots lines
% plot(x, low_y, 'r', 'LineWidth', 1);
% hold on;
% plot(x, high_y, 'r', 'LineWidth', 1);
%legend('YawRate','median filtered','95% high','95% low');

end