function fcn_plotRawYawFromPosition(d,fig_number,fig_name)

x = (1:length(d.Yaw_deg_from_position(:,1)))';

yaw_change_due_to_navMode = d.EmptyVector;
yaw_change_due_to_navMode(d.navMode_bad_indices) = d.Yaw_deg_from_position(d.navMode_bad_indices);

% Calculate median filtered values for yaw
yaw_median = medfilt1(d.Yaw_deg_from_position,7,'truncate');

% Add bounds to the position-based yaw
highest_expected_yaw_angles_in_deg = yaw_median + 2*d.Yaw_deg_from_position_Sigma;
lowest_expected_yaw_angles_in_deg = yaw_median - 2*d.Yaw_deg_from_position_Sigma;


% % Plot the RAW yaw angle results from position
% h_fig = figure(56747844);
% set(h_fig,'Name','Raw_yaw_angle_in_deg_from_position');
%
% plot(Yaw_deg_from_position,'b');
% hold on;
% plot(raw_yaw_angles_in_deg_from_velocity,'g');
% plot(yaw_change_due_to_navMode,'k.');
%



% % Average adjacent values
% average_yaw_angles_in_deg = ...
%     (raw_yaw_angles_in_deg_from_position(1:end-2) + ...
%     raw_yaw_angles_in_deg_from_position(3:end))/2;
% average_yaw_angles_in_deg = [raw_yaw_angles_in_deg_from_position(1,:);
%     average_yaw_angles_in_deg;
%     raw_yaw_angles_in_deg_from_position(end,:)];
%
% %
% highest_expected_yaw_angles_in_deg_ave = average_yaw_angles_in_deg + 2*sigmas_in_Yaw;
% lowest_expected_yaw_angles_in_deg_ave = average_yaw_angles_in_deg - 2*sigmas_in_Yaw;




% To debug this, plot the yaw angle results
h_fig = figure(fig_number);
set(h_fig,'Name',fig_name);
%p1 = subplot(2,1,1);

plot(d.Yaw_deg_from_position,'b');
hold on;
plot(yaw_change_due_to_navMode,'k.');
plot(yaw_median,'c');
legend('Yaw from position','drop-out points','median filtered');
fill_plot(x,lowest_expected_yaw_angles_in_deg,highest_expected_yaw_angles_in_deg)

% 
% plot(highest_expected_yaw_angles_in_deg,'r');
% plot(lowest_expected_yaw_angles_in_deg,'r');
% plot(highest_expected_yaw_angles_in_deg_ave,'c');
% plot(lowest_expected_yaw_angles_in_deg_ave,'c');

grid on;
xlabel('Index [unitless]') %set  x label
ylabel('Yaw angle [deg]') % set y label
title('Plot of raw yaw angle calculations');


% % Plot the navMode
% p2 = subplot(2,1,2);
% plot(navMode);
% linkaxes([p1,p2],'x')
end

function fill_plot(x,low_y,high_y)
plot(x, low_y, 'r', 'LineWidth', 1);
hold on;
plot(x, high_y, 'r', 'LineWidth', 1);
% x2 = [x, fliplr(x)];
% inBetween = [curve1, fliplr(curve2)];
% patch([X fliplr(X)].',[low_y fliplr(high_y)].','r')
% fill(x2, inBetween, 'g');
% % ,'FaceAlpha',.3,'EdgeAlpha',.3)
end