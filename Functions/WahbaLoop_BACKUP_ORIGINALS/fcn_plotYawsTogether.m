function fcn_plotYawsTogether(d,fig_number,fig_name)


% To debug this, plot the yaw angle results
h_fig = figure(fig_number);
clf;
set(h_fig,'Name',fig_name);

plot(d.Yaw_deg,'k','Linewidth',2);
hold on;
plot(d.Yaw_deg_from_velocity,'b','Linewidth',2);
plot(d.Yaw_deg_from_position,'c','Linewidth',2);

if all(isnan(d.Yaw_deg))
    legend('Yaw from GPS (none)','Yaw from velocity','Yaw from position');
else
    legend('Yaw from GPS','Yaw from velocity','Yaw from position');
end

grid on;
xlabel('Index [unitless]') %set  x label
ylabel('Yaw angle [deg]') % set y label
title('Plot of yaw angles from all sources together');

end