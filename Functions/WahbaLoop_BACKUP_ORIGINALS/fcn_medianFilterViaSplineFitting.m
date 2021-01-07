function fcn_medianFilterViaSplineFitting(xEast,xEast_increment,xEast_increment_pred)

threshold = 10;
flag_plot_variances = 1;
pvalue = 0.00004;
pvalue_tight = pvalue*10;

%%
% Start a spline fit

% % This is the new way
% %spline_fit = csaps(xEast,xEast_increment,pvalue,[]);
% 
% time = (0:0.05:length(xEast(:,1))-1)';
% slope_initial = [0; 0];
% slope_final = [0;0];
% x = time;
% y = [slope_initial, [xEast'; xEast_increment'], slope_final];
% pp = spline(x,y);
% 
% %yy = ppval(pp, linspace(0.5*pi,2*pi,101));
% 
% yy = ppval(pp, linspace(min(x),max(x),101));
% figure(3333);
% 
% plot(yy(1,:),yy(2,:),'-b',y(1,2:5),y(2,2:5),'or')
% axis equal

% % This is the old way
% spline_fit = spaps(xEast,xEast_increment,pvalue);
% points = fnplt(spline_fit,'r');

% This gave a curvy plot, but not what we were expecting
% indices = 1:length(xEast);
% spline_fit = spaps(indices,[xEast';xEast_increment'],pvalue);
% points = fnplt(spline_fit,'r');

spline_fit = spaps(cumsum(xEast_increment_pred),[xEast';xEast_increment'],pvalue);
points = fnplt(spline_fit,'r');

% This is the old way
% spline_fit = spaps(xEast,xEast_increment,pvalue);
% points = fnplt(spline_fit,'r');

% Calculate the distances from the spline for each data point
spline_distance_squared = 0*xEast;
closest_indices = 0*xEast;
for i=1:length(xEast)
    [min_dist,ind_min_dist] = min((points(1,:)-xEast(i,1)).^2 + (points(2,:)-xEast_increment(i,1)).^2);
    spline_distance_squared(i,1) = min_dist;
    closest_indices(i,1) = ind_min_dist;
end

% TO DO: check that closest_indices is sorted (e.g. issorted)

% Remove extreme outliers
distances_smoothed = medfilt1(spline_distance_squared,7);
distances_range = std(distances_smoothed);
bad_indices = find(spline_distance_squared > threshold*distances_range);
xEast_clean = xEast;
xEast_increment_clean = xEast_increment;
spline_indices_to_borrow = closest_indices(bad_indices,1);
xEast_clean(bad_indices,1) = points(1,spline_indices_to_borrow);
xEast_increment_clean(bad_indices,1) = points(2,spline_indices_to_borrow);
   
figure(38383);
clf
hold on;
plot(xEast,xEast_increment,'k')
plot(xEast,xEast_increment_pred,'m')
plot(points(1,:),points(2,:),'r');

plot(xEast(bad_indices),xEast_increment(bad_indices),'ro')
plot(xEast_clean,xEast_increment_clean,'b');
temp = ylim;
ylim([max(-1.5,temp(1,1)) min(1.5,temp(1,2))]);
title(sprintf('Mean bias: %f',mean(xEast_increment_pred - xEast_increment)));

if 1==flag_plot_variances
    figure(575757);
    clf;
    hold on;
    plot(spline_distance_squared,'r');
    plot(distances_smoothed,'b');
    threshold_line = ones(length(spline_distance_squared(:,1)),1)*threshold*distances_range;
    plot(threshold_line,'g');
    ylim([0 1.2*threshold*distances_range]);
end

%% Repeat the process now on cleaned data
% First, grab the data
xEast = xEast_clean;
xEast_increment = xEast_increment_clean;

% Start a spline fit
spline_fit = csaps(xEast,xEast_increment,pvalue_tight,[]);
points = fnplt(spline_fit,'r');

% Calculate the distances from the spline for each data point
spline_distance_squared = 0*xEast;
closest_indices = 0*xEast;
for i=1:length(xEast)
    [min_dist,ind_min_dist] = min((points(1,:)-xEast(i,1)).^2 + (points(2,:)-xEast_increment(i,1)).^2);
    spline_distance_squared(i,1) = min_dist;
    closest_indices(i,1) = ind_min_dist;
end

% Remove extreme outliers
distances_smoothed = medfilt1(spline_distance_squared,7);
distances_range = std(distances_smoothed);
bad_indices = find(spline_distance_squared > threshold*distances_range);
xEast_clean = xEast;
xEast_increment_clean = xEast_increment;
spline_indices_to_borrow = closest_indices(bad_indices,1);
xEast_clean(bad_indices,1) = points(1,spline_indices_to_borrow);
xEast_increment_clean(bad_indices,1) = points(2,spline_indices_to_borrow);
   
figure(3837);
clf
hold on;
plot(xEast,xEast_increment,'k')
plot(points(1,:),points(2,:),'r');
plot(xEast(bad_indices),xEast_increment(bad_indices),'ro')
plot(xEast_clean,xEast_increment_clean,'b');
yrange_sigma = distances_range.^0.5;
ymin = mean(xEast_increment_clean) - 20*yrange_sigma;
ymax = mean(xEast_increment_clean) + 20*yrange_sigma;
ylim([ymin ymax]);

% Plot the distance-squared, so the threshold can be seen
if 1==flag_plot_variances
    figure(5757);
    clf;
    hold on;
    plot(spline_distance_squared,'r');
    plot(distances_smoothed,'b');
    threshold_line = ones(length(spline_distance_squared(:,1)),1)*threshold*distances_range;
    plot(threshold_line,'g');
    ylim([0 1.2*threshold*distances_range]);
end

