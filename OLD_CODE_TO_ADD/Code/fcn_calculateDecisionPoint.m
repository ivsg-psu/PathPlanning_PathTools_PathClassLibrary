function [x_decision, y_decision] = fcn_calculateDecisionPoint(x1, y1,s1, x2, y2, s2, k)
% fcn_calculateDecisionPoint This function takes in two average path data
% and calculates the point where a user is most likely to make a decision
% Inputs: x1 ---- x coordinates of path 1
%         y1 ---- y coordinates of path 1
%         s1 ---- standard deviation of the probability distribution of user
%         input at each x-y point of path 1
%         x2 ---- x coordinates of path 2
%         y2 ---- y coordinates of path 2
%         s2 ---- standard deviation of the probability distribution of user
%         input at each x-y point of path 2
% Output: x_decision ---- x coordinate of the decision point
%         y_decision ---- y coordinate of the decision point
if nargin < 7
    k = 1;
end

x_decision = [];
y_decision = [];
s_coef = (0.5:0.5:5)';
[PVIP_x, PVIP_y, cluster_idx] = fcn_obtainPVIPs(x1, y1, s1, x2, y2, s2, k);
% for int 1 and 2
%PVIP_x(2) = [];
%PVIP_y(2) = [];
%cluster_idx(2) = [];
length(PVIP_x);
hold on
fcn_plotPathWithVarianceShading(x1, y1, s1, 1, 'r')
fcn_plotPathWithVarianceShading(x2, y2, s2, 1, 'g')
xmin = min(PVIP_x) - 50;
xmax = max(PVIP_x) + 50;
ymin = min(PVIP_y) - 50;
ymax = max(PVIP_y) + 50;
for i = 1:k % k is the number of clusters
    % plot PVIP points
    pvip = plot(PVIP_x(cluster_idx == i), PVIP_y(cluster_idx == i), 'ro');
    
    % fit and plot regression line through PVIP
    p = polyfit(PVIP_x(cluster_idx == i), PVIP_y(cluster_idx == i),1);
    f = polyval(p, xmin:xmax);
    reg = plot(xmin:xmax, f, 'b');
    
    % snap fit PVIP to the regression line
    x_snapfit = (PVIP_x(cluster_idx == i)./p(1) + PVIP_y(cluster_idx == i) - p(2)) ./ (p(1) + 1/p(1));
    y_snapfit = -x_snapfit./p(1) + PVIP_x(cluster_idx == i)./p(1) + PVIP_y(cluster_idx == i);
    snap = plot(x_snapfit, y_snapfit, 'y*');
    
    % Define an arbitary origin
    x_origin = 1053; 
    y_origin = p(1) * x_origin + p(2);
    
    %origin = plot(x_origin, y_origin, 'ko', 'MarkerSize',10);
    % Calculate r sigma
    r_sigma = ((x_snapfit - x_origin).^2 + (y_snapfit - y_origin).^2).^0.5;
    %rad = plot(s_coef, r_sigma(1:10),'ko');
    length(s_coef)
    length(r_sigma)
    
    coef = polyfit(s_coef,r_sigma(1:10),1);
    regression = polyval(coef,s_coef);
    %reg2 = plot(s_coef,regression);
    r_origin = coef(2);
    sigma_dist = abs(regression - r_origin);
    
    % Plot decision making point and decision making region
    dy = r_origin * sin(atan(p(1)));
    dx = r_origin * cos(atan(p(1)));
    x_decision = [x_decision; x_origin + dx];
    y_decision = [y_decision; y_origin + dy];
    dec = plot(x_decision(i), y_decision(i), 'ro', 'MarkerFaceColor','red', 'MarkerSize',10);
    r_decision = sigma_dist .* s_coef;
    theta = linspace(0,2*pi,100);
    for j = 1:4
        % plot(x_decision + r_decision(j)*cos(theta), y_decision + r_decision(j)*sin(theta),'k');
        circle = fill(x_decision(i) + r_decision(j)*cos(theta), y_decision(i) + r_decision(j)*sin(theta), 'k');
        set(circle,'facealpha',0.1);
    end
end
    %axis([xmin xmax ymin ymax]);
    xlabel("x [m]")
    ylabel("y [m]")
    legend([pvip, reg, snap, dec, circle], "Path Variance Intersection Points (PVIPs)",...
        "Regression Line through PVIPs", "Snap-fit PVIPs on Regression Line", ...
        "Decision-Making Point", "Decision-Mkaing Region")
    %{
    % step 1
    subplot(3,3,1)
    grid on
    hold on
    fcn_plotPathWithVarianceShading(x1, y1, s1, 1, 'r')
    fcn_plotPathWithVarianceShading(x2, y2, s2, 1, 'g')    
    pt = plot(PVIP_x(cluster_idx == i), PVIP_y(cluster_idx == i), 'ko');
    legend(pt, "Path Variance Intersection Points (PVIPs)")
    axis([xmin xmax ymin ymax]);
    xlabel("x [m]")
    ylabel("y [m]")
    title("Step 1")
    hold off
    
    % step 2
    subplot(2,2,1)
    grid on
    hold on
    fcn_plotPathWithVarianceShading(x1, y1, s1, 1, 'r')
    fcn_plotPathWithVarianceShading(x2, y2, s2, 1, 'g')    
    plot(PVIP_x(cluster_idx == i), PVIP_y(cluster_idx == i), 'bo');
    p = polyfit(PVIP_x(cluster_idx == i), PVIP_y(cluster_idx == i),1);
    f = polyval(p, xmin:xmax);
    reg = plot(xmin:xmax, f, 'k');
    axis([xmin xmax ymin ymax]);
    xlabel("x [m]")
    ylabel("y [m]")
    title("Step 2")
    legend([reg, pt],"Regression Line through PVIPs" , "Path Variance Intersection Points (PVIPs)")
    hold off
    
    % step 3
    subplot(3,3,3)
    grid on
    hold on
    fcn_plotPathWithVarianceShading(x1, y1, s1, 1, 'r')
    fcn_plotPathWithVarianceShading(x2, y2, s2, 1, 'g')
    x_snapfit = (PVIP_x(cluster_idx == i)./p(1) + PVIP_y(cluster_idx == i) - p(2)) ./ (p(1) + 1/p(1));
    y_snapfit = -x_snapfit./p(1) + PVIP_x(cluster_idx == i)./p(1) + PVIP_y(cluster_idx == i);
    snap = plot(x_snapfit, y_snapfit, 'bo');
    plot(xmin:xmax, f, 'k');
    axis([xmin xmax ymin ymax]);
    xlabel("x [m]")
    ylabel("y [m]")
    title("Step 3")
    legend([snap, reg], "Snap-fit PVIPs on Regression Line", "Regression Line through PVIPs")
    hold off
    
    % step 4
    subplot(3,3,4)
    % Define an arbitary origin
    x_origin = 1000; 
    y_origin = p(1) * x_origin + p(2);
    origin = plot(x_origin, y_origin, 'ko', 'MarkerSize',12);
    
    grid on
    hold on
    fcn_plotPathWithVarianceShading(x1, y1, s1, 1, 'r')
    fcn_plotPathWithVarianceShading(x2, y2, s2, 1, 'g')
    
    snap = plot(x_snapfit, y_snapfit, 'bo');
    plot(xmin:xmax, f, 'k');
    axis([xmin xmax ymin ymax]);
    xlabel("x [m]")
    ylabel("y [m]")
    title("Step 4")
    legend([origin, snap, reg], "Arbitrary Origin", "Snap-fit PVIPs on Regression Line", "Regression Line through PVIPs")
    hold off
    
    % step 5
    subplot(3,3,5)
    grid on
    hold on
    r_sigma = ((x_snapfit - x_origin).^2 + (y_snapfit - y_origin).^2).^0.5;
    rad = plot(s_coef, r_sigma(1:10),'ko');
    coef = polyfit(s_coef,r_sigma(1:10),1)
    regression = polyval(coef,s_coef);
    reg2 = plot(s_coef,regression);
    r_origin = coef(2);
    sigma_dist = abs(regression - r_origin);
    xlabel('Standard Deviation Sigma Variances')
    ylabel('Radius from Arbitrary Origin to Regression Line')
    title("Step 5")
    legend([rad, reg2],"PVIP's Radii", "PVIP's Radii Regression Line")
    
    % step 6
    subplot(3,3,6)
    grid on
    hold on
    rad = plot(s_coef, r_sigma(1:10),'ko');
    reg2 = plot(s_coef,regression);
    xlabel('Standard Deviation Sigma Variances')
    ylabel('Radius from Arbitrary Origin to Regression Line')
    title("Step 5")
    legend([rad, reg2],"PVIP's Radii", "PVIP's Radii Regression Line")
    
    % step 8 Calculate the Decision Region Radii
    subplot(3,3,8)
    grid on
    hold on
    fcn_plotPathWithVarianceShading(x1, y1, s1, 1, 'r')
    fcn_plotPathWithVarianceShading(x2, y2, s2, 1, 'g')
    dy = r_origin * sin(atan(p(1)));
    dx = r_origin * cos(atan(p(1)));
    x_decision = x_origin + dx;
    y_decision = y_origin + dy;
    dec = plot(x_decision, y_decision, 'ko');
    r_decision = sigma_dist .* s_coef;
    theta = linspace(0,2*pi,100);
    for j = 1:4
        % plot(x_decision + r_decision(j)*cos(theta), y_decision + r_decision(j)*sin(theta),'k');
        circle = fill(x_decision + r_decision(j)*cos(theta), y_decision + r_decision(j)*sin(theta), 'k');
        set(circle,'facealpha',0.1);
    end
    axis([xmin xmax ymin ymax]);
    legend([dec, circle], "Decision Point", "Decision-Making Region")
    title("Step 8")
    %}
end

%{
function [x_decision, y_decision] = fcn_calculateDecisionPoint(x1, y1,s1, x2, y2, s2, k)
% fcn_calculateDecisionPoint This function takes in two average path data
% and calculates the point where a user is most likely to make a decision
% Inputs: x1 ---- x coordinates of path 1
%         y1 ---- y coordinates of path 1
%         s1 ---- standard deviation of the probability distribution of user
%         input at each x-y point of path 1
%         x2 ---- x coordinates of path 2
%         y2 ---- y coordinates of path 2
%         s2 ---- standard deviation of the probability distribution of user
%         input at each x-y point of path 2
% Output: x_decision ---- x coordinate of the decision point
%         y_decision ---- y coordinate of the decision point
if nargin < 7
    k = 1;
end

x_decision = [];
y_decision = [];
s_coef = (0.5:0.5:5)';
[PVIP_x, PVIP_y, cluster_idx] = fcn_obtainPVIPs(x1, y1, s1, x2, y2, s2, k);

figure
for i = 1:k
    subplot(2,3,1)
    hold on
    pt = plot(PVIP_x(cluster_idx == i), PVIP_y(cluster_idx == i), 'o');
    
    subplot(2,3,2)
    hold on
    plot(PVIP_x(cluster_idx == i), PVIP_y(cluster_idx == i), 'o');
    p = polyfit(PVIP_x(cluster_idx == i), PVIP_y(cluster_idx == i),1);
    f = polyval(p,PVIP_x(cluster_idx == i));
    reg = plot(PVIP_x(cluster_idx == i), f);
    
    subplot(2,3,3)
    hold on
    x_snapfit = (PVIP_x(cluster_idx == i)./p(1) + PVIP_y(cluster_idx == i) - p(2)) ./ (p(1) + 1/p(1));
    y_snapfit = -x_snapfit./p(1) + PVIP_x(cluster_idx == i)./p(1) + PVIP_y(cluster_idx == i);
    snap = plot(x_snapfit, y_snapfit, 'o');
    plot(PVIP_x(cluster_idx == i), f);
    
    subplot(2,3,4)
    hold on
    % Define an arbitary origin
    x_origin = 1000; 
    y_origin = p(1) * x_origin + p(2);
    r_sigma = ((x_snapfit - x_origin).^2 + (y_snapfit - y_origin).^2).^0.5;
    rad = plot(s_coef, r_sigma(1:10),'o');
    coef = polyfit(s_coef,r_sigma(1:10),1)
    regression = polyval(coef,s_coef);
    reg2 = plot(s_coef,regression);
    r_origin = coef(2);
    sigma_dist = abs(regression - r_origin);
    
    
    subplot(2,3,5)
    hold on
    dy = r_origin * sin(atan(p(1)));
    dx = r_origin * cos(atan(p(1)));
    x_decision = x_origin + dx;
    y_decision = y_origin + dy;
    dec = plot(x_decision, y_decision, 'o');
    r_decision = sigma_dist .* s_coef;
    theta = linspace(0,2*pi,100);
    for j = 1: length(r_decision)
        % plot(x_decision + r_decision(j)*cos(theta), y_decision + r_decision(j)*sin(theta),'k');
        circle = fill(x_decision + r_decision(j)*cos(theta), y_decision + r_decision(j)*sin(theta), 'k');
        set(circle,'facealpha',0.1);
    end
end

% step 1
subplot(2,3,1)
fcn_plotPathWithVarianceShading(x1, y1, s1, 1, 'r')
fcn_plotPathWithVarianceShading(x2, y2, s2, 1, 'g')
legend(pt, "Path Variance Intersection Points (PVIPs)")
axis([min(PVIP_x)-10 max(PVIP_x)+10 min(PVIP_y)-10 max(PVIP_y)+10]);
xlabel("x [m]")
ylabel("y [m]")
title("Step 1")
hold off

% step 2
subplot(2,3,2)
fcn_plotPathWithVarianceShading(x1, y1, s1, 1, 'r')
fcn_plotPathWithVarianceShading(x2, y2, s2, 1, 'g')
axis([min(PVIP_x)-10 max(PVIP_x)+10 min(PVIP_y)-10 max(PVIP_y)+10]);
xlabel("x [m]")
ylabel("y [m]")
title("Step 2")
legend([reg, pt],"Regression Line through PVIPs" , "Path Varia nce Intersection Points (PVIPs)")
hold off

% step 3
subplot(2,3,3)
fcn_plotPathWithVarianceShading(x1, y1, s1, 1, 'r')
fcn_plotPathWithVarianceShading(x2, y2, s2, 1, 'g')
axis([min(PVIP_x)-10 max(PVIP_x)+10 min(PVIP_y)-10 max(PVIP_y)+10]);
xlabel("x [m]")
ylabel("y [m]")
title("Step 3")
legend([snap, reg], "Snap-fit PVIPs on Regression Line", "Regression Line through PVIPs")
hold off

% step 4
subplot(2,3,4)
xlabel('Standard Deviation Sigma Variances')
ylabel('Radius from Arbitrary Origin to Regression Line')
title("Step 4")
legend([rad, reg2],"PVIP's Radii", "PVIP's Radii Regression Line")


% step 5 Identify the Center of the Decision Region]
subplot(2,3,5)
fcn_plotPathWithVarianceShading(x1, y1, s1, 1, 'r')
fcn_plotPathWithVarianceShading(x2, y2, s2, 1, 'g')
axis([min(PVIP_x)-10 max(PVIP_x)+10 min(PVIP_y)-10 max(PVIP_y)+10]);
legend([dec, circle], "Decision Point", "Radius of Decision")
title("Step 5")

end

%}