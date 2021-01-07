% Variables used to set up this simulation
closest_space = 0.001; % in meters
furthest_space = 5;    % in meters  (this comes from 100 m/s (200 mph) at 20 Hz)
numDistances = 200;    % How many distances to consider between spaces
one_sigma = 0.01; % Standard deviation in position measurement, in meters
numPoints = 10000; % How many points to use for the Monte Carlo sim

% Define a range of distances to test
distances = logspace(log10(closest_space),log10(furthest_space),numDistances);

% Initialize the standard deviations
std_devs = 0*distances;

for i_distance = 1:numDistances
    d = distances(i_distance); % distance between test points, in meters
    
    pos1 = one_sigma*randn(numPoints,2);  % Generate XY data as Nx2 random matrix
    pos2 = one_sigma*randn(numPoints,2) + [ones(numPoints,1)*d, zeros(numPoints,1)];  % Generate XY data as Nx2 random matrix
    
    delta_pos = pos2-pos1;
    angles_in_deg = atan2d(delta_pos(:,2),delta_pos(:,1));
    std_angle = std(angles_in_deg);
    std_devs(i_distance) = std_angle;
    
    % Plot them both
    figure(565657);
    clf;
    subplot(2,1,1);
    plot(pos1(:,1),pos1(:,2),'r.');
    hold on;
    plot(pos2(:,1),pos2(:,2),'b.');
    axis equal;
    title(sprintf('Spacing is: %f',d));
    
    % Plot the histogram
    subplot(2,1,2);
    hist(angles_in_deg,360);
    xlim([-180 180]);
    hold on;
    temp = ylim;
    plot(-1*2*[std_angle;std_angle],temp','r-');
    plot(2*[std_angle;std_angle],temp','r-');
    
    
    pause(0.01);
end


figure(56756);
clf;
loglog(distances,std_devs);
hold on;
xlabel('S-coordinate distance [m]');
ylabel('Expected one-sigma variance in yaw [deg]');

% Try a fit - note that this is just a low-pass filter in the spatial
% domain!
tau = 0.0085;
mag = 100;
first_order_fit = mag*(1 + (distances/tau).^2).^-0.5;
loglog(distances,first_order_fit);

wn = 0.008;
zeta = 0.8;
mag = 100;
mag_num = ((wn^2).^2+(wn*distances).^2).^0.5;
mag_den = ((wn^2 - distances.^2).^2 + (2*zeta*wn*distances).^2).^0.5;
second_order_fit = mag*(mag_num./mag_den);
loglog(distances,second_order_fit);


figure(55675);
clf;
plot(distances,std_devs);
hold on;
xlabel('S-coordinate distance [m]');
ylabel('Expected one-sigma variance in yaw [deg]');
plot(distances,second_order_fit);

%% Redo plots with fitting function
figure(1616161);
clf;
loglog(distances,std_devs);
hold on;
xlabel('S-coordinate distance [m]');
ylabel('Expected one-sigma variance in yaw [deg]');

[predicted_sigmaDeg] = fcn_predictYawSigma(distances);
loglog(distances,predicted_sigmaDeg,'r.');


figure(5656);
clf;
plot(distances,std_devs);
hold on;
xlabel('S-coordinate distance [m]');
ylabel('Expected one-sigma variance in yaw [deg]');
plot(distances,predicted_sigmaDeg,'r.');
