% This is a script to simulate the variance to expect from a velocity-based
% estimate of GPS.

% Versons: 2019_11_16 - updated to allow for general calculation of
% velocities, for variable sigmas (not just the fixed ones).

% The results show a curve fit of
% the expected variance, in degrees, for various speeds

%% To start, run the simulation with a 1-unit sigma...
% If speed/sigma is constant, then results should be the same across
% situations.

% Variables used to set up this simulation
lowest_speed = 0.01; % in meters
fastest_speed = 1000;    % in meters  (this comes from 100 m/s (200 mph) at 20 Hz)
numSpeeds = 200;    % How many distances to consider between spaces
one_sigma = 1; % Standard deviation in velocity measurement, in m/s
numPoints = 10000; % How many points to use for the Monte Carlo sim

% Define a range of speeds to test
speeds = logspace(log10(lowest_speed),log10(fastest_speed),numSpeeds);

% Initialize the standard deviations
std_devs = 0*speeds;

for i_speed = 1:numSpeeds
    speed = speeds(i_speed)*(2)^0.5; % speed in m/s
    
    speed1 = one_sigma*randn(numPoints,1) + ones(numPoints,1)*speed; 
    speed2 = one_sigma*randn(numPoints,1) + ones(numPoints,1)*speed; 
    
    angles_in_deg = atan2d(speed1(:,1),speed2(:,1))-45;
    std_angle = std(angles_in_deg);
    std_devs(i_speed) = std_angle;
    
    % Plot them both
    figure(565657);
    clf;
    subplot(2,1,1);
    plot(speed1(:,1),speed2(:,1),'r.');
    axis equal;
    title(sprintf('Speed is: %f',speed));
    
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
loglog(speeds,std_devs,'b');
hold on;
xlabel('S-coordinate speed/sigma_speed [unitless]');
ylabel('Expected one-sigma variance in yaw [deg]');
legend('Sim result');

% Try a fit - note that this is just a low-pass filter in the spatial
% domain!
tau = 0.4*one_sigma;
mag = 100;
first_order_fit = mag*(1 + (speeds/tau).^2).^-0.5;
loglog(speeds,first_order_fit,'g');
legend('Sim result','First-order fit');

% Try a 2nd order fit to show that it works better
wn = 0.3;
zeta = 0.8;
mag = 100;
mag_num = ((wn^2).^2+(wn*speeds).^2).^0.5;
mag_den = ((wn^2 - speeds.^2).^2 + (2*zeta*wn*speeds).^2).^0.5;
second_order_fit = mag*(mag_num./mag_den);
loglog(speeds,second_order_fit,'r');
legend('Sim result','First-order fit','Second-order fit');

%% Repeat results, but with different sigma now (to show it still works)

% Variables used to set up this simulation
lowest_speed = 0.01; % in meters
fastest_speed = 1000;    % in meters  (this comes from 100 m/s (200 mph) at 20 Hz)
numSpeeds = 200;    % How many distances to consider between spaces
one_sigma = 0.05; % Standard deviation in velocity measurement, in m/s
numPoints = 10000; % How many points to use for the Monte Carlo sim

% Define a range of speeds to test
speeds = logspace(log10(lowest_speed),log10(fastest_speed),numSpeeds);

% Initialize the standard deviations
std_devs = 0*speeds;

for i_speed = 1:numSpeeds
    speed = speeds(i_speed)*(2)^0.5; % speed in m/s
    
    speed1 = one_sigma*randn(numPoints,1) + ones(numPoints,1)*speed; 
    speed2 = one_sigma*randn(numPoints,1) + ones(numPoints,1)*speed; 
    
    angles_in_deg = atan2d(speed1(:,1),speed2(:,1))-45;
    std_angle = std(angles_in_deg);
    std_devs(i_speed) = std_angle;
    
    % Plot them both
    figure(565657);
    clf;
    subplot(2,1,1);
    plot(speed1(:,1),speed2(:,1),'r.');
    axis equal;
    title(sprintf('Speed is: %f',speed));
    
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
loglog(speeds,std_devs,'b');
hold on;
xlabel('S-coordinate speed/sigma_speed [unitless]');
ylabel('Expected one-sigma variance in yaw [deg]');
legend('Sim result');

% Try a 2nd order fit to show that it works better
wn = 0.3*one_sigma;
zeta = 0.8;
mag = 100;
mag_num = ((wn^2).^2+(wn*speeds).^2).^0.5;
mag_den = ((wn^2 - speeds.^2).^2 + (2*zeta*wn*speeds).^2).^0.5;
second_order_fit = mag*(mag_num./mag_den);
loglog(speeds,second_order_fit,'r');
legend('Sim result','First-order fit','Second-order fit');




%% Do final plots
figure(55675);
clf;
plot(speeds,std_devs,'b');
hold on;
xlabel('S-coordinate speed [m/s]');
ylabel('Expected one-sigma variance in yaw [deg]');
plot(speeds,second_order_fit,'r');
legend('Sim result','Second-order fit');

%% Redo plots with fitting function
figure(1616161);
clf;
loglog(speeds,std_devs);
hold on;
xlabel('S-coordinate speeds [m/s]');
ylabel('Expected one-sigma variance in yaw [deg]');

[predicted_sigmaDeg] = fcn_DataClean_predictYawSigmaFromVelocity(speeds,one_sigma);
loglog(speeds,predicted_sigmaDeg,'r.');
legend('Sim result','Second-order fit');


figure(5656);
clf;
plot(speeds,std_devs);
hold on;
xlabel('S-coordinate speeds [m/s]');
ylabel('Expected one-sigma variance in yaw [deg]');
plot(speeds,predicted_sigmaDeg,'r.');

%% Finally, show the effects of speed




