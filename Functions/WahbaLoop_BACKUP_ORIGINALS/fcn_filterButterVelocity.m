function FilteredYaw = fcn_filterButterVelocity(MergedData,cleanAndTimeAlignedData)

%10-27-2019 edited by liming, design the filter

%% Grab data
%data1 =  MergedData.Velocity.Velocity_Average;
data1 = MergedData.velMagnitude.Center;
time11 = cleanAndTimeAlignedData.Clocks.targetTimeVector_GPS{5};
time1 = time11-time11(1);

%% Design the filter

% FFT
Fs = 20;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 20000;             % Length of signal

Y = fft(data1);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

figure(2)
plot(f,P1) 
title('Velocity Mag Single-Sided Amplitude Spectrum')
xlabel('f (Hz)')
ylabel('|P1(f)|')
xlim([0 1])
grid on;



Fc = 0.3; % Cutting frequency(Hz)
T = 1/Fs; % Sampling period

wn= Fc/(Fs/2); % Normalized cutting frequency
[b,a]=butter(2,wn); % Find parameters for 2nd order butterworth filter
[b2,a2] = butter(4,wn); % Find parameters for 4th order butterworth filter


filtered_vel = filter(b,a,data1); % Filtered signal using 2nd butterworth filter
y2 = filter(b2,a2,data1); % Filtered signal using 4th butterworth filter

figure(3)
hold on 
plot(time1,data1,'r'); % plot mixture of signal and noise
plot(time1,filtered_vel,'b'); % plot filtered signal with 2nd order butter filter
plot(time1,y2,'k'); % plot filtered signal with 4th butter filter
legend('Original signal','2nd Butterworth Filtered Signal','4th Butterworth Filtered Signal')

title('Butter filter result')

grid on


error = data1 - filtered_vel;

FilteredYaw.Center = filtered_vel;
FilteredYaw.Upper  = filtered_vel +2*std(error);
FilteredYaw.Lower  = filtered_vel -2*std(error);
FilteredYaw.centiSeconds = 5;

disp('STOPPED HERE');

%% Show the fits?
if 1==0
    plot(t, yaw_calc3,'b');
    legend('Median Yaw',...
        'yawPred',...
        'y=scale*yawPred+a*t+b*t^2',...
        'y=yaw_offset+scale*yawPred+a*t+b*t^2');
end

if 1==0
    figure(24464);
    clf;
    
    p1 = subplot(2,1,1);
    hold on;
    plot(t_full, y_full,'k');
    plot(t_full, yaw_calc3,'b');
    plot(t_full, yaw_calc4,'r');
    
    legend('Median Yaw',...
        'y=yaw_offset+scale*yawPred+a*t+b*t^2',...
        'window regression');
    
    p2 = subplot(2,1,2);
    plot(t,error,'k');
    xlabel('Time [sec]');
    ylabel('Error [deg]');
    
    linkaxes([p1 p2],'x');
    
    figure(7566);
    clf;
    hold on;
    plot(yaw_offsets);
    xlabel('index');
    ylabel('yaw_offset values from regression');
    
    
    figure(74747);
    clf;
    hold on;
    plot(scales);
    xlabel('index');
    ylabel('scale values from regression');
    
    figure(343434);
    clf;
    hold on;
    plot(as);
    xlabel('index');
    ylabel('a values from regression');
    
    figure(766768);
    clf;
    hold on;
    plot(bs);
    xlabel('index');
    ylabel('b values from regression');
    
end

