function fcn_plotArtificialYawRateFromYaw(MergedData,cleanAndTimeAlignedData)
% Create artificial yaw rate
data = MergedData.Yaw.Yaw_Median;
time = cleanAndTimeAlignedData.targetTimeVector_GPS{5};

[b,a] = butter(2,10/50);
dataFiltered = filtfilt(b,a,data);
diff_yaw = [0; diff(dataFiltered)];
diff_t   = [0.05; diff(time)];
yaw_rate = diff_yaw./diff_t;
figure(254);
hold on;
plot(time - time(1,1), yaw_rate,'b');
end