function Merged = fcn_mergeAllVelocitySources(d)

%% Set up the data
y = d.Encoder_RearWheels.velMagnitude;
x = d.GPS_Novatel.velMagnitude(:,1);
Novatel_resampled = reshape([x, x, x, x, x]',5*length(x(:,1)),1);
Novatel_resampled = Novatel_resampled(1:length(y(:,1)),1);


data = [y, Novatel_resampled];



%% Find averages
C = mean(data,2);

Merged.Velocity_Average = C;

end