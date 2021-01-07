function MergedYawRate = fcn_mergeAllYawRateSources(d)

%% Set up the data
data1 = d.IMU_ADIS.ZGyro;
data2 = d.IMU_Novatel.ZGyro;

% Check cross-correllation over 200 milliseconds (sample rate is 10 ms);
[cross_correlation,lags] = xcorr(data1,data2,40);

if 1==0
    figure(363563);
    plot(lags,cross_correlation);
end
[~,max_correllation_index] = max(cross_correlation);
index_offset = lags(max_correllation_index);

% Shift the ADIS data, as I trust the GPS more
N = length(data1(:,1));
if index_offset<0  % Shift foward in time
    index_offset   = -index_offset;
    source_indices = 1:(N-index_offset);
    dest_indices   = (index_offset+1):N;
else  % Shift backward in time
    dest_indices   = 1:(N-index_offset);
    source_indices = (index_offset+1):N;
end
data1(dest_indices) = data1(source_indices);

if 1==0
    % Check the result
    [cross_correlation,lags] = xcorr(data1,data2,40);

    figure(363563);
    plot(lags,cross_correlation);
    [~,max_correllation_index] = max(cross_correlation);
    index_offset = lags(max_correllation_index);
end

data = [data1, data2];



%% Find averages
C = mean(data,2);

MergedYawRate.YawRate_Average = C;

end