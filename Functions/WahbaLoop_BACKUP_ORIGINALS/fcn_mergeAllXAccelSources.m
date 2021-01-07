function MergedXAccel = fcn_mergeAllXAccelSources(d)

flag_do_debug = 1;

%% Set up the data
data1 = d.IMU_ADIS.YAccel;
data2 = d.IMU_Novatel.YAccel;

data1b = d.IMU_ADIS.YAccel;
data2b = d.IMU_Novatel.YAccel;

mag1 = (data1.^2 + data1b.^2).^0.5;
mag2 = (data2.^2 + data2b.^2).^0.5;


if 1==flag_do_debug
    
    figure(3633);
    clf;
    hold on;
    plot(data1,'r');
    plot(data1b,'r--');
    plot(data2,'b');
    plot(data2b,'b--');
    
    figure(4646);
    clf;
    hold on;
    plot(data1,'r');
    plot(data2,'b');
    
    figure(4646333);
    clf;
    hold on;
    plot(mag1,'r');
    plot(mag2,'b');
    

    
end

% Check cross-correllation over 200 milliseconds (sample rate is 10 ms);
[cross_correlation,lags] = xcorr(data1,data2,40);

if 1==flag_do_debug
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
    index_offset = lags(max_correllation_index); %#ok<NASGU>
end

data = [data1, data2];



%% Find averages
C = mean(data,2);

MergedXAccel.XAccel_Average = C;

end