function fcn_plotArtificialPositionFromIncrementsAndVelocity(MergedData,cleanAndTimeAlignedData)

% Create artificial yaw angle
data1 = cleanAndTimeAlignedData.GPS_Hemisphere.xEast;
time1 = cleanAndTimeAlignedData.Clocks.targetTimeVector_GPS{5};

data2 = MergedData.xEast_increments.Center;
time2 = MergedData.Clocks.targetTimeVector_GPS{5};

xEast_calc = cumsum(data2)+data1(1,1);

figure(24464);
clf;
hold on;
plot(time1 - time1(1,1), data1,'k');
plot(time2 - time2(1,1), xEast_calc,'r');
legend('GPS_Hemisphere xEast','Calculated xEast from Increments');


% Now attempt to calibrate one to the other using the model
% x_true = scale*x_pred + a*t + b*t^2.

% First, create a y vector to represent yaw_true
% y = interp1(time1, data1,time2);
% t = (time2 - time2(1,1));
y = data1;
t = (time2 - time2(1,1));

% Next, set up the regression
X = [xEast_calc t t.^2];
m = X'*X\X'*y;
scale = m(1);
a = m(2);
b = m(3);

% Finally, test the regreesion
xEast_calc2 = scale*xEast_calc + a*t + b*t.^2;
plot(time2 - time2(1,1), xEast_calc2,'m');
legend('GPS_Hemisphere xEast','Calculated xEast from Increments','y=scale*xEast_pred+a*t+b*t^2');

%% Second attempt
% Now attempt to calibrate one to the other using the model
% yaw_true = yaw_offset + scale*y_pred + a*t + b*t^2.

% First, create a y vector to represent yaw_true (Skip this, as we did it
% earlier)
% y = interp1(time1, data1,time2);

% Next, set up the regression
X = [ones(length(xEast_calc(:,1)),1),... 
    xEast_calc,... 
    t,...
    t.^2];
%m = X'*X\X'*y;
m = lscov(X,y);
xEast_offset = m(1);
scale = m(2);
a = m(3);
b = m(4);

% Finally, test the regreesion
xEast_calc3 = xEast_offset + scale*xEast_calc + a*t + b*t.^2;
plot(t, xEast_calc3,'b');
legend('GPS_Hemisphere xEast',...
    'Calculated xEast from Increments',...
    'y=scale*xEast_pred+a*t+b*t^2',...
    'y=yaw_offset+scale*yawPred+a*t+b*t^2');

%% Finally, we implement a rolling window regression
num_seconds = 10;
N_window = round(num_seconds/0.01);

% First, fill the data with prior fit
xEast_calc4 = xEast_calc3;

% Create vectors to store the results, and intialize them
ones_vector = ones(length(xEast_calc),1);
xEast_offsets = xEast_offset * ones_vector;
scales      = scale      * ones_vector;
as          = a          * ones_vector;
bs          = b          * ones_vector;
t_full      = t;
y_full      = y;
xEast_calc_full = xEast_calc;
m_full      = m;
is_singular = ones_vector;

deltaT = 0.01;
local_t = (0:deltaT:(deltaT*N_window))';
local_t = local_t - local_t(end);
local_t_squared = local_t.^2;
local_indices = (-N_window:0)';

% Next, loop through the raw data,
for i_start = (N_window+1):length(xEast_calc_full)
    
    % Grab the current time and yaw_calc for i_start
    t_now = t_full(i_start,1);
    yaw_calc_now = xEast_calc_full(i_start,1);
    
    % Fill in window of indices
    window_indices = i_start+local_indices;
    
    % Fill in local yaw_calc and time vectors
    xEast_calc = xEast_calc_full(window_indices,1);  
    y        =        y_full(window_indices,1);
    
    % Fill data using just this window
    X = [ones(length(xEast_calc(:,1)),1),...
        xEast_calc,...
        local_t,...s
        local_t_squared];

    m = lscov(X,y);
    % Do regression if matrix is not singular

    % Matrix doesn't become singular (much), so not bothering with the code
    % that follows
    %     if rank(X'*X)==4
    %         %m = X'*X\X'*y;
    %         m = lscov(X,y);
    %         is_singular(i_start,1) = 0;
    %         fprintf(1,' ');
    %     else
    %         m = m_full;
    %         fprintf(1,'Singular ')
    %     end
    
    % Extract results
    xEast_offset = m(1);
    scale = m(2);
    a = m(3);
    b = m(4);
    
    % Predict the yaw at this intex from the regression. Since t=0 in local
    % time, we only need to use the first 2 terms.
    xEast_newcalc_now = xEast_offset + scale*yaw_calc_now;
    
    if mod(i_start,100)==0
        fprintf('Eval %d of %d  ',i_start,length(xEast_calc_full));
        fprintf('Pred: %f Actual: %f error: %f\n',...
            xEast_newcalc_now,...
            xEast_calc3(i_start),...
            abs(xEast_newcalc_now - xEast_calc3(i_start)));
    end
    
    % Debugging:    
    if 1==1
        if mod(i_start,2)==0
            local_xEast_pred = xEast_offset + scale*xEast_calc + a*local_t + b*local_t_squared;
            
            range = 100;
            figure(444);
            plot(local_t(end-range:end),y(end-range:end),'r',...
                local_t(end-range:end),local_xEast_pred(end-range:end),'b');
            pause(0.01);
        end
    end    
    % Save results
    xEast_calc4(i_start,1)    = xEast_newcalc_now;
    xEast_offsets(i_start,1)  = xEast_offset;
    scales(i_start,1)       = scale;
    as(i_start,1)           = a;
    bs(i_start,1)           = b;
end

% Show the results

figure(24464);
plot(t_full, xEast_calc4,'c');
legend('Median Yaw',...
    'yawPred',...
    'y=scale*yawPred+a*t+b*t^2',...
    'y=yaw_offset+scale*yawPred+a*t+b*t^2',...
    'window regression');



