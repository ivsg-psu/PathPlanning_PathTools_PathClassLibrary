function fcn_plotArtificialYawFromYawRate(MergedData,cleanAndTimeAlignedData)
% Create artificial yaw angle
data1 = MergedData.Yaw.Yaw_Median;
time1 = cleanAndTimeAlignedData.targetTimeVector_GPS{5};

data2 = MergedData.YawRate.YawRate_Average;
time2 = cleanAndTimeAlignedData.targetTimeVector_GPS{1};

yaw_calc = cumsum(data2)+data1(1,1);

figure(24464);
clf;
hold on;
plot(time1 - time1(1,1), data1,'k');
plot(time2 - time2(1,1), yaw_calc,'r');
legend('Median Yaw','Calculated Yaw from YawRate');


% Now attempt to calibrate one to the other using the model
% yaw_true = scale*y_pred + a*t + b*t^2.

% First, create a y vector to represent yaw_true
y = interp1(time1, data1,time2);
t = (time2 - time2(1,1));

% Next, set up the regression
X = [yaw_calc t t.^2];
m = X'*X\X'*y;
scale = m(1);
a = m(2);
b = m(3);

% Finally, test the regreesion
yaw_calc2 = scale*yaw_calc + a*t + b*t.^2;
plot(time2 - time2(1,1), yaw_calc2,'m');
legend('Median Yaw','yaw_pred','y=scale*yaw_pred+a*t+b*t^2');

%% Second attempt
% Now attempt to calibrate one to the other using the model
% yaw_true = yaw_offset + scale*y_pred + a*t + b*t^2.

% First, create a y vector to represent yaw_true (Skip this, as we did it
% earlier)
% y = interp1(time1, data1,time2);

% Next, set up the regression
X = [ones(length(yaw_calc(:,1)),1),... 
    yaw_calc,... 
    t,...
    t.^2];
%m = X'*X\X'*y;
m = lscov(X,y);
yaw_offset = m(1);
scale = m(2);
a = m(3);
b = m(4);

% Finally, test the regreesion
yaw_calc3 = yaw_offset + scale*yaw_calc + a*t + b*t.^2;
plot(t, yaw_calc3,'b');
legend('Median Yaw',...
    'yawPred',...
    'y=scale*yawPred+a*t+b*t^2',...
    'y=yaw_offset+scale*yawPred+a*t+b*t^2');

%% Finally, we implement a rolling window regression
num_seconds = 10;
N_window = round(num_seconds/0.01);

% First, fill the data with prior fit
yaw_calc4 = yaw_calc3;

% Create vectors to store the results, and intialize them
ones_vector = ones(length(yaw_calc),1);
yaw_offsets = yaw_offset * ones_vector;
scales      = scale      * ones_vector;
as          = a          * ones_vector;
bs          = b          * ones_vector;
t_full      = t;
y_full      = y;
yaw_calc_full = yaw_calc;
m_full      = m;
is_singular = ones_vector;

local_t = (0:deltaT:(deltaT*N_window))';
local_t = local_t - local_t(end);
local_t_squared = local_t.^2;

% Next, loop through the raw data,
for i_start = (N_window+1):(length(yaw_calc_full)-N_window)
    fprintf('Eval %d of %d  ',i_start,(length(yaw_calc_full)-N_window));
    % Grab the current time and yaw_calc for i_start
    t_now = t_full(i_start,1);
    yaw_calc_now = yaw_calc_full(i_start,1);
    
    % Fill in window of indices
    window_indices = (i_start - N_window):(i_start + N_window);
    
    % Fill in local yaw_calc and time vectors
    yaw_calc = yaw_calc_full(window_indices,1);  
    y        =        y_full(window_indices,1);
    
    % Fill data using just this window
    X = [ones(length(yaw_calc(:,1)),1),...
        yaw_calc,...
        local_t,...s
        local_t_squared];
    
    % Do regression if matrix is not singular
    if rank(X'*X)==4
        %m = X'*X\X'*y;
        m = lscov(X,y);
        is_singular(i_start,1) = 0;
        fprintf(1,'\n');
    else
        m = m_full;
        fprintf(1,'Singular\n')
    end
    yaw_offset = m(1);
    scale = m(2);
    a = m(3);
    b = m(4);
    
    % Predict the yaw at this intex from the regression
    yaw_calc_now = yaw_offset +...
        scale*yaw_calc_now +...
        a*t_now + b*t_now.^2;
    
    % Save results
    yaw_calc4(i_start,1)    = yaw_calc_now;
    yaw_offsets(i_start,1)  = yaw_offset;
    scales(i_start,1)       = scale;
    as(i_start,1)           = a;
    bs(i_start,1)           = b;
end

% Show the results
plot(t_full, yaw_calc4,'c');
legend('Median Yaw',...
    'yawPred',...
    'y=scale*yawPred+a*t+b*t^2',...
    'y=yaw_offset+scale*yawPred+a*t+b*t^2',...
    'window regression');



