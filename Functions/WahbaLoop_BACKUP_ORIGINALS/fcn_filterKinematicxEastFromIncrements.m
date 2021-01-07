function FilteredxEast = fcn_filterKinematicxEastFromIncrements(MergedData)

flag_print_progress_to_screen = 1;
flag_make_live_animation_during_regression = 1;

%% Grab data
data1 = MergedData.xEast.Center;
time1 = MergedData.Clocks.targetTimeVector_GPS{5};

data2 = MergedData.xEast_increments.Center;
time2 = MergedData.Clocks.targetTimeVector_GPS{5};

% Create a y vector to represent yaw_true 
% y = interp1(time1, data1,time2);
y = data1;
t = (time2 - time2(1,1));



%% Integrate the yaw angles
xEast_cumulative = cumsum(data2);


%% Global Fit
% Now attempt to calibrate one to the other using the model
% yaw_true = yaw_offset + scale*y_pred + a*t + b*t^2.


% Next, set up the regression
X = [ones(length(xEast_cumulative(:,1)),1),...
    xEast_cumulative,...
    t,...
    t.^2];
%m = X'*X\X'*y;
m = lscov(X,y);
xEast_offset = m(1);
scale = m(2);
a = m(3);
b = m(4);

% Apply the regression
xEast_calc3 = xEast_offset + scale*xEast_cumulative + a*t + b*t.^2;



%% Finally, we implement a rolling window regression
num_seconds = 100;
N_window = round(num_seconds/0.05);

% First, fill the data with prior fit
xEast_calc4 = xEast_calc3;

% Create vectors to store the results, and intialize them
ones_vector = ones(length(xEast_cumulative),1);
xEast_offsets = xEast_offset * ones_vector;
scales      = scale      * ones_vector;
as          = a          * ones_vector;
bs          = b          * ones_vector;
t_full      = t;
y_full      = y;

deltaT = 0.05;
local_t = (0:deltaT:(deltaT*N_window))';
local_t = local_t - local_t(end);
local_t_squared = local_t.^2;
local_indices = (-N_window:0)';

% Next, loop through the raw data,
for i_start = (N_window+1):length(xEast_cumulative)
    

    % Fill in window of indices
    window_indices = i_start+local_indices;
    
    % Fill in local yaw_calc and time vectors
    xEast_local = xEast_cumulative(window_indices,1);
    y_local   =        y_full(window_indices,1);
    
    % Fill data using just this window
    X = [ones(length(xEast_local(:,1)),1),...
        xEast_local,...
        local_t,...
        local_t_squared];
    
    % Do regression if matrix is not singular
    m = lscov(X,y_local);
       
    % Extract results
    xEast_offset = m(1);
    scale = m(2);
    a = m(3);
    b = m(4);
    
    % Predict the yaw at this intex from the regression. Since t=0 in local
    % time, we only need to use the first 2 terms.
    
    % Grab the current yaw_calc for i_start
    xEast_calc_now = xEast_cumulative(i_start,1);
    
    % Do the prediction
    xEast_newcalc_now = xEast_offset + scale*xEast_calc_now;
    
    % Save results
    xEast_calc4(i_start,1)    = xEast_newcalc_now;
    xEast_offsets(i_start,1)  = xEast_offset;
    scales(i_start,1)       = scale;
    as(i_start,1)           = a;
    bs(i_start,1)           = b;

    % Show progress?
    if 1 == flag_print_progress_to_screen
        if mod(i_start,1000)==0
            fprintf('Eval %d of %d  ',i_start,length(xEast_cumulative));
            fprintf('Pred: %f Actual: %f error: %f\n',...
                xEast_newcalc_now,...
                xEast_calc3(i_start),...
                abs(xEast_newcalc_now - xEast_calc3(i_start)));
        end
    end
    
    % Debugging plots via live animation?
    if 1== flag_make_live_animation_during_regression
        if mod(i_start,40)==0
            local_xEast_pred = xEast_offset + scale*xEast_local + a*local_t + b*local_t_squared;
            
            figure(444);
            plot(local_t,y_local,'r',...
                local_t,local_xEast_pred,'b');
            pause(0.01);
        end
    end
    
end

%%  Show the results
if 1==0
    figure(24464);
    clf;
    hold on;
    plot(time1 - time1(1,1), data1,'k');
    plot(time2 - time2(1,1), xEast_cumulative,'r');
    legend('Median Yaw','Calculated Yaw from YawRate');

    % Now attempt to calibrate one to the other using the model
    % yaw_true = scale*y_pred + a*t + b*t^2.
    
    % First, create a y vector to represent yaw_true
    y = interp1(time1, data1,time2);
    t = (time2 - time2(1,1));
    
    % Next, set up the regression
    X = [xEast_cumulative t t.^2];
    m = X'*X\X'*y;
    scale = m(1);
    a = m(2);
    b = m(3);
    
    % Test and show the regreesion
    xEast_calc2 = scale*xEast_cumulative + a*t + b*t.^2;
    plot(time2 - time2(1,1), xEast_calc2,'m');
    legend('Median Yaw','yaw_pred','y=scale*yaw_pred+a*t+b*t^2');
end

error = y_full - xEast_calc4;

FilteredxEast.Center = xEast_calc4;
FilteredxEast.Upper  = xEast_calc4 +2*std(error);
FilteredxEast.Lower  = xEast_calc4 -2*std(error);
FilteredxEast.centiSeconds = 5;

%% Show the fits?
if 1==0
    plot(t, xEast_calc3,'b');
    legend('Median Yaw',...
        'yawPred',...
        'y=scale*yawPred+a*t+b*t^2',...
        'y=yaw_offset+scale*yawPred+a*t+b*t^2');
end

if 1==flag_print_progress_to_screen
    figure(24464);
    clf;
    
    p1 = subplot(2,1,1);
    hold on;
    plot(t_full, y_full,'k');
    plot(t_full, xEast_calc3,'b');
    plot(t_full, xEast_calc4,'r');
    
    legend('Median Yaw',...
        'y=yaw_offset+scale*yawPred+a*t+b*t^2',...
        'window regression');
    
    p2 = subplot(2,1,2);
    plot(t,error,'k');
    xlabel('Time [sec]');
    ylabel('Error [m]');
    
    linkaxes([p1 p2],'x');
    
    figure(7566);
    clf;
    hold on;
    plot(xEast_offsets);
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

