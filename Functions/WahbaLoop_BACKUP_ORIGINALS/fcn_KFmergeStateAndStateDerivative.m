function [x_kf,sigma_x] = fcn_KFmergeStateAndStateDerivative(t_x1,x1,x1_Sigma,t_x1dot,x1dot,x1dot_Sigma,nameString)

% Update history
% 2019_11_21 First write of code by sbrennan@psu.edu
% 2019_11_21 Edited to remove yaw-specific behavior


% To-do:
% (as of 2019_11_21) - allow conditional updating rather than time resamplming...

%% Set the flags
flag_plot_open_loop_fit = 0;
flag_plot_results = 0;
flag_do_debug = 1;

%% Let the user know what we are doing
if flag_do_debug
    % Grab function name
    st = dbstack;
    namestr = st.name;
    
    % Show what we are doing
    fprintf(1,'\nWithin function: %s\n',namestr);
    fprintf(1,'\tPerforming KF on variable: %s\n',nameString);
    
    fprintf(1,'\tLength of time for x1: %d\n',length(t_x1));
    fprintf(1,'\tLength of x1: %d\n',length(x1));
    fprintf(1,'\tLength of x1_Sigma: %d\n',length(x1_Sigma));
    
    fprintf(1,'\tLength of time for x1dot: %d\n',length(t_x1dot));
    fprintf(1,'\tLength of x1dot: %d\n',length(x1dot));
    fprintf(1,'\tLength of x1dot_Sigma: %d\n',length(x1dot_Sigma));
    
end



%% Prepare the signals to ensure they have the same time base
% Shift the time (this does not affect the KF)
t_x1 = t_x1 - t_x1(1,1);
t_x1dot = t_x1dot - t_x1dot(1,1);

% Calculate dt to use
centiSeconds_x1 = round(mean(diff(t_x1))*100);
centiSeconds_x1dot = round(mean(diff(t_x1dot))*100);

if centiSeconds_x1 <= centiSeconds_x1dot
    t = t_x1;
    dt = centiSeconds_x1*0.01;
    
    % interpolate the velocity data
    % format is: Vq = interp1(X,V,Xq,METHOD,EXTRAPVAL)
    x1_resampled = x1;
    x1_Sigma_resampled = x1_Sigma;
    x1dot_resampled = interp1(t_x1dot,x1dot,t_x1,'linear','extrap');
    x1dot_Sigma_resampled = interp1(t_x1dot,x1dot_Sigma,t_x1,'nearest','extrap');

else
    t = t_x1dot;
    dt = centiSeconds_x1dot*0.01;
    
    % interpolate the position data
    % format is: Vq = interp1(X,V,Xq,METHOD,EXTRAPVAL)
    x1_resampled = interp1(t_x1,x1,t_x1dot,'linear','extrap');
    x1_Sigma_resampled = interp1(t_x1,x1_Sigma,t_x1dot,'linear','extrap');
    x1dot_resampled = x1dot;
    x1dot_Sigma_resampled = x1dot_Sigma;
end


%% Open loop plotting of dynamics, to see if they fit? (they do!)
if 1==flag_plot_open_loop_fit
    % Calculate position variable from velocity variable using integration
    x1_calc = cumsum(x1dot_resampled)*dt + x1_resampled(1,1);
    x1_Sigma_calc = cumsum(x1dot_Sigma_resampled) + x1_Sigma_resampled(1,1);
    figure(24464);
    clf;
    hold on;
    plot(t, x1_resampled,'k');
    plot(t, x1_resampled+2*x1_Sigma_resampled,'Color',[0.25 0.25 0.25]);
    plot(t, x1_resampled-2*x1_Sigma_resampled,'Color',[0.25 0.25 0.25]);
    plot(t, x1_calc,'r');
    %     plot(t, x1_calc + 2*x1_Sigma_calc,'Color',[1 0.25 0.25]);
    %     plot(t, x1_calc - 2*x1_Sigma_calc,'Color',[1 0.25 0.25]);
    ylim([min(x1_resampled)-4*max(x1_Sigma_resampled) max(x1_resampled)+4*max(x1_Sigma_resampled)]);
    legend('State','State + 2sigma','State - 2sigma','Calculated State from pure Integration');
    %legend('State','State + 2sigma','State - 2sigma','Calculated State from pure Integration','Plus 2 sigma','Minus 2 sigma');
end

%% Kalman filter results using dynamic example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Motion equations %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the state vector with the first values in each state
Xk_prev = [x1_resampled(1,1); x1dot_resampled(1,1)];

% Define motion equation to be: Xk = Phi*Xk_prev + Noise, 
% that is Xk(n) = Xk(n-1) + Vk(n-1) * dt
% Phi represents the dynamics of the system: it is the motion equation's A
% matrix
Phi = [1 dt;
       0  1];

% Initialize the error matrix (or the confidence matrix): P states whether we should
% give more weight to the new measurement or to the model estimate
P = [x1_Sigma_resampled(1,1).^2             0;
                 0 x1dot_Sigma_resampled(1,1).^2];

% M is the measurement matrix.
M = [1 0; 0 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Kalman iteration %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Buffers for later display
Nsamples = length(t);
Xk_buffer = zeros(2,Nsamples);
Xk_buffer(:,1) = Xk_prev;
Sigma_buffer = zeros(Nsamples,1);
Sigma_buffer(1,1) = P(1,1).^0.5;
Z_buffer = zeros(2,Nsamples);

for k=1:Nsamples

    % Q is the process noise covariance. It represents the amount of
    % uncertainty in the model. In our case, we assume that the model is
    % near perfect in the position calculation (only numerical errors
    % allowed in the integration - e.g. 0.0000001), but that the velocity
    % term adds uncertainty. Generally, Q may vary during the process to
    % reflect changing uncertainty in the velocity variable.
    Q = [(0.00000001).^2 0;
        0 x1dot_Sigma_resampled(k,1).^2];
    
    % R is the measurement noise covariance. Generally R will
    % vary between samples.
    R = [x1_Sigma_resampled(k,1).^2 0; 0 x1dot_Sigma_resampled(k,1).^2];
    
    % Z is the measurement vector. In our
    % case, Z = TrueData + RandomGaussianNoise(defined by R)
    Z = [x1_resampled(k,1);x1dot_resampled(k,1)];
    
    % Kalman iteration
    P1 = Phi*P*Phi' + Q;
    S = M*P1*M' + R;

    % K is Kalman gain. If K is large, more weight goes to the measurement.
    % If K is low, more weight goes to the model prediction.
    K = P1*M'*inv(S); %#ok<MINV>
    P = P1 - K*M*P1;

    Xk = Phi*Xk_prev + K*(Z-M*Phi*Xk_prev);
    
    % Save the buffer results
    Z_buffer(:,k) = Z;
    Xk_buffer(:,k) = Xk;
    Sigma_buffer(k) = P(1,1).^0.5;
    

    % For the next iteration
    Xk_prev = Xk;
end

%% Plot resulting graphs 
x_kf = Xk_buffer(1,:)';
sigma_x = Sigma_buffer(:,1);

if 1 == flag_plot_results
    % plot raw measurement data for position:
    figure(4646547);
    clf;
    hold on
    grid on
    
    plot(t,x1_resampled,'r');
    plot(t,x_kf,'b');
    plot(t,x_kf+sigma_x*2,'Color',[0.5 0.5 1]);
    plot(t,x_kf-sigma_x*2,'Color',[0.5 0.5 1]);
    legend('Original data','KF filtered','High','Low');
    xlabel('Time');
    ylabel('Position variable');
        
end



return


%%
% KALMANF - updates a system state vector estimate based upon an
%           observation, using a discrete Kalman filter.
%
% Version 1.0, June 30, 2004
%
% This tutorial function was written by Michael C. Kleder
%
% INTRODUCTION
%
% Many people have heard of Kalman filtering, but regard the topic
% as mysterious. While it's true that deriving the Kalman filter and
% proving mathematically that it is "optimal" under a variety of
% circumstances can be rather intense, applying the filter to
% a basic linear system is actually very easy. This Matlab file is
% intended to demonstrate that.
%
% An excellent paper on Kalman filtering at the introductory level,
% without detailing the mathematical underpinnings, is:
% "An Introduction to the Kalman Filter"
% Greg Welch and Gary Bishop, University of North Carolina
% http://www.cs.unc.edu/~welch/kalman/kalmanIntro.html
%
% PURPOSE:
%
% The purpose of each iteration of a Kalman filter is to update
% the estimate of the state vector of a system (and the covariance
% of that vector) based upon the information in a new observation.
% The version of the Kalman filter in this function assumes that
% observations occur at fixed discrete time intervals. Also, this
% function assumes a linear system, meaning that the time evolution
% of the state vector can be calculated by means of a state transition
% matrix.
%
% USAGE:
%
% s = kalmanf(s)
%
% "s" is a "system" struct containing various fields used as input
% and output. The state estimate "x" and its covariance "P" are
% updated by the function. The other fields describe the mechanics
% of the system and are left unchanged. A calling routine may change
% these other fields as needed if state dynamics are time-dependent;
% otherwise, they should be left alone after initial values are set.
% The exceptions are the observation vectro "z" and the input control
% (or forcing function) "u." If there is an input function, then
% "u" should be set to some nonzero value by the calling routine.
%
% SYSTEM DYNAMICS:
%
% The system evolves according to the following difference equations,
% where quantities are further defined below:
%
% x = Ax + Bu + w  meaning the state vector x evolves during one time
%                  step by premultiplying by the "state transition
%                  matrix" A. There is optionally (if nonzero) an input
%                  vector u which affects the state linearly, and this
%                  linear effect on the state is represented by
%                  premultiplying by the "input matrix" B. There is also
%                  gaussian process noise w.
% z = Hx + v       meaning the observation vector z is a linear function
%                  of the state vector, and this linear relationship is
%                  represented by premultiplication by "observation
%                  matrix" H. There is also gaussian measurement
%                  noise v.
% where w ~ N(0,Q) meaning w is gaussian noise with covariance Q
%       v ~ N(0,R) meaning v is gaussian noise with covariance R
%
% VECTOR VARIABLES:
%
% s.x = state vector estimate. In the input struct, this is the
%       "a priori" state estimate (prior to the addition of the
%       information from the new observation). In the output struct,
%       this is the "a posteriori" state estimate (after the new
%       measurement information is included).
% s.z = observation vector
% s.u = input control vector, optional (defaults to zero).
%
% MATRIX VARIABLES:
%
% s.A = state transition matrix (defaults to identity).
% s.P = covariance of the state vector estimate. In the input struct,
%       this is "a priori," and in the output it is "a posteriori."
%       (required unless autoinitializing as described below).
% s.B = input matrix, optional (defaults to zero).
% s.Q = process noise covariance (defaults to zero).
% s.R = measurement noise covariance (required).
% s.H = observation matrix (defaults to identity).
%
% NORMAL OPERATION:
%
% (1) define all state definition fields: A,B,H,Q,R
% (2) define intial state estimate: x,P
% (3) obtain observation and control vectors: z,u
% (4) call the filter to obtain updated state estimate: x,P
% (5) return to step (3) and repeat
%
% INITIALIZATION:
%
% If an initial state estimate is unavailable, it can be obtained
% from the first observation as follows, provided that there are the
% same number of observable variables as state variables. This "auto-
% intitialization" is done automatically if s.x is absent or NaN.
%
% x = inv(H)*z
% P = inv(H)*R*inv(H')
%
% This is mathematically equivalent to setting the initial state estimate
% covariance to infinity.
%
% SCALAR EXAMPLE (Automobile Voltimeter):
%
% % Define the system as a constant of 12 volts:
% clear s
% s.x = 12;
% s.A = 1;
% % Define a process noise (stdev) of 2 volts as the car operates:
% s.Q = 2^2; % variance, hence stdev^2
% % Define the voltimeter to measure the voltage itself:
% s.H = 1;
% % Define a measurement error (stdev) of 2 volts:
% s.R = 2^2; % variance, hence stdev^2
% % Do not define any system input (control) functions:
% s.B = 0;
% s.u = 0;
% % Do not specify an initial state:
% s.x = nan;
% s.P = nan;
% % Generate random voltages and watch the filter operate.
% tru=[]; % truth voltage
% for t=1:20
%    tru(end+1) = randn*2+12;
%    s(end).z = tru(end) + randn*2; % create a measurement
%    s(end+1)=kalmanf(s(end)); % perform a Kalman filter iteration
% end
% figure
% hold on
% grid on
% % plot measurement data:
% hz=plot([s(1:end-1).z],'r.');
% % plot a-posteriori state estimates:
% hk=plot([s(2:end).x],'b-');
% ht=plot(tru,'g-');
% legend([hz hk ht],'observations','Kalman output','true voltage',0)
% title('Automobile Voltimeter Example')
% hold off

function s = kalmanf(s)

% set defaults for absent fields:
if ~isfield(s,'x'); s.x=nan*z; end
if ~isfield(s,'P'); s.P=nan; end
if ~isfield(s,'z'); error('Observation vector missing'); end
if ~isfield(s,'u'); s.u=0; end
if ~isfield(s,'A'); s.A=eye(length(x)); end
if ~isfield(s,'B'); s.B=0; end
if ~isfield(s,'Q'); s.Q=zeros(length(x)); end
if ~isfield(s,'R'); error('Observation covariance missing'); end
if ~isfield(s,'H'); s.H=eye(length(x)); end

if isnan(s.x)
    % initialize state estimate from first observation
    if diff(size(s.H))
        error('Observation matrix must be square and invertible for state autointialization.');
    end
    s.x = inv(s.H)*s.z; %#ok<MINV>
    s.P = inv(s.H)*s.R*inv(s.H');  %#ok<MINV>
else
    
    % This is the code which implements the discrete Kalman filter:
    
    % Prediction for state vector and covariance:
    s.x = s.A*s.x + s.B*s.u;
    s.P = s.A * s.P * s.A' + s.Q;
    
    % Compute Kalman gain factor:
    K = s.P*s.H'*inv(s.H*s.P*s.H'+s.R); %#ok<MINV>
    
    % Correction based on observation:
    try
        s.x = s.x + K*(s.z-s.H*s.x);
    catch
        disp('debug here');
    end
    
    s.P = s.P - K*s.H*s.P;
    
    % Note that the desired result, which is an improved estimate
    % of the sytem state vector x and its covariance P, was obtained
    % in only five lines of code, once the system was defined. (That's
    % how simple the discrete Kalman filter is to use.) Later,
    % we'll discuss how to deal with nonlinear systems.
    
end

return


function fcn_plotVarianceBand(x,low_y,high_y)
% See: https://www.mathworks.com/matlabcentral/fileexchange/58262-shaded-area-error-bar-plot
% options.color_area = [128 193 219]./255;    % Blue theme

if 1==1  % This one looks best, but is memory intensive
    % Plotting the result
    x_vector = [x', fliplr(x')];
    y_vector = [high_y',fliplr(low_y')];
    patch = fill(x_vector, y_vector,[128 193 219]./255);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', 0.5);
else % Less memory intensive way is here - it just plots lines
    plot(x, low_y, 'r', 'LineWidth', 1);
    hold on;
    plot(x, high_y, 'r', 'LineWidth', 1);
    %legend('Data','median filtered','95% high','95% low');
end
return
