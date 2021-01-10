function Encoder_RearWheels = fcn_DataClean_loadRawData_Encoder_RearWheels(d,GPS_Novatel,data_source,flag_do_debug)

% This function is used to load the raw data collected with the Penn State Mapping Van.
% This is the Encoder_RearWheels data
% Input Variables:
%      d = raw data from Encoder_RearWheels(format:struct)
%      GPS_Novatel = GPS_Novatel data (format:struct)
%      data_source = the data source of the raw data, can be 'mat_file' or 'database'(format:struct)
%
% Returned Results:
%      Encoder_RearWheels
% Author: Liming Gao
% Created Date: 2020_12_07
%
%
% Updates:
%
% To do lists:
% 1.
%
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% the field name from mat_file is different from database, so we process
% them seperately
if strcmp(data_source,'mat_file')
    Encoder_RearWheels.ROS_Time             = d.Time';
    Encoder_RearWheels.centiSeconds         = 1; % This is sampled every 1 ms
    Encoder_RearWheels.Npoints              = length(Encoder_RearWheels.ROS_Time(:,1));
    Encoder_RearWheels.EmptyVector          = fcn_DataClean_fillEmptyStructureVector(Encoder_RearWheels); % Fill in empty vector (this is useful later)
    %Encoder_RearWheels.GPS_Time             = Encoder_RearWheels.EmptyVector;
    Encoder_RearWheels.deltaT_ROS           = mean(diff(Encoder_RearWheels.ROS_Time));
    %Encoder_RearWheels.deltaT_GPS           = mean(diff(Encoder_RearWheels.GPS_Time));
    Encoder_RearWheels.CountsL              =  d.CountsL';
    Encoder_RearWheels.CountsR              = d.CountsR';
    Encoder_RearWheels.AngularVelocityL     = d.AngularVelocityL';
    Encoder_RearWheels.AngularVelocityR     = d.AngularVelocityR';
    Encoder_RearWheels.DeltaCountsL         = d.DeltaCountsL';
    Encoder_RearWheels.DeltaCountsR         = d.DeltaCountsR';
    %Encoder_RearWheels.DeltaCountsR        = [0; diff(Encoder_RearWheels.CountsR)];
    
    
elseif strcmp(data_source,'database')
    Encoder_RearWheels.ROS_Time             = d.time;
    Encoder_RearWheels.centiSeconds         = 1; % This is sampled every 1 ms
    Encoder_RearWheels.Npoints              = length(Encoder_RearWheels.ROS_Time(:,1));
    Encoder_RearWheels.EmptyVector          = fcn_DataClean_fillEmptyStructureVector(Encoder_RearWheels); % Fill in empty vector (this is useful later)
    %Encoder_RearWheels.GPS_Time             = Encoder_RearWheels.EmptyVector;
    Encoder_RearWheels.deltaT_ROS           = mean(diff(Encoder_RearWheels.ROS_Time));
    %Encoder_RearWheels.deltaT_GPS           = mean(diff(Encoder_RearWheels.GPS_Time));
    Encoder_RearWheels.CountsL              =  d.left_counts;
    Encoder_RearWheels.CountsR              = d.right_counts;
    Encoder_RearWheels.AngularVelocityL     = d.left_angular_velocity;
    Encoder_RearWheels.AngularVelocityR     = d.right_angular_velocity;
    Encoder_RearWheels.DeltaCountsL         = d.left_delta_counts;
    Encoder_RearWheels.DeltaCountsR         = d.right_delta_counts;
    %Encoder_RearWheels.DeltaCountsR        = [0; diff(Encoder_RearWheels.CountsR)];
else
    error('Please indicate the data source')
end


% Calculate the wheel radius, on average
t = GPS_Novatel.ROS_Time;
V = GPS_Novatel.velMagnitude;
t_enc = Encoder_RearWheels.ROS_Time;  % encoder time
w = abs(Encoder_RearWheels.AngularVelocityR);
V_enc = interp1(t,V,t_enc,'nearest','extrap'); % velocity in encoder time
Encoder_RearWheels.RadiusAveR_in_meters = w'*w\(w'*V_enc);

% Use the radius to find the velocity
Encoder_RearWheels.VelocityR            = Encoder_RearWheels.RadiusAveR_in_meters*abs(Encoder_RearWheels.AngularVelocityR);

% Calculate the standard deviation in velocity prediction
error = Encoder_RearWheels.VelocityR - V_enc;
% For debugging
% figure; hist(error,10000);
Encoder_RearWheels.VelocityR_Sigma      = std(error);
Encoder_RearWheels.velMagnitude         = Encoder_RearWheels.VelocityR;
Encoder_RearWheels.velMagnitude_Sigma   = Encoder_RearWheels.VelocityR_Sigma;

%%%% The following was to check to see if tire compressibility was a
%%%% factor...
% % Calculate the wheel radius, with compressibility of tire, k
% % Model is: vel = (r_wheel + v*yawrate*k)*w
% t_yawrate = rawData.IMU_Novatel.ROS_Time;
% yawrate = rawData.IMU_Novatel.ZGyro;
% yawrate_enc = interp1(t_yawrate,yawrate,t_enc,'nearest','extrap');
% w_new = [w V_enc.*yawrate_enc];
% solution = w_new'*w_new\(w_new'*V_enc);
% radius = solution(1);
% k = solution(2);
%
% Encoder_RearWheels.VelocityR_with_compression = (radius + k*V_enc.*yawrate_enc).*abs(Encoder_RearWheels.AngularVelocityR);
% error2 = Encoder_RearWheels.VelocityR_with_compression - V_enc;
% std(error2)
% % For debugging
% figure(4); hist(error2,10000);

if 1==0  % For debugging - to compare steering input to yaw rate (should be linear)
    figure(46464);
    hold on;
    plot(t_enc,Encoder_RearWheels.VelocityR,'r');
    plot(t,V,'k');
    
    %%% The following plots were used to show that the counts are not right
    %     figure(34547);
    %
    %     subplot(2,1,1);
    %     plot(t,V,'k');
    %     p1 = gca;
    %
    %     subplot(2,1,2);
    %     %     plot(...
    %     %         t,Encoder_RearWheels.AngularVelocityR,'r',...
    %     %         t,Encoder_RearWheels.AngularVelocityL,'b');
    %     plot(...
    %         t2,w,'r');
    %     %     plot(...
    %     %         t,Encoder_RearWheels.DeltaCountsR,'r');
    %     p2 = gca;
    %     linkaxes([p1,p2],'x')
end

% Close out the loading process
if flag_do_debug
    % Show what we are doing
    % Grab function name
    st = dbstack;
    namestr = st.name;
    fprintf(1,'\nFinished processing function: %s\n',namestr);
end

return
