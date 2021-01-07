function [clean_xEast,clean_yNorth] = fcn_Kalman_filter(xEast,yNorth,rawTime,velMagnitude,Yaw_Angle,navMode,StdDevResid)
%KALMAN_FILTER Summary of this function goes here
%  use the velocity and yaw angle to calculate the x- and y-increments
%   Detailed explanation goes here
% 1. xEast and yNorth are the measured EN data 
% 2. rawTime is the time 
% 3. navMode is 6 means the sensor can receive correction signal, 
% 4. velMagnitude is the magnitude of velocity 
% 5. Yaw_Angle is the yaw angle 
% 6. StdDevResid is the GPS measurement standard deviation 

if nargin==5
    navMode = 5; 
    StdDevResid = 0.3;
end

clean_xEast = xEast;
clean_yNorth = yNorth;
Delta_time =[0 diff(rawTime)];

p_xEast = xEast;
p_yNorth = yNorth;
estimate_uncertainty = zeros(1,length(xEast));
measurement_uncertainty = zeros(1,length(xEast));
Kalman_gain = zeros(1,length(xEast));
estimate_uncertainty(1) =0.3^2;

for i =2:length(xEast)-1
    
    if ((navMode(i)== 6)) %6 means the GPS with RTK 
         
        %there is no filter in this mode 
         %std_gps_measurement = 0.02;
         std_gps_measurement = StdDevResid(i);
         p_xEast(i) = xEast(i); 
         p_yNorth(i) =yNorth(i);
         estimate_uncertainty(i) =std_gps_measurement^2;
         
    elseif ((navMode(i)~= 6)) %
         
         %input 
         std_gps_measurement = 0.3;  %
         %std_gps_measurement = StdDevResid(i);
         measurement_uncertainty(i) = std_gps_measurement^2;
         
         %update
         Kalman_gain(i) = estimate_uncertainty(i-1) /(estimate_uncertainty(i-1) + measurement_uncertainty(i));
         p_xEast(i) = p_xEast(i-1) + Kalman_gain(i) *(xEast(i)-p_xEast(i-1));
         p_yNorth(i) = p_yNorth(i-1) + Kalman_gain(i) *(yNorth(i)-p_yNorth(i-1));
         clean_xEast(i) = p_xEast(i);
         clean_yNorth(i) = p_yNorth(i);
         
         estimate_uncertainty(i) = (1+Kalman_gain(i))*estimate_uncertainty(i-1);
         if (estimate_uncertainty(i)  > 100)
             estimate_uncertainty(i) =100;
         end
         %predict
%          p_xEast(i+1) = p_xEast(i) + clean_velEast(i)*Delta_time(i+1);
%          p_yNorth(i+1) = p_yNorth(i) + clean_velNorth(i)*Delta_time(i+1);
         
           p_xEast(i+1) = p_xEast(i) + velMagnitude(i-1)*cosd(Yaw_Angle(i-1))*Delta_time(i+1);
           p_yNorth(i+1) = p_yNorth(i) + velMagnitude(i-1)*sind(Yaw_Angle(i-1))*Delta_time(i+1);
       
         estimate_uncertainty(i+1) =  estimate_uncertainty(i)+0.01;
  
     end
end

end

