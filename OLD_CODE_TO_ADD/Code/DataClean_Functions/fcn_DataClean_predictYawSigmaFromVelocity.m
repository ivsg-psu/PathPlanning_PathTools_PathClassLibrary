function [predicted_sigmaDeg] = fcn_DataClean_predictYawSigmaFromVelocity(speeds,variance)
% This predicts the variance expected due to yaw angle errors in DGPS data,
% assuming 0.5 meter/sec standard deviation in velocity if in normal GPS,
% or 0.005 if differential GPS mode. 
%
% Written by S. Brennan
% See script_findYawSigmaFromVelocity.m for derivation

% Revision history: 
% 2019_10_01 - first write of code by S.Brennan
% 2019_10_20 - edited code to include real (0.5 m/s) variance as seen in
% data.
% 2019_11_16 - edited code to include any variance 

% if isequal(variance,0.5)
%     wn = 0.2;
%     zeta = 0.8;
%     mag = 70;
% elseif isequal(variance,0.005)
%     wn = 0.008;
%     zeta = 0.8;
%     mag = 19;
% else
%     error('unknown variance for velocity calculations.');
% end

wn = 0.3*variance;
zeta = 0.8;
mag = 100;
    
mag_num = ((wn.^2).^2+(wn.*speeds).^2).^0.5;
mag_den = ((wn.^2 - speeds.^2).^2 + (2*zeta*wn.*speeds).^2).^0.5;
predicted_sigmaDeg = mag*(mag_num./mag_den);

end

