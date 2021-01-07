function [predicted_sigmaDeg] = fcn_predictYawSigma(distances)
% This predicts the variance expected due to yaw angle errors in DGPS data,
% assuming 0.01 meter standard deviation
% Written by S. Brennan
% See script_findRelationBetweenMotionAndYawError.m for derivation

wn = 0.008;
zeta = 0.8;
mag = 100;
mag_num = ((wn^2).^2+(wn*distances).^2).^0.5;
mag_den = ((wn^2 - distances.^2).^2 + (2*zeta*wn*distances).^2).^0.5;
predicted_sigmaDeg = mag*(mag_num./mag_den);

end

