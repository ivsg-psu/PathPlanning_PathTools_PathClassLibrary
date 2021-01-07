function velTrans = fcn_velocityTransform(velArray, angularPos, angularVel, pos, varargin)
% Finds velocity of hemisphere-{2} in reference frame-{0} 
% given velocity/orientation of novatel-{1} in reference frame-{0},  
% angular velocity of novatel-{1} in novatel frame-{1}, and 
% position of hemisphere-{2} wrt novatel-{1} in novatel frame-{1}
%
% pos = 0.0254*[21, 0, 0] Hemisphere WRT Novatel-GPS in Novatel-GPS
% pos = 0.0254*[8.625, 0, 0] Hemisphere WRT Novatel-IMU in Novatel-IMU
% pos = 0.0254*[18, -1, 0] Hemisphere WRT ADIS-IMU in ADIS-IMU
%
% 'velArray' contains velocity of {1} in {0}
% 'angularPos' contains orientation of {1} in {0} [ROLL, PITCH, YAW]
% 'angularVel' contains angular velocity of {1} in {1}
% 'pos' contains position of {2} wrt {1} in {1}
%
% Author: Srivenkata Satya Prasad Maddipatla
% Date: 26th Oct, 2019
% Working of gyroscope needs to be confirmed


%% Error checking
% Check to see if data is an integer type - if so, it will cause problems
% when doing floating-point calculations and give odd results.
if isinteger(velArray)
    error(message('MATLAB:var:integerClass'));
end

% Check to see if the number of arguments is correct
if nargin ~= 4 && nargin ~= 5
    error('Invalid number of input arguments: Expecting minimum four arguments. Fifth argument should be deg or no entry');
end

% Check to see if the matrices are of same size
if ~isequal(size(angularPos),size(angularVel))
    error('Inconsistent dimensions on input arguments. angularPos must have same size as angularVel.');
end

% Check to see if the input arguments do not have NaN
if max(max(isnan(velArray))) || max(max((isnan(angularPos)))) || max(max((isnan(angularVel)))) || max((isnan(pos)))
    error('NaN data detected in input arguments. All inputs must be of numeric values.');
end


%% Data processing
if (1 == length(varargin))
    if strcmp('deg',varargin{1})
        rollArray   = (pi/180)*angularPos(:,1);
        pitchArray  = (pi/180)*angularPos(:,2);
        yawArray    = (pi/180)*angularPos(:,3);

        angularVel  = (pi/180)*angularVel;
    end
else
    rollArray   = angularPos(:,1);
    pitchArray  = angularPos(:,2);
    yawArray    = angularPos(:,3);
    
end

% cosine and sine of roll
croll   = cos(rollArray);
sroll   = sin(rollArray);
% cosine and sine of pitch
cpitch  = cos(pitchArray);
spitch  = sin(pitchArray);
% cosine and sine of yaw
cyaw    = cos(yawArray);
syaw    = sin(yawArray);

% rij forms the element in i th row and j th column
% of rotation matrix
r11 = cyaw .* cpitch;
r12 = (cyaw .* spitch .* sroll) - (syaw .* croll);
r13 = (cyaw .* spitch .* croll) + (syaw .* sroll);
r21 = syaw .* cpitch;
r22 = (syaw .* spitch .* sroll) + (cyaw .* croll);
r23 = (syaw .* spitch .* croll) - (cyaw .* sroll);
r31 = -spitch;
r32 = cpitch .* sroll;
r33 = cpitch .* croll;

% Create a position array equal to the size of the desired output array
posArray    = pos .* ones(size(angularVel));
% Cross product of angular velocity of {1} in {1} and position of {2} wrt {1} in {1}
crossProd   = cross(angularVel,posArray);
xCross      = crossProd(:,1);
yCross      = crossProd(:,2);
zCross      = crossProd(:,3);

% Transforming cross product from {1} to {0}
crossProd   = [r11 .* xCross + r12 .* yCross + r13 .* zCross,...
               r21 .* xCross + r22 .* yCross + r23 .* zCross,...
               r31 .* xCross + r32 .* yCross + r33 .* zCross];

if (3 == size(velArray,2))
    velTrans = velArray + crossProd;
else
    velTrans = velArray + crossProd(:,1:2);
end

end