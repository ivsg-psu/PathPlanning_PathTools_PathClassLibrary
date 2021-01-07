function pointsTrans = fcn_positionTransform(pointArray, angularPos, pos, varargin)
% Finds position of hemisphere-{2} in reference frame-{0} 
% given pose of novatel-{1} in reference frame-{0} and
% position of hemisphere-{2} wrt novatel-{1} in novatel frame-{1}.
%
% pos = 0.0254*[21, 0, 0] Hemisphere WRT Novatel-GPS in Novatel-GPS
% pos = 0.0254*[8.625, 0, 0] Hemisphere WRT Novatel-IMU in Novatel-IMU
% pos = 0.0254*[18, -1, 0] Hemisphere WRT ADIS-IMU in ADIS-IMU
%
% 'pointArray' contains position of {1} in {0}
% 'angularPos' contains orientation of {1} in {0} [ROLL, PITCH, YAW]
% 'pos' contains position of {2} wrt {1} in {1}
%
% Author: Srivenkata Satya Prasad Maddipatla
% Date: 26th Oct, 2019


%% Error checking
% Check to see if data is an integer type - if so, it will cause problems
% when doing floating-point calculations and give odd results.
if isinteger(pointArray)
    error(message('MATLAB:var:integerClass'));
end

% Check to see if the number of arguments is correct
if nargin ~= 3 && nargin ~= 4
    error('Invalid number of input arguments: Expecting minimum three arguments. Fourth argument must be deg or no entry');
end

% Check to see if the matrices are of same size
if ~isequal(size(pointArray),size(angularPos))
    error('Inconsistent dimensions on input arguments. pointArray must have same size as angularPos.');
end

% Check to see if the input arguments do not have NaN
if max(max(isnan(pointArray))) || max(max((isnan(angularPos)))) || max((isnan(pos)))
    error('NaN data detected in input arguments. All inputs must be of numeric values.');
end

%% Data Processing
if (1 == length(varargin))
    if strcmp('deg',varargin{1})
        rollArray   = (pi/180)*angularPos(:,1);
        pitchArray  = (pi/180)*angularPos(:,2);
        yawArray    = (pi/180)*angularPos(:,3);
    else
        error('Fouth argument must be deg.');
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
% of transformation matrix
r11 = cyaw .* cpitch;
r12 = (cyaw .* spitch .* sroll) - (syaw .* croll);
r13 = (cyaw .* spitch .* croll) + (syaw .* sroll);
r21 = syaw .* cpitch;
r22 = (syaw .* spitch .* sroll) + (cyaw .* croll);
r23 = (syaw .* spitch .* croll) - (cyaw .* sroll);
r31 = -spitch;
r32 = cpitch .* sroll;
r33 = cpitch .* croll;

% find position of {2} given the position of {1}, orientation of {1}
% and position of {2} from {1} in {1}
pointsTrans = pointArray + [r11*pos(1)+r12*pos(2)+r13*pos(3),...
                            r21*pos(1)+r22*pos(2)+r23*pos(3),...
                            r31*pos(1)+r32*pos(2)+r33*pos(3)];

end