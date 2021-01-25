function point_XYZ = fcn_GPS_enu2xyz(point_ENU, reference_LLA)
% fcn_GPS_enu2xyz.m
% transforms a point in ENU coordinate system to ECEF coordinate system. 
% This is written to test the GPS class.
%
% FORMAT:
%   point_XYZ = fcn_GPS_enu2xyz(point_ENU, reference_LLA)
%
% INPUTS:
%   point_ENU: a point as 1x3 vector in ENU coordinate system
%   reference_LLA: a reference point as 1x3 vector in Geodetic coordinate
%   system
%
% OUTPUTS:
%   point_XYZ: a point as 1x3 vector in ECEF coordinate system
%
% EXAMPLES:
%   See the script: script_test_fcn_GPS_enu2xyz.m for a full test suite.
%
% This function was written on 2021_01_14 by Satya Prasad
% Questions or comments? szm888@psu.edu

% Revision history:
%   2021_01_14:
%       - wrote the code

flag_do_debug = 0; % Flag to plot the results for debugging
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'Starting function: %s, in file: %s\n', st(1).name, st(1).file);
end

%% check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _       
%  |_   _|                 | |      
%    | |  _ __  _ __  _   _| |_ ___ 
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |                  
%              |_| 
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_check_inputs
    % Are there the right number of inputs?
    if 2 ~= nargin
        error('Incorrect number of input arguments.')
    end
    
    % INTPUT: point_ENU
    % Check the size of inputs
    if 1 ~= size(point_ENU,1) || 3 ~= size(point_ENU,2)
        error('Input(point_ENU) must be a 1x3 vector.')
    end
    
    % Check the type and validity of inputs
    if ~isnumeric(point_ENU) || any(isnan(point_ENU))
        error('Input(point_ENU) must be numeric data.')
    end
    
    % INTPUT: reference_LLA
    % Check the size of inputs
    if 1 ~= size(reference_LLA,1) || 3 ~= size(reference_LLA,2)
        error('Input(reference_LLA) must be a 1x3 vector.')
    end

    % Check the type and validity of inputs
    if ~isnumeric(reference_LLA) || any(isnan(reference_LLA))
        error('reference_LLA must be numeric data.')
    end

    % Check the domain of inputs (latitude and longitude)
    if ((-90.0 > reference_LLA(1,1)) || (90.0 < reference_LLA(1,1)) || ...
            (-180.0 > reference_LLA(1,2)) || (180.0 < reference_LLA(1,2)))
        error('WGS lat or WGS lon are out of range');
    end
end

%% Convert from ENU to ECEF coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining parameters
factor_deg2rad = pi/180; % multiplying factor to convert degrees into radians

% Transformation process
% Read input data into separate variables
reference_latitude  = reference_LLA(1,1);
reference_longitude = reference_LLA(1,2);

% Rotate the ENU vector to ECEF coordinate system
slat = sin(factor_deg2rad*reference_latitude); % sine of latitude
clat = cos(factor_deg2rad*reference_latitude); % cosine of latitude

slon = sin(factor_deg2rad*reference_longitude); % sine of longitude
clon = cos(factor_deg2rad*reference_longitude); % cosine of longitude

rotation_matrix = eye(3);
rotation_matrix(1,1) = -slon;
rotation_matrix(1,2) = -slat * clon;
rotation_matrix(1,3) = clat * clon;
rotation_matrix(2,1) = clon;
rotation_matrix(2,2) = -slat * slon;
rotation_matrix(2,3) = clat * slon;
rotation_matrix(3,2) = clat;
rotation_matrix(3,3) = slat;

diffXYZ = (rotation_matrix * point_ENU')';

reference_XYZ = fcn_GPS_lla2xyz(reference_LLA); % Transform reference_LLA into ECEF coordinates
point_XYZ = diffXYZ + reference_XYZ; % Shift the origin from reference to center of earth

%% Any debugging?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _                 
%  |  __ \     | |                
%  | |  | | ___| |__  _   _  __ _ 
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n', st(1).name, st(1).file);
end
end
