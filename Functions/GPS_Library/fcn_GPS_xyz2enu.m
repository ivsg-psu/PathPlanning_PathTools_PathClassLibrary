function point_ENU = fcn_GPS_xyz2enu(point_XYZ, reference_LLA)
% fcn_GPS_xyz2enu.m
% transforms a point in ECEF coordinate system to ENU coordinate system. 
% This is written to test the GPS class.
%
% FORMAT:
%   point_ENU = fcn_GPS_xyz2enu(point_XYZ, reference_LLA)
%
% INPUTS:
%   point_XYZ: a point as 1x3 vector in ECEF coordinate system
%   reference_LLA: a reference point as 1x3 vector in Geodetic coordinate
%   system
%
% OUTPUTS:
%   point_ENU: a point as 1x3 vector in ENU coordinate system
%
% EXAMPLES:
%   See the script: script_test_fcn_GPS_xyz2enu.m for a full test suite.
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
    
    fcn_GPS_checkInputsToFunctions(point_XYZ, 'point_ECEF')
    fcn_GPS_checkInputsToFunctions(reference_LLA, 'point_LLA')
end

%% Convert from ECEF to ENU coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining constants
factor_deg2rad = pi/180; % multiplying factor to convert degrees into radians

% Transformation process
% Read input data into separate variables
reference_latitude  = reference_LLA(1,1);
reference_longitude = reference_LLA(1,2);

% Rotate the diffxyz vector (often short) about Z-axis followed by X-axis
% to ENU frame
slat = sin(factor_deg2rad*reference_latitude); % sine of latitude
clat = cos(factor_deg2rad*reference_latitude); % cosine of latitude

slon = sin(factor_deg2rad*reference_longitude); % sine of longitude
clon = cos(factor_deg2rad*reference_longitude); % cosine of longitude

rotation_matrix = eye(3);
rotation_matrix(1,1) = -slon;
rotation_matrix(1,2) = clon;
rotation_matrix(2,1) = -slat * clon;
rotation_matrix(2,2) = -slat * slon;
rotation_matrix(2,3) = clat;
rotation_matrix(3,1) = clat * clon;
rotation_matrix(3,2) = clat * slon;
rotation_matrix(3,3) = slat;

reference_XYZ = fcn_GPS_lla2xyz(reference_LLA); % Transform reference_LLA into ECEF coordinates
diffXYZ = point_XYZ - reference_XYZ; % Shift the origin of ECEF to the reference

point_ENU = (rotation_matrix * diffXYZ')';

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
