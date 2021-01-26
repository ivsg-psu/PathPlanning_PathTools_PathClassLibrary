function point_XYZ = fcn_GPS_lla2xyz(point_LLA)
% fcn_GPS_lla2xyz.m
% transforms a point in Geodetic coordinate system to ECEF coordinate
% system. This is written to test the GPS class.
%
% FORMAT:
%   point_XYZ = fcn_GPS_lla2xyz(point_LLA)
%
% INPUTS:
%   point_LLA: a point as 1x3 vector in Geodetic coordinate system
%
% OUTPUTS:
%   point_XYZ: a point as 1x3 vector in ECEF coordinate system
%
% EXAMPLES:
%   See the script: script_test_fcn_GPS_lla2xyz.m for a full test suite.
%
% This function was written on 2021_01_14 by Satya Prasad
% Questions or comments? szm888@psu.edu

% Revision history:
%   2021_01_14:
%       - wrote the code
% 2021_01_25:
%       - Added function to check inputs

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
    if 1 ~= nargin
        error('Incorrect number of input arguments.')
    end
    
    fcn_GPS_checkInputsToFunctions(point_LLA, 'point_LLA')
end

%% Convert from Geodetic to ECEF coordinates
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
GPS.semi_major_axis = 6378137; % semi-major axis of the earth [meters]
GPS.flattening = 1/298.257223563; % flattening of the earth
GPS.first_eccentricity_squared = (2 - GPS.flattening) * GPS.flattening; % square of the first eccentricity of the ellipsoid

% Read input data into separate variables
latitude  = point_LLA(1,1);
longitude = point_LLA(1,2);
altitude  = point_LLA(1,3);

% Transformation process
slat = sin(factor_deg2rad*latitude); % sine of latitude
clat = cos(factor_deg2rad*latitude); % cosine of latitude

slon = sin(factor_deg2rad*longitude); % sine of longitude
clon = cos(factor_deg2rad*longitude); % cosine of longitude

primal_vertical_radius_of_curvature = GPS.semi_major_axis / sqrt(1 - GPS.first_eccentricity_squared * slat * slat); % primal vertical radius of curvature [meters]

% first column is X, second column is Y, and third column is Z in ECEF
point_XYZ = [ (primal_vertical_radius_of_curvature + altitude) * clat * clon,...
              (primal_vertical_radius_of_curvature + altitude) * clat * slon,...
              ((1 - GPS.first_eccentricity_squared) * primal_vertical_radius_of_curvature + altitude) * slat ];

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
