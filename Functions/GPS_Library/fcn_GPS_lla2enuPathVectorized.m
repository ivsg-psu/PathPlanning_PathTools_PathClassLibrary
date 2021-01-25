function path_ENU = fcn_GPS_lla2enuPathVectorized(path_LLA, reference_LLA, varargin)
% fcn_GPS_lla2xyzPathVectorized.m
% transforms a point in Geodetic coordinate system to ENU coordinate
% system. This is written to test the GPS class.
%
% FORMAT:
%   path_ENU = fcn_GPS_lla2enuPathVectorized(path_LLA, reference_LLA, varargin)
%
% INPUTS:
%   path_LLA: a point as Nx3 vector in Geodetic coordinate system
%
% OUTPUTS:
%   path_ENU: a point as Nx3 vector in ENU coordinate system
%
% EXAMPLES:
%   See the script: script_test_fcn_GPS_lla2enuPath.m for a full test suite.
%
% This function was written on 2021_01_25 by Satya Prasad
% Questions or comments? szm888@psu.edu

% Revision history:
%   2021_01_25:
%       - wrote the code

flag_do_debug = 0; % Flag to plot the results for debugging
flag_do_plots = 0; % Flag to plot the final results
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
    if 2 > nargin || 3 < nargin
        error('Incorrect number of input arguments')
    end
    
    % INPUT: path_LLA
    % Check the size of inputs
    if 1 > size(path_LLA,1) || 3 ~= size(path_LLA,2)
        error('Input(path_LLA) must be a Nx3 vector.')
    end
    
    % Check the type and validity of inputs
    if ~isnumeric(path_LLA) || any(isnan(path_LLA),'all')
        error('Input(path_LLA) must be numeric data.')
    end
    
    % Check the domain of inputs (latitude and longitude)
    if (any(-90.0 > path_LLA(:,1)) || any(90.0 < path_LLA(:,1)) || ...
            any(-180.0 > path_LLA(:,2)) || any(180.0 < path_LLA(:,2)))
        error('WGS lat or WGS lon are out of range');
    end
    
    % INPUT: reference_LLA
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

%% Check for variable argument inputs (varargin)

% Does user want to show the plots?
if 3 == nargin
    fig_num = varargin{1};
    figure(fig_num);
    flag_do_plots = 1;
else
    if flag_do_debug
        fig = figure;
        fig_num = fig.Number;
        flag_do_plots = 1;
    end
end

%% Convert from Geodetic to ENU coordinates
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
latitude  = [reference_LLA(1,1); path_LLA(:,1)];
longitude = [reference_LLA(1,2); path_LLA(:,2)];
altitude  = [reference_LLA(1,3); path_LLA(:,3)];

% Transformation process
slat = sin(factor_deg2rad*latitude); % sine of latitude
clat = cos(factor_deg2rad*latitude); % cosine of latitude

slon = sin(factor_deg2rad*longitude); % sine of longitude
clon = cos(factor_deg2rad*longitude); % cosine of longitude

primal_vertical_radius_of_curvature = GPS.semi_major_axis ./ sqrt(1 - GPS.first_eccentricity_squared * slat.^2); % primal vertical radius of curvature [meters]

% ECEF coordinate system
path_XYZ = [ (primal_vertical_radius_of_curvature + altitude) .* clat .* clon, ...
             (primal_vertical_radius_of_curvature + altitude) .* clat .* slon, ...
             ((1 - GPS.first_eccentricity_squared) * primal_vertical_radius_of_curvature + altitude) .* slat ];

% First element corresponds to reference_LLA
reference_XYZ = path_XYZ(1,:); % reference point for ENU coordinate system
diffxyz = path_XYZ(2:end,:) - reference_XYZ; % Shift the origin of ECEF to the reference

rotation_matrix = eye(3);
rotation_matrix(1,1) = -slon(1);
rotation_matrix(1,2) = clon(1);
rotation_matrix(2,1) = -slat(1) * clon(1);
rotation_matrix(2,2) = -slat(1) * slon(1);
rotation_matrix(2,3) = clat(1);
rotation_matrix(3,1) = clat(1) * clon(1);
rotation_matrix(3,2) = clat(1) * slon(1);
rotation_matrix(3,3) = slat(1);

path_ENU = [ sum(rotation_matrix(1,:).*diffxyz, 2), ...
             sum(rotation_matrix(2,:).*diffxyz, 2), ...
             sum(rotation_matrix(3,:).*diffxyz, 2) ];

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

if flag_do_plots
    figure(fig_num)
    clf
    tiledlayout(2,1)
    
    % Tile 1
    nexttile
    geoplot(path_LLA(:,1), path_LLA(:,2))
    geobasemap colorterrain
    
    % Tile 2
    nexttile
    plot(path_ENU(:,1), path_ENU(:,2))
    grid on
    axis equal
    xlabel('East (m)')
    ylabel('North (m)')
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n', st(1).name, st(1).file);
end
end