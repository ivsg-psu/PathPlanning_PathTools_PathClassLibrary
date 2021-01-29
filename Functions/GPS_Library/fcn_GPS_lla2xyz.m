function path_XYZ = fcn_GPS_lla2xyz(path_LLA, varargin)
% fcn_GPS_lla2xyz.m
% transforms a path(s) in Geodetic coordinate system to ECEF coordinate
% system. This is written to test the GPS class.
%
% FORMAT:
%   path_XYZ = fcn_GPS_lla2xyz(path_LLA, fig_num)
%
% INPUTS:
%   path_LLA: a path(s) as Nx3 vector in Geodetic coordinate system
%   varargin: figure number for debugging or plotting
%
% OUTPUTS:
%   path_XYZ: a path(s) as Nx3 vector in ECEF coordinate system
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
%   2021_01_28:
%       - Vectorized the function

flag_do_debug     = 0; % Flag to plot the results for debugging
flag_do_plots     = 0; % Flag to plot the final results
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
    if 1 > nargin || 2 < nargin
        error('Incorrect number of input arguments.')
    end
    
    fcn_GPS_checkInputsToFunctions(path_LLA, 'path_LLA')
end

%% Check for variable argument inputs (varargin)

% Does user want to show the plots?
if 2 == nargin
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
factor_deg2rad      = pi/180; % multiplying factor to convert degrees into radians
GPS.semi_major_axis = 6378137; % semi-major axis of the earth [meters]
GPS.flattening      = 1/298.257223563; % flattening of the earth
GPS.first_eccentricity_squared = (2-GPS.flattening) * GPS.flattening; % square of the first eccentricity of the ellipsoid

% Read input data into separate variables
latitude  = path_LLA(:,1);
longitude = path_LLA(:,2);
altitude  = path_LLA(:,3);

% Transformation process begins from here
slat = sin(factor_deg2rad*latitude); % sine of latitude
clat = cos(factor_deg2rad*latitude); % cosine of latitude

slon = sin(factor_deg2rad*longitude); % sine of longitude
clon = cos(factor_deg2rad*longitude); % cosine of longitude

primal_vertical_radius_of_curvature = GPS.semi_major_axis ./ sqrt(1 - GPS.first_eccentricity_squared * slat.^2); % primal vertical radius of curvature [meters]

% first column is X, second column is Y, and third column is Z in ECEF
path_XYZ = [ (primal_vertical_radius_of_curvature+altitude) .* clat .* clon,...
             (primal_vertical_radius_of_curvature+altitude) .* clat .* slon,...
             ((1-GPS.first_eccentricity_squared) * primal_vertical_radius_of_curvature + altitude) .* slat ];

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
    plot(path_XYZ(:,1), path_XYZ(:,2))
    grid on
    axis equal
    xlabel('X-direction (m)')
    ylabel('Y-direction (m)')
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n', st(1).name, st(1).file);
end
end
