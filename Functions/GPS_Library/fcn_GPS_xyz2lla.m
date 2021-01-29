function path_LLA = fcn_GPS_xyz2lla(path_XYZ, varargin)
% fcn_GPS_xyz2lla.m
% transforms a path(s) in ECEF coordinate system to Geodetic coordinate
% system. This is written to test the GPS class.
%
% FORMAT:
%   path_LLA = fcn_GPS_xyz2lla(path_XYZ, fig_num)
%
% INPUTS:
%   path_XYZ: a path(s) as Nx3 vector in ECEF coordinate system
%
% OUTPUTS:
%   path_LLA: a path(s) as Nx3 vector in Geodetic coordinate system
%
% EXAMPLES:
%   See the script: script_test_fcn_GPS_xyz2lla.m for a full test suite.
%
% This function was written on 2021_01_14 by Satya Prasad
% Questions or comments? szm888@psu.edu
% Reference: http://read.pudn.com/downloads559/sourcecode/others/2303011/gps_matlab/Wgsxyz2lla.m__.htm

% Revision history:
%   2021_01_14:
%       - wrote the code
%   2021_01_25:
%       - Added function to check inputs
%   2021_01_28:
%       - Added a for loop so that the function will work for multiple
%       points as well

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
if flag_check_inputs == 1
    % Are there the right number of inputs?
    if  1 > nargin || 2 < nargin
        error('Incorrect number of input arguments.')
    end
    
    fcn_GPS_checkInputsToFunctions(path_XYZ, 'path_ECEF')
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

%% Convert from ECEF to Geodetic coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define some parameters
factor_rad2deg = 180/pi; % multiplying factor to convert radians into degrees
GPS.semi_major_axis = 6378137; % semi-major axis of the earth [meters]
GPS.flattening = 1/298.257223563; % flattening of the earth
GPS.first_eccentricity_squared = (2 - GPS.flattening) * GPS.flattening; % square of the first eccentricity of the ellipsoid

path_LLA = NaN(size(path_XYZ, 1), 3); % intialize variable to store the path in LLA

% loop over all the points
for i = 1:size(path_XYZ, 1)
    % Longitude estimation
    if (0.0 == path_XYZ(i,1)) && (0.0 == path_XYZ(i,2))
        longitude = 0.0;
    else
        longitude = factor_rad2deg * atan2(path_XYZ(i,2), path_XYZ(i,1));
    end

    % Latitude and altitude estimation
    if (0.0 == path_XYZ(i,1)) && (0.0 == path_XYZ(i,2)) && (0.0 == path_XYZ(i,3))
        error('WGS XYZ is located at the center of the earth')
    else
        % Make initial estimates of latitude and altitude based on spherical 
        % earth assumption
        rho_squared = path_XYZ(i,1)^2 + path_XYZ(i,2)^2;
        rho = sqrt(rho_squared);
        est_lat = atan2(path_XYZ(i,3), rho);
        est_alt = sqrt(rho_squared + path_XYZ(i,3)^2) - GPS.semi_major_axis;
        % initialize errors/residuals on rho and z
        error_rho = 1.0;
        error_z   = 1.0;

        % Newton's method iteration on est_lat and est_alt makes the
        % residuals on rho and z progressively smaller.
        % Loop is implemented as a 'while' instead of a 'do' to simplify
        % porting to MATLAB
        while (1e-6 < abs(error_rho)) || (1e-6 < abs(error_z))
            slat = sin(est_lat); % sine of latitude
            clat = cos(est_lat); % cosine of latitude

            q = 1 - GPS.first_eccentricity_squared * slat * slat;
            primal_vertical_radius_of_curvature = GPS.semi_major_axis / sqrt(q); % primal vertical radius of curvature [meters]
            drdl = primal_vertical_radius_of_curvature * GPS.first_eccentricity_squared * slat * clat / q; % d(primal_vertical_radius_of_curvature)/d(est_lat)

            % error in rho and z (Z-coordinates in ECEF)
            error_rho = (primal_vertical_radius_of_curvature + est_alt) * clat - rho;
            error_z   = ((1 - GPS.first_eccentricity_squared) * ...
                primal_vertical_radius_of_curvature + est_alt) * slat - path_XYZ(i,3);

            %   --   Find Jacobian               -- --      --
            %   |  d(error_rho)/d(est_lat)  d(error_rho)/d(est_alt) | |  a  b  |
            %   |                                                   |=|        |
            %   |  d(error_z)/d(est_lat)    d(error_z)/d(est_alt)   | |  c  d  |
            %   --                               -- --      --
            jacobian_a = drdl * clat - (primal_vertical_radius_of_curvature + est_alt) * slat;
            jacobian_b = clat;
            jacobian_c = (1 - GPS.first_eccentricity_squared) * drdl * slat + ...
                         ((1 - GPS.first_eccentricity_squared) * primal_vertical_radius_of_curvature + est_alt) * clat;
            jacobian_d = slat;

            determinant_inverse = 1.0 / (jacobian_a * jacobian_d - jacobian_b * jacobian_c);
            est_lat = est_lat - determinant_inverse * (jacobian_d*error_rho - jacobian_b*error_z);
            est_alt = est_alt - determinant_inverse * (-jacobian_c*error_rho + jacobian_a*error_z);
        end

        latitude = factor_rad2deg * est_lat;
        altitude = est_alt;
    end

    path_LLA(i,:) = [latitude, longitude, altitude];
end

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
    plot(path_XYZ(:,1), path_XYZ(:,2))
    grid on
    axis equal
    title('Path in ECEF Coordinate System')
    xlabel('X-direction (m)')
    ylabel('Y-direction (m)')
    
    % Tile 2
    nexttile
    geoplot(path_LLA(:,1), path_LLA(:,2))
    geobasemap colorterrain
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n', st(1).name, st(1).file);
end
end
