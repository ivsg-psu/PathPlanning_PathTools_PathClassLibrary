function path_XYZ = fcn_GPS_enu2xyz(path_ENU, reference_LLA, varargin)
% fcn_GPS_enu2xyz.m
% transforms a path(s) in ENU coordinate system to ECEF coordinate system. 
% This is written to test the GPS class.
%
% FORMAT:
%   path_XYZ = fcn_GPS_enu2xyz(path_ENU, reference_LLA, fig_num)
%
% INPUTS:
%   path_ENU: a path(s) as 1x3 vector in ENU coordinate system
%   reference_LLA: a reference point as 1x3 vector in Geodetic coordinate
%   system
%   varargin: figure number for debugging or plotting
%
% OUTPUTS:
%   path_XYZ: a path(s) as Nx3 vector in ECEF coordinate system
%
% EXAMPLES:
%   See the script: script_test_fcn_GPS_enu2xyz.m for a full test suite.
%
% Dependencies:
%       fcn_GPS_checkInputsToFunctions
%       fcn_GPS_lla2xyz
%
% This function was written on 2021_01_14 by Satya Prasad
% Questions or comments? szm888@psu.edu

% Revision history:
%   2021_01_14:
%       - wrote the code
%   2021_01_25:
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
    if 2 > nargin || 3 < nargin
        error('Incorrect number of input arguments.')
    end
    
    fcn_GPS_checkInputsToFunctions(path_ENU, 'path_ENU')
    fcn_GPS_checkInputsToFunctions(reference_LLA, 'point_LLA')
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

diffXYZ = [ sum(rotation_matrix(1,:).*path_ENU, 2), ...
            sum(rotation_matrix(2,:).*path_ENU, 2), ...
            sum(rotation_matrix(3,:).*path_ENU, 2) ];

reference_XYZ = fcn_GPS_lla2xyz(reference_LLA); % Transform reference_LLA into ECEF coordinates
path_XYZ      = diffXYZ + reference_XYZ; % Shift the origin from reference to center of earth

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
    plot(path_ENU(:,1), path_ENU(:,2))
    grid on
    axis equal
    title('Path in ENU Coordinate System')
    xlabel('East (m)')
    ylabel('North (m)')
    
    % Tile 2
    nexttile
    plot(path_XYZ(:,1), path_XYZ(:,2))
    grid on
    axis equal
    title('Path in ECEF Coordinate System')
    xlabel('X-direction (m)')
    ylabel('Y-direction (m)')
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n', st(1).name, st(1).file);
end
end
