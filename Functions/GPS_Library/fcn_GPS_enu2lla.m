function path_LLA = fcn_GPS_enu2lla(path_ENU, reference_LLA, varargin)
% fcn_GPS_enu2llaPath.m
% transforms a path(s) in ENU coordinate system to Geodetic coordinate
% system. This is written to test the GPS class.
%
% FORMAT:
%   path_LLA = fcn_GPS_enu2llaPath(path_ENU, reference_LLA, varargin)
%
% INPUTS:
%   path_ENU: a path(s) as Nx3 vector in ENU coordinate system
%   reference_LLA: a reference point as 1x3 vector in Geodetic coordinate
%   system
%
% OUTPUTS:
%   path_LLA: a path(s) as Nx3 vector in Geodetic coordinate system
%
% EXAMPLES:
%   See the script: script_test_fcn_GPS_enu2lla.m for a full test suite.
%
% Dependencies:
%       fcn_GPS_checkInputsToFunctions
%       fcn_GPS_enu2xyz
%       fcn_GPS_xyz2lla
%
% This function was written on 2021_01_14 by Satya Prasad
% Questions or comments? szm888@psu.edu

% Revision history:
%   2021_01_14:
%       - wrote the code
% 2021_01_25:
%       - Added function to check inputs
% 2021_01_28:
%       - For loop is removed as it is moved to fcn_GPS_xyz2lla.M

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

%% Convert from ENU to Geodetic coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
point_XYZ = fcn_GPS_enu2xyz(path_ENU, reference_LLA); % ENU to ECEF transformation
path_LLA  = fcn_GPS_xyz2lla(point_XYZ); % ECEF to Geodetic transformation

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
    geoplot(path_LLA(:,1), path_LLA(:,2))
    geobasemap colorterrain
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n', st(1).name, st(1).file);
end
end
