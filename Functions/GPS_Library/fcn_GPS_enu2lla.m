function path_LLA = fcn_GPS_enu2lla(path_ENU, reference_LLA, varargin)
% fcn_GPS_enu2lla.m
% transforms a path(s) in ENU coordinate system to Geodetic coordinate
% system. This is written to test the GPS class.
%
% FORMAT:
%   path_LLA = fcn_GPS_enu2lla(path_ENU, reference_LLA, varargin)
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
% This function was written on 2021_01_14 by Satya Prasad
% Questions or comments? szm888@psu.edu

% Revision history:
%   2021_01_14:
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
    
    % INPUT: path_ENU
    % Check the size of inputs
    if 1 > size(path_ENU,1) || 3 ~= size(path_ENU,2)
        error('Input(path_ENU) must be a Nx3 vector.')
    end
    
    % Check the type and validity of inputs
    if ~isnumeric(path_ENU) || any(isnan(path_ENU),'all')
        error('Input(path_ENU) must be numeric data.')
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
path_LLA = NaN(size(path_ENU,1),3);

for i = 1:size(path_ENU,1)
    point_XYZ = fcn_GPS_enu2xyz(path_ENU(i,:), reference_LLA); % ENU to ECEF transformation
    path_LLA(i,:) = fcn_GPS_xyz2lla(point_XYZ); % ECEF to Geodetic transformation
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
    plot(path_ENU(:,1), path_ENU(:,2))
    grid on
    axis equal
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
