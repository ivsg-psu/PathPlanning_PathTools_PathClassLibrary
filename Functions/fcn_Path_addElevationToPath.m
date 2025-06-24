function elevated_path = fcn_Path_addElevationToPath(path, ...
                         reference_elevated_path, varargin)
% fcn_Path_addElevationToPath
% Adds elevation to the 2-dimensional 'path' based on the nearest neighbors
% in a 3-dimesional 'reference_elevated_path'
% 
% FORMAT:
%
%      elevated_path = fcn_Path_addElevationToPath(path, ...
%                      reference_elevated_path, (fig_num))
%
% INPUTS:
%
%      path: a Nx2 vector of [X Y] path points, where N is the number of
%      points on the path, N >= 2.
%      reference_elevated_path: a Nx3 vector of [X Y Z] path points, where 
%      N is the number of points on the path, N >= 2.
%
%      (optional inputs)
%
%      figure_number: figure number where results are plotted
%
% OUTPUTS:
%
%      elevated_path: a Nx3 vector of [X Y Z] path points, where N is the 
%      number of points on the path, N >= 2.
%
% EXAMPLES:
% 
% See the script: script_test_fcn_Path_addElevationToPath for a full test 
% suite.
%
% DEPENDENCIES:
%
%     fcn_Path_snapPointToPathViaVectors
%     fcn_DebugTools_checkInputsToFunctions
%
% This function was written on 2021_03_06 by Satya Prasd
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2021_03_06 - S. Brennan
% -- wrote the code
% 2022_01_07 - S. Brennan
% -- updated header to fix input definitions
% 2022_08_20 - S. Brennan
% -- allow empty figure argument to avoid plotting
% 2024_03_14 - S. Brennan
% -- defaulted plotting view to 3D
% 2025_06_23 - S. Brennan
% -- Updated debugging and input checks

% TO-DO
% 2025_06_23 - S. Brennan
% -- Remove the for loop after fcn_Path_findSensorHitOnWall is implemented
% -- Transition code to use fcn_Path_findSensorHitOnWall instead

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==3 && isequal(varargin{end},-1))
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_PATHCLASS_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_PATHCLASS_FLAG_CHECK_INPUTS");
    MATLABFLAG_PATHCLASS_FLAG_DO_DEBUG = getenv("MATLABFLAG_PATHCLASS_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_PATHCLASS_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_PATHCLASS_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_PATHCLASS_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_PATHCLASS_FLAG_CHECK_INPUTS);
    end
end

% flag_do_debug = 1;

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 999978; %#ok<NASGU>
else
    debug_fig_num = []; %#ok<NASGU>
end

%% check input arguments?
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
if 0==flag_max_speed
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(2,3);

        % Check the data inputs
        fcn_DebugTools_checkInputsToFunctions(path, 'path');
        fcn_DebugTools_checkInputsToFunctions(reference_elevated_path, 'elevated_path');
    end
end

% Does user want to show the plots?
flag_do_plots = 0; % Default is to NOT show plots
if (0==flag_max_speed) && (3 == nargin) 
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        fig_num = temp;
        figure(fig_num);
        flag_do_plots = 1;
    end
else
    if flag_do_debug
        fig = figure;  
        fig_num = fig.Number;
        flag_do_plots = 1;
    end
end

%% Add elevation to the path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Npoints = length(path(:,1)); % number of points in the path
% Initialize vectors to store
% 1) indices of nearest neighbor that is behind the points of path
% 2) indices of nearest neighbor that is ahead of the points of path
% 3) percentage of travel of point on the path
vector_first_path_point_index   = NaN(Npoints,1);
vector_second_path_point_index  = NaN(Npoints,1);
vector_percent_along_length     = NaN(Npoints,1);
closest_point_on_reference_path = NaN(Npoints,3);

Ndimensions = length(path(1,:));

for i = 1:Npoints
    %     [closest_point_on_reference_path(i,:), ~, vector_first_path_point_index(i), ...
    %         vector_second_path_point_index(i), vector_percent_along_length(i)] ...
    %         = fcn_Path_snapPointToPathViaVectors(path(i,:), reference_elevated_path(:,[1,2]));
    
    if 2 == Ndimensions
        query_point = [path(i,:) 0];
    elseif 3 == Ndimensions
        query_point = path(i,:);
    else
        error('Number of columns in the path must be 2 or 3, e.g. for 2D points or 3D points');
    end

    [closest_point_on_reference_path(i,:), ~, vector_first_path_point_index(i), ...
        vector_second_path_point_index(i), vector_percent_along_length(i)] ...
        = fcn_Path_snapPointToPathViaVectors(query_point, reference_elevated_path);
end
% % estimate the elevation as an average
% elevation = (1-vector_percent_along_length).*reference_elevated_path(vector_first_path_point_index,3)+...
%              vector_percent_along_length.*reference_elevated_path(vector_second_path_point_index,3);
% elevated_path = [path, elevation];

elevated_path = [path, closest_point_on_reference_path(:,3)];

%% Plot the results (for debugging)?
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
    grid on
    hold on
    % Plot the path
    plot3(reference_elevated_path(:,1),reference_elevated_path(:,2),...
          reference_elevated_path(:,3),'r.-','Linewidth',2,'Markersize',20);
    plot3(path(:,1),path(:,2),zeros(Npoints,1),'b.-','Linewidth',2,'Markersize',20);
    plot3(elevated_path(:,1),elevated_path(:,2),elevated_path(:,3),'g.-',...
        'Linewidth',2,'Markersize',20);
    plot3(...
        closest_point_on_reference_path(:,1),...
        closest_point_on_reference_path(:,2),...
        closest_point_on_reference_path(:,3),'k.',...
        'Linewidth',2,'Markersize',20);
    
    
    for ith_point = 1:length(path)
        % Plot dotted lines from the path to elevated path
        plot3(...
            [path(ith_point,1)   elevated_path(ith_point,1)],...
            [path(ith_point,2)   elevated_path(ith_point,2)],...
            [0                   elevated_path(ith_point,3)],...
            'k--','Linewidth',2);
        
        % Plot dotted lines from the reference to elevated path
        plot3(...
            [closest_point_on_reference_path(ith_point,1)   elevated_path(ith_point,1)],...
            [closest_point_on_reference_path(ith_point,2)   elevated_path(ith_point,2)],...
            [closest_point_on_reference_path(ith_point,3)   elevated_path(ith_point,3)],...
            'k--','Linewidth',2);
    end
    
    legend('Reference Path', 'Input Path', 'Elevated Path')
    xlabel('xEast')
    ylabel('yNorth')
    zlabel('zUp')

    view(3);
end % Ends the flag_do_debug if statement

end % Ends the main function


%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§

