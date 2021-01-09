function yaw_angles_in_radians = fcn_Path_calcYawFromPathSegments(Path,varargin)
% fcn_Path_calcYawFromPathSegments
% Calculates the yaw angles of the path segments. If there are N points
% in the path, there are N-1 segments and thus N-1 angles. 
%
% The method used for calculating the path angles is to use the first
% segment of the path to calculate the angle of this segment, using the
% arctangent function (atan2). For any subsequent segments, if any, the
% change in angles is calculated using the dot product and cross product,
% and these angle changes are then cumulatively added to determine yaw
% angles thereafter. This method with the dot product and the cross product
% allows us to avoid errors caused by angle rollover that occurs with the
% arctangent functions. Thus, a path that circles clockwise two times will
% have a yaw angle that starts somewhere between -pi and pi, and decreases
% by abou 4 pi, or two revolutions. Note that the angle convention for yaw
% follows Cartesian definitions, namely yaw is positive when angles
% increase counter-clockwise, and yaw angles start at 0 at the x-axis.
%
% FORMAT: 
%
%       yaw_angles_in_radians = ...
%       fcn_Path_calcYawFromPathSegments(Path,(fig_num))
%
% INPUTS:
%
%      Path: an N x 2 vector with [X Y] data in each row. N must be >= 2.
%
%      (OPTIONAL INPUTS)
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      yaw_angles_in_radians: the N-1 yaw angles in radians
%
% DEPENDENCIES
%
%      fcn_Path_calcDiffAnglesBetweenPathSegments
%
% EXAMPLES:
%      
%       See the script:
%       script_test_fcn_Path_calcYawFromPathSegments.m for a full
%       test suite.
%
% This function was written on 2021_01_06 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%     2021_01_06
%       -- wrote the code
%     2021_01_06
%       -- minor comment clean-ups


flag_do_debug = 0; % Flag to plot the results for debugging
flag_do_plots = 0; % Flag to plot the final results
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
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
    if nargin < 1 || nargin > 2
        error('Incorrect number of input arguments')
    end
    
    % Check the Path variables        
    fcn_Path_checkInputsToFunctions(Path, 'path');

end

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

%% Main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

starting_angle = atan2(Path(2,2)-Path(1,2),Path(2,1)-Path(1,1));

Nangles = length(Path(:,1)) - 1;
if Nangles>1
    diff_angles = fcn_Path_calcDiffAnglesBetweenPathSegments(Path);
    cumulative_angles = [0; cumsum(diff_angles)];
else
    cumulative_angles = 0;
end
yaw_angles_in_radians = cumulative_angles + starting_angle;


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
    figure(fig_num);
    
    hold on;
    grid on;
    xlabel('Index');
    ylabel('Yaw Angle [deg]');

    % Plot the angle differences
    plot(yaw_angles_in_radians*180/pi,'.-','Linewidth',3,'Markersize',25);

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end
end



