function diff_angles = fcn_Path_calcDiffAnglesBetweenPathSegments(Path,varargin)
% fcn_Path_calcDiffAnglesBetweenSegments
% Calculates the change in angles between segments. If there are N points
% in the path, there are N-1 segments and thus N-2 angles between segments.
%
% Note that this method uses the dot product and the cross product to avoid
% errors caused by angle rollover.
%
% FORMAT: 
%
%       diff_angles = fcn_Path_calcDiffAnglesBetweenPathSegments(Path,(fig_num))
%
% INPUTS:
%
%      Path: an N x 2 vector with [X Y] data in each row. N must be >= 3.
%
%      (OPTIONAL INPUTS)
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      diff_angles: tan (N-2) x 1 vector of the change in angles in radians
%
% EXAMPLES:
%      
%       See the script:
%       script_test_fcn_Path_calcDiffAnglesBetweenPathSegments.m for a full
%       test suite.
%
% This function was written on 2021_01_03 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%     2021_01_03
%     -- first writing of the code
%     2021_01_06
%     -- fixed typos in the comments
%     2021_01_07
%     -- fixed typos in the comments


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
    fcn_Path_checkInputsToFunctions(Path, 'paths');

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

%% Solve for the circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note that this method uses the dot product and the cross product to avoid
% errors caused by angle rollover. One could use the arctan method of
% calculating angles, which works in general but fails for paths that point
% straight to the left due to the cross-over point for the atan2
% calculation. For example, the following code does NOT always work:
%
%     % angles = atan2(path_average(2:end,2)-path_average(1:end-1,2),path_average(2:end,1)-path_average(1:end-1,1));
%     % angles = [angles; angles(end)];  % Pad the last point twice
%     % diff_angles2 = abs(diff(angles));

a_vector = Path(2:end-1,:)-Path(1:end-2,:);
b_vector = Path(3:end,:)-Path(2:end-1,:);
a_dot_b = sum(a_vector.*b_vector,2); % Do the dot product
a_cross_b = crossProduct(a_vector,b_vector);

a_mag = sum(a_vector.^2,2).^0.5;
b_mag = sum(b_vector.^2,2).^0.5;

% Calculate the change in angle using the dot product
%diff_angles = [acos(a_dot_b./(a_mag.*b_mag)); 0];
diff_angles = sign(a_cross_b).*[acos(a_dot_b./(a_mag.*b_mag))];


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
    clf;
    hold on;
    grid on;
    xlabel('Index');
    ylabel('Angle [deg]');

    % Plot the angle differences
    plot(diff_angles*180/pi,'k.-','Linewidth',3,'Markersize',25);

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end
end

%% Calculate cross products
function result = crossProduct(v,w)
result = v(:,1).*w(:,2)-v(:,2).*w(:,1);
end

