function clean_path = fcn_Path_cleanPathFromForwardBackwardJogs...
    (path_with_jogs,varargin)
% Finds and removes situations where the path is jumping forward and
% backward. This is detected by finding situations where the angle between
% segments is more than a threshold (currently pi/4), and then taking these
% segments, and the one before and after, and removing them. It then
% re-scans the path (up to 3 times) to again check for these situations.
%
% FORMAT: 
%
%     clean_path = fcn_Path_cleanPathFromForwardBackwardJogs...
%     (path_with_jogs, (fig_num));
%
% INPUTS:
%
%      path_with_jogs: a paths type consisting of (N x 2) array with N>=3
%
%      (OPTIONAL INPUTS)
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%      clean_path: the resulting path
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_Path_calcDiffAnglesBetweenPathSegments
%
% EXAMPLES:
%      
%     See the script: 
%     script_test_fcn_Path_cleanPathFromForwardBackwardJogs
%     for a full test suite.
%
% This function was written on 2021_01_09 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%     
%     2021_01_09:
%     -- wrote the code originally 

flag_do_debug = 0; % Flag to show the results for debugging
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

    % Check the data input
    fcn_DebugTools_checkInputsToFunctions(path_with_jogs, 'paths');
        
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

%% Main code starts here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iteration_count = 1;
flag_average_is_good = 0;

while (0==flag_average_is_good)  && (iteration_count<=3)
    % Calculate angle changes between points
    diff_angles = fcn_Path_calcDiffAnglesBetweenPathSegments(path_with_jogs);
    
    % Find outliers
    outliers = find(abs(diff_angles)>pi/4);
    
    % Are there any back/forth jogs?
    if ~isempty(outliers)
        % ID adjacent points
        outliers = unique([outliers;outliers+1;outliers-1]);
        outliers = min(outliers,length(path_with_jogs)-1);
        outliers = max(1,outliers);
        
        % Create a set of indices we will save
        indices = (1:length(path_with_jogs(:,1)))';
        indices(outliers+1) = 0;
        
        % Save the clean path
        clean_path = path_with_jogs(indices~=0,:);        
        
        % For debugging
        if flag_do_debug
            figure(8888);
            clf;
            hold on;
            grid on;
            grid minor;
            
            x_indices = 1:length(diff_angles);
            plot(x_indices,diff_angles,'k-');
            plot(x_indices(outliers), diff_angles(outliers),'ro');
        end
        
    else % No back/forth jogs left!
        flag_average_is_good = 1;
        clean_path = path_with_jogs;
    end
    
    % Show results for debugging?
    if flag_do_debug
        figure(777777);
        clf;
        hold on;
        grid on;
        grid minor;
        axis equal;
        
        plot(path_with_jogs(:,1),path_with_jogs(:,2),'k.-');
        plot(path_with_jogs(outliers+1,1),path_with_jogs(outliers+1,2),'ro');
        plot(clean_path(:,1),clean_path(:,2),'b-');
    end
    
    % Increment the iteration count
    iteration_count = iteration_count + 1;
    
    % Reset the path average for the next round
    path_with_jogs = clean_path;
end


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
    
    figure(fig_num);
    clf;
    hold on;
    grid on;
    grid minor;
    axis equal;
    
    plot(path_with_jogs(:,1),path_with_jogs(:,2),'k.-');
    plot(path_with_jogs(outliers+1,1),path_with_jogs(outliers+1,2),'ro');
    plot(clean_path(:,1),clean_path(:,2),'b-');
    
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end
end % End of function



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


