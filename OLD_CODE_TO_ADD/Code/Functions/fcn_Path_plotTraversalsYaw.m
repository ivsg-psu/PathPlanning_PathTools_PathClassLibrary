function h = fcn_Path_plotTraversalsYaw(data,varargin)
% fcn_Path_plotTraversalsYaw
% Plots the Yaw angles of all traversals existing in a data structure
%
% FORMAT: 
%
%       h = fcn_Path_plotTraversalsYaw(data,varargin)
%
% INPUTS:
%
%      data: a structure containing subfields of station and Yaw in the
%      following form
%           data.traversal{i_traversal}.Station
%           data.traversal{i_traversal}.Yaw
%      Note that i_path denotes an array of traversals. Each traversal will
%      be plotted separately.
%
% OUTPUTS:
%
%      h: a handle to the resulting figure
%
% EXAMPLES:
%      
%       See the script: script_test_fcn_Path_plotPathYaw.m for a full test
%       suite. 
%
% This function was written on 2020_11_12 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%      2020_11_12 
%      -- first wrote the code
%      2021_01_06
%      -- added functions for input checking
%      -- renamed function for traversals
%      -- fixed error in yaw plotting (yaw is shorter than station!)
%      2021_01_07
%      -- fixed error in yaw plotting


flag_do_debug = 0; % Flag to plot the results for debugging
flag_this_is_a_new_figure = 1; % Flag to check to see if this is a new figure
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'Starting function: %s, in file: %s\n',st(1).name,st(1).file);
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
    fcn_Path_checkInputsToFunctions(data, 'traversals');

end

% Does user want to show the plots?
if 2 == nargin
    fig_num = varargin{1};
    figure(fig_num);
    flag_this_is_a_new_figure = 0;
else    
    fig = figure;
    fig_num = fig.Number;
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
 
figure(fig_num);

% Check to see if hold is already on. If it is not, set a flag to turn it
% off after this function is over so it doesn't affect future plotting
flag_shut_hold_off = 0;
if ~ishold
    flag_shut_hold_off = 1;
    hold on
end

NumTraversals = length(data.traversal);
h = zeros(NumTraversals,1);
for i_path= 1:NumTraversals
    h(i_path) = plot(data.traversal{i_path}.Station(1:end-1,1),data.traversal{i_path}.Yaw*180/pi);
end

% Shut the hold off?
if flag_shut_hold_off
    hold off;
end

% Add labels? 
if flag_this_is_a_new_figure == 1
    title('Station vs Yaw')
    xlabel('Station [m]')
    ylabel('Yaw [degree]')
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
if flag_do_debug
    % Nothing in here yet
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end
end

