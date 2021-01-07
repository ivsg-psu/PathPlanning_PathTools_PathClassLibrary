function pathSXY = fcn_Path_convertXYtoSXY(X,Y)
% fcn_Path_convertXYtoSXY
% Takes XY positions and creates a SXY path where S is the path distance
%
% FORMAT: 
%
%     traversal = fcn_Path_convertXYtoSXY(X,Y)
%
% INPUTS:
%
%     X: an N x 1 vector of X positions
%     Y: an N x 1 vector of Y positions
%
% OUTPUTS:
%     pathSXY: a [N x 3] vector containing:
%         S: an N x 1 vector of Stations
%         X: an N x 1 vector of X positions
%         Y: an N x 1 vector of Y positions
%
% EXAMPLES:
%      
%     See the script:
%     script_test_fcn_Path_convertXYtoSXY.m for a full
%     test suite. 
%
% This function was written on 2020_11_16 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%     2020_11_16:
%     - wrote the code


flag_do_debug = 0; % Flag to plot the results for debugging
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

Npoints = length(X(:,1));
if flag_check_inputs == 1
    % Are there the right number of inputs?
    if nargin < 2 || nargin > 2
        error('Incorrect number of input arguments')
    end
    
    if length(Y(:,1))~= Npoints
        error('Y input must have same number of rows as X');
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

path = [X, Y];
path_dist = sum(diff(path).^2,2).^0.5;
path_S = [0; cumsum(path_dist)];
pathSXY = [path_S path];


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
