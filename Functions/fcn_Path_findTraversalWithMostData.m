function  index_of_longest = fcn_Path_findTraversalWithMostData(data)
% fcn_Path_findTraversalWithMostData.m
% finds the traversal index with the most amount of data (determined as the
% most elements in the X array)
%
% FORMAT: 
%
%      index_of_longest = fcn_Path_findTraversalWithMostData(data)
%
% INPUTS:
%
%      data: a structure containing subfields of X in the
%      following form
%           data.traversal{i_path}.X
%      Note that i_path denotes an array of paths. Each path will be
%      compared separately
%
% OUTPUTS:
%
%      index: the index of the traversal with the most data
%
% EXAMPLES:
%      
%       See the script: script_test_fcn_Path_findTraversalWithMostData.m
%       for a full test suite. 
%
% This function was written on 2020_11_12 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%       2020_11_12 
%       - wrote the code
%       2021_01_02
%       - added more checks to traversal type
%     2022_01_06:
%     -- fixed comments in header


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

if flag_check_inputs == 1
    % Are there the right number of inputs?
    if nargin > 1
        error('Incorrect number of input arguments')
    end
    
    % Check the data input
    fcn_Path_checkInputsToFunctions(data, 'traversals');

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
 
% choose the initial reference path, choose the path with most data points 
data_length = zeros(1,length(data.traversal));
for i_path = 1:length(data.traversal)
    data_length(i_path) = length(data.traversal{i_path}.X);
end
[~,index_of_longest] = max(data_length);


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

