
function output_string = fcn_DebugTools_addStringToEnd(input_string,value_to_add,varargin)
%% fcn_DebugTools_addStringToEnd
% Adds information to a string. The input is the starter string.
% The value to add can be a cell array, string, or numeric. 
%
% If it is a cell array, the first index of the cell is appended. An
% optional index value can be given to specify a different cell.
%
% If it is a string or character, then it is appended
%
% If it is numeric, the numeric value is formatted to be "pretty" to read
% using fcn_AutoExam_number2string functionality
%
% FORMAT:
%
%      output_string = fcn_DebugTools_addStringToEnd(input_string,value_to_add,(index))
%
% INPUTS:
%
%      input_string: the string to start with
%
% OUTPUTS:
%
%      value_to_add: the numeric, string, or cell value to add
%
%      (OPTIONAL INPUTS)
%      index: the index of a call array, if passed as the value_to_add, to use
%
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
%
% EXAMPLES:
%
%     See the script: script_test_fcn_DebugTools_addStringToEnd
%     for a full test suite.
%
% This function was written on 2022_11_14 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2022_11_14:
% -- wrote the code originally by copying out of old Exam2 code
% 2023_01_17:
% -- added code to the DebugTools repo
% -- Add test scripts

% TO DO
% -- Add input argument checking

flag_do_debug = 0; % Flag to show the results for debugging
flag_do_plots = 0; % % Flag to plot the final results
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

if flag_check_inputs
    % Are there the right number of inputs?
    narginchk(2,3);
end

if nargin>2
    index_value = varargin{1};
else
    index_value = 1;
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


if iscell(value_to_add)
    output_string = cat(2,input_string,' ',value_to_add{index_value});
elseif isstring(value_to_add) || ischar(value_to_add)
    output_string = cat(2,input_string,' ',value_to_add);    
elseif isnumeric(value_to_add)
    output_string = cat(2,input_string,' ',fcn_DebugTools_number2string(value_to_add));
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

    % Nothing to do here


end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends main function




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
