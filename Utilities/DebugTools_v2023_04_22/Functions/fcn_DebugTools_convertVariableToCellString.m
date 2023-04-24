function cellString = fcn_DebugTools_convertVariableToCellString(input_variable)
%% fcn_DebugTools_convertVariableToCellString
% Takes a variable of almost any type, converts it into a single string
% that is then put into a cell format. Allows multi-cell entries to be
% saved as just one string entry, separated by spaces.
%
% FORMAT:
%
%      cellString = fcn_DebugTools_convertVariableToCellString(input_variable)
%
% INPUTS:
%
%      input_variable: a numeric, character, or cell array
%
% OUTPUTS:
%
%      cellString: a string created from the inputs, pushed into a single cell
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_DebugTools_number2string
%
% EXAMPLES:
%
%     See the script: script_test_fcn_DebugTools_convertVariableToCellString
%     for a full test suite.
%
% This function was written on 2022_11_14 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2022_11_14:
% -- wrote the code originally by copying out of old Exam2 code
% 2023_01_18
% -- Add test scripts
% -- Add input argument checking
% -- Better comments

% TO DO
% -- (none)

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
    narginchk(1,1);

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

if iscell(input_variable)
    if length(input_variable)==1
        if ischar(input_variable{1})            
            cellString = input_variable{1};
        elseif isnumeric(input_variable{1})
            cellString = fcn_DebugTools_number2string(input_variable{1});
        else
            error('entered here in fcn_AutoExam_convertVariableToString -  a variable type was passed in that is not recognized');
        end
        cellString = {cellString};
    else
        % Build an answer we can store
        temp = ''; % Initialize an empty string
        for ith_cell = 1:length(input_variable)-1
            result = fcn_DebugTools_convertVariableToCellString(input_variable{ith_cell});
            temp = cat(2,temp,char(result),', ');
        end
        result = fcn_DebugTools_convertVariableToCellString(input_variable{length(input_variable)});
        temp = cat(2,temp,char(result)); % Last concatenation does not have space or comma at end
        cellString = {temp};
    end
elseif isnumeric(input_variable)
    cellString = {fcn_DebugTools_number2string(input_variable)};
elseif ischar(input_variable)
    cellString = {input_variable};
else
    error('entered here in fcn_AutoExam_convertVariableToString - this should not happen - why?');
end
if ~iscell(cellString)
    error('Should not have this happen!');
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
