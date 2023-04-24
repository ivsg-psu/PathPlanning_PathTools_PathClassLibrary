function [Flag_StringsMatch] = fcn_DebugTools_doStringsMatch(...
    InputAnswer,...
    CorrectAnswer)

% fcn_DebugTools_doStringsMatch
% Checks if an input string matches the "correct answer" string.
% This handles as well situations where the input, such as 'a',
% is a member of possible answers 'abd'. The function returns true if they
% match, false if not. It will also reject "spam" inputs such as 'aaa' for
% inputs 'abc'.
%
% FORMAT:
%
%     [Flag_StringsMatch] = fcn_DebugTools_doStringsMatch(...
%         InputAnswer,...
%         CorrectAnswer)
%
% INPUTS:
%
%      InputAnswer: a string that is to be tested
%
%      CorrectAnswer: the string that contains the input answer
%
% OUTPUTS:
%
%      flag_stringsMatch
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
%     See the script: script_test_fcn_DebugTools_doStringsMatch
%     for a full test suite.
%
% This function was written on 2022_12_09 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2022_12_09:
% -- wrote the code originally
% 2023_01_17:
% -- Moved out of the AutoExam codeset, into DebugTools
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
    narginchk(2,2);

    %     % Check the query_path input, 2 or 3 columns, 1 or more rows
    %     fcn_DebugTools_checkInputsToFunctions(query_path, '2or3column_of_numbers',[1 2]);
    %
    %     % Check the zone_center input, 2 or 3 columns, 1 row
    %     fcn_DebugTools_checkInputsToFunctions(zone_center, '2or3column_of_numbers',[1 1]);
    %
    %     % Check the zone_radius input, 1 column, 1 row
    %     fcn_DebugTools_checkInputsToFunctions(zone_radius, 'positive_1column_of_numbers',[1 1]);


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


try
    % Cannot let the student give more answers than allowed. If this is the
    % case, then automatically wrong. If the strings are same length, then
    % check if any of them match.
    if length(InputAnswer) > length(CorrectAnswer) % Are the student answers longer than the problem allows? If so, wrong
        Flag_StringsMatch = false;
    elseif length(InputAnswer)==1 && any(lower(CorrectAnswer) == lower(InputAnswer)) 
        % This format is used for Multiple Choice: ('a', 'b', or 'c'). If so, the
        % student's answer will be a single character and the correct
        % answer will be a string that may be more than one character,
        % since more than one can be right. If so, check: do the strings match? If so, right!
        % NOTE: this also marks correct answers such as "y" when answer is
        % "yes", or "n" when answer is "no".
        Flag_StringsMatch = true;
    elseif strcmpi(InputAnswer,CorrectAnswer) % Are strings non-zero length, but match exactly? If so, right!
        Flag_StringsMatch = true;
    else % The strings are same length, but do not match so this is wrong
        Flag_StringsMatch = false;
    end
catch
    error('A string condition was found that threw an error!')
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

%% fcn_DebugTools_doStringsMatch




