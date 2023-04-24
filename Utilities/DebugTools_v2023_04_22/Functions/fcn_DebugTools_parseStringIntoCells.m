function [outputCells,outputCellsString] = fcn_DebugTools_parseStringIntoCells(inputString)
%% fcn_DebugTools_parseStringIntoCells(inputString)
% Takes an input string, and distributes the outputs into cells.
%
% FORMAT:
%
%      [outputCells,outputCellsString] = fcn_DebugTools_parseStringIntoCells(inputString)
%
% INPUTS:
%
%      inputString: the value to start with
%
% OUTPUTS:
%
%
%      outputCells: an array of cells representing the parsing of the input
%      string
%
%      outputCellsString: a string created from the inputs, pushed into a cell
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
%     See the script: script_test_fcn_DebugTools_parseStringIntoCells
%     for a full test suite.
%
% This function was written on 2022_11_14 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2022_11_14:
% -- wrote the code originally by copying out of old Exam2 code
% 2023_01_18:
% -- migrated code out of AutoExam and into DebugTools
% -- Added improved comments
% -- Add test scripts
% -- Add input argument checking

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

    % Check the inputString input, make sure is character type
    fcn_DebugTools_checkInputsToFunctions(inputString, '_of_chars');


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

% Check for complex string cases in answer

% Parse answer by commas OR by spaces
expression = '[\s\,]'; % Split by comma (\c) or any whitespace character (\s)
outputRegexp = regexp(inputString,expression,'split','forceCellOutput');
if ~isempty(outputRegexp)
    outputCellArray = outputRegexp{1};
else
    outputCellArray = '';    
end


good_output = [];
for ith_cell = 1:length(outputCellArray)
    if ~isempty(outputCellArray{ith_cell})
        good_output = [good_output,{outputCellArray{ith_cell}}]; %#ok<AGROW> 
    end
end

outputCells = good_output;

% Construct a string of correct answers - this is useful because sometimes
% the student selects "C" out of answers "B" "C", and it could get counted
% wrong unless we check each student answer is checked against the entire
% list
outputCellsString = [];
for ith_part = 1:length(outputCells)
    outputCellsString = [outputCellsString,lower(outputCells{ith_part})]; %#ok<AGROW>
end

% 
% % Check to see if student answer  is a number. If so, print student answer
% % to string number (no decimal place) if the correct answer is a string.
% StudentAnswer = regexprep(StudentAnswerIn,' +',' '); % Remove extra spaces, if any, in the student answer
% cleanedAnswer_string = StudentAnswer;
% if ischar(CorrectAnswer) && (~ischar(cleanedAnswer_string))
%     cleanedAnswer_string = sprintf('%.0f',StudentAnswer);
% elseif isnumeric(cleanedAnswer_string) && length(StudentAnswer) == 1
%     cleanedAnswer_string = sprintf('%f',StudentAnswer);
% end
% 
% % Check to see if there are spaces in the correct answer
% if length(CorrectAnswer)>1
%     cleanedAnswer_cellarray = regexp(char(cleanedAnswer_string),'\s','split');
% else
%     cleanedAnswer_cellarray = {cleanedAnswer_string};
% end
% 



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
