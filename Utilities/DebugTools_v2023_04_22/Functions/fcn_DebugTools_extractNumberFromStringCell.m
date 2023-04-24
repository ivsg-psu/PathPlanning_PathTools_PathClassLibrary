function outputNumberCell = fcn_DebugTools_extractNumberFromStringCell(inputStringCell)
%% fcn_DebugTools_extractNumberFromStringCell
% Takes an input cell that typically contains a string, and extracts the
% first number it finds. Useful to convert student entries such as {'1
% turtle'} into data, e.g. {'1'}
%
% FORMAT:
%
%      number = fcn_DebugTools_extractNumberFromStringCell(inputStringCell)
%
% INPUTS:
%
%      inputStringCell: a cell type containing a string that may include
%      numbers and text
%
% OUTPUTS:
%
%      numberCell: a cell type containing only the first number part of the input
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
%
% EXAMPLES:
%
%     See the script: script_test_fcn_DebugTools_extractNumberFromStringCell
%     for a full test suite.
%
% This function was written on 2022_11_14 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2022_11_14:
% -- wrote the code originally by copying out of old Exam2 code

% TO DO
% -- Add test scripts
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
    narginchk(1,1);

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

% Initialize output
outputNumberCell = {};

% Convert input to characters
StudentAnswer_string = char(inputStringCell);

if ~isempty(StudentAnswer_string)

    % We use regular expression, regexp, to extract the number portion of
    % the string. The challenge is that the number can be simple, such as
    % "1" or difficult such as "1,234.50".
    % The following works for any number with commas or decimals - see
    % https://www.mathworks.com/matlabcentral/answers/34548-regexp-help
    % See also: https://www.mathworks.com/matlabcentral/fileexchange/48930-interactive-regular-expression-tool
    expression = '(\d+,)*\d+(\.\d*)?'; 
    match = regexp(StudentAnswer_string,expression,'match','forceCellOutput');

    % Did we find anything?
    if ~isempty(match)
        numberCell = match{1};

        % Did we get more than one cell? If so, only keep the first hit
        if length(numberCell)>1
            numberCell = {numberCell{1}};
        end

        % Did we get an empty result? If not, check for leading decimals,
        % zeros and negative signs
        if ~isempty(numberCell)
            % Pull out just the string
            numberCell_string = numberCell{1};

            % Find the start index of the match
            startIndex = strfind(StudentAnswer_string,numberCell_string);
            startIndex = startIndex(1); % Keep only the first one

            % Do an error check
            if isempty(startIndex)
                error('An expression was returned from regexp that is not found in strfind');
            end

            % Is there a leading decimal point? If so, add this back into
            % the number.
            if startIndex~=1 && (StudentAnswer_string(startIndex-1)=='.')
                numberCell_string = cat(2,'.',numberCell_string);

                % Find the start index of the match
                startIndex = strfind(StudentAnswer_string,numberCell_string);
                startIndex = startIndex(1); % Keep only the first one
            end
            

            % Is there a leading negative sign? If so, note it as we have
            % to add it later.
            if startIndex~=1 && (StudentAnswer_string(startIndex-1)=='-')
                flag_is_negative = 1;
            else
                flag_is_negative = 0;
            end

            % After this point, the number should only be positive, and only of
            % characters 0 to 9, commas, and decimals.

            % To remove leading zeros, we search for 0 in the front of a
            % number, and start the number if the character is NOT zero.
            % However, there's one special condition where we keep a
            % leading zero, and that's when the next character is a
            % decimal, as in the following: '0.4'. We go ahead and strip
            % the zero in this special condition, but we need to add it
            % later. See steps that follow.
            starting_index = 1;
            flag_number_started = 0;
            for ith_character = 1:(length(numberCell_string)-1)
                if numberCell_string(ith_character)~='0' && (flag_number_started==0)
                    starting_index = ith_character;
                    flag_number_started = 1;
                    %                 elseif numberCell_string(ith_character)=='0' && numberCell_string(ith_character+1)=='.' &&  (flag_number_started==0)
                    %                     % This is the special case of numbers such as '0.4'
                    %                     starting_index = ith_character;
                    %                     flag_number_started = 1;
                end
            end
            if flag_number_started==0
                starting_index = length(numberCell_string);
            end
            good_numberCell_string = numberCell_string(starting_index:end);

            % Add a leading zero back in?
            if good_numberCell_string(1)=='.'
                good_numberCell_string = cat(2,'0',good_numberCell_string);
            end

            % Add a leading negative sign back in?
            if flag_is_negative
                good_numberCell_string = cat(2,'-',good_numberCell_string);
            end
            outputNumberCell = {good_numberCell_string};
        end % Ends check if Numbercell is empty
    end % Ends check if match is empty
end % Ends check if string is empty


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

