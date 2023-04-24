function stringNumber = fcn_DebugTools_number2string(number)

% fcn_DebugTools_number2string
% Converts a number into a string format where the printed string is "nice"
% looking, e.g. small numbers have just a few decimal places, large numbers
% don't have decimal places.
%
% FORMAT:
%
%      stringNumber = fcn_DebugTools_number2string(number)
%
% INPUTS:
%
%      number: a number to be converted
%
% OUTPUTS:
%
%      stringNumber: a string representation of a number
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
%
% EXAMPLES:
%
%     See the script: script_test_fcn_DebugTools_number2string
%     for a full test suite.
%
% This function was written on 2022_11_14 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
%
% 2022_11_14:
% -- wrote the code originally by copying out of old Exam2 code
% 2021_12_12:
% -- first write of the script
% 2023_02_17
% -- copied code out of AutoExam and into DebugTools
% -- Add test scripts
% -- Add input argument checking

% TO DO
% -- vectorize code, allowing matrices of numbers?

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

    % Check the number input, 1 columns, 1 rows
    if ~isempty(number)
        fcn_DebugTools_checkInputsToFunctions(number, '1column_of_numbers',[1 1]);
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
try
    if ~isempty(number)
        if abs(number)<(100*eps) % Is it zero?
            stringNumber = sprintf('%1.0f',0);
        elseif all(round(number)==number) && abs(number)<1000 % Is it a small integer?
            stringNumber = sprintf('%.0d',number);
        elseif abs(number)<1 % is it a very small number? If so, print it
            stringNumber = sprintf('%f',number);
        elseif abs(number)<10 % is it a float less than 10? If so, print 2 decimals
            stringNumber = sprintf('%.2f',number);
        elseif abs(number)<100 % is it a float less than 100? If so, print 1 decimals
            stringNumber = sprintf('%.1f',number);
        else % it is a very large number, print NO decimals
            stringNumber = sprintf('%.0f',number);
        end
    else
        stringNumber = sprintf(' ');
    end
catch
    error('debug here');
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
