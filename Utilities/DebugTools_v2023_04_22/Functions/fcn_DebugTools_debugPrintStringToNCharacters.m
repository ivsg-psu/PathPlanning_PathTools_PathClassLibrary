function fixed_length_string = fcn_DebugTools_debugPrintStringToNCharacters(input_sequence,N)
% fcn_DebugTools_debugPrintStringToNCharacters
% Given a string and an integer N representing the number of characters to
% keep or pad, creates a new string of exactly length N by cropping the
% string or padding it (on the right) with spaces.
%
% FORMAT:
%
%      fixed_length_string = ...
%      fcn_DebugTools_debugPrintStringToNCharacters(input_sequence,N)
%
% INPUTS:
%
%      input_sequence: a string type of variable length
%
%      N: an integeter saying how long the output string should be
%
% OUTPUTS:
%
%      a string of length N
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
% See the script: script_test_fcn_DebugTools_debugPrintStringToNCharacters
% for a full test suite.
%
% This function was written on 2021_12_12 by S. Brennan from same function
% in PathTools repo.
% Questions or comments? sbrennan@psu.edu

% Revision history:
%      2021_12_12:
%      -- first write of the code
%      2023_01_16:
%      -- added input checking


flag_do_debug = 0; % Flag to debug the results
flag_do_plot = 0;  % Flag to plot the results
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
    narginchk(2,3);
    
    % Check the input_sequence input
    fcn_DebugTools_checkInputsToFunctions(input_sequence, '_of_chars');

    % Check the N input   
    fcn_DebugTools_checkInputsToFunctions(N, '1column_of_integers',1);

end


% % Does user want to show the plots?
% if 3 == nargin
%     temp = varargin{1};
%     if ~isempty(temp)
%         figure(fig_num);
%         flag_do_plot = 1;
%     end
% else
%     if flag_do_debug
%         fig = figure;
%         fig_num = fig.Number;
%         flag_do_plot = 1;
%     end
% end

%% Start of main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



input_string = char(input_sequence); % Convert to a character array
if strlength(input_string) == N % It is just right!
    fixed_length_string = input_string;
elseif strlength(input_string) > N  % Make it shorter?
    fixed_length_string = input_string(1:N);
else % Make it longer
    fixed_length_string = pad(input_string,N);
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
if 1 == flag_do_plot
    %     figure(fig_num);
    %     clf;
    %     hold on;
    %     grid on;
    %
    %        
    
end % Ends the flag_do_plot if statement

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends the function

