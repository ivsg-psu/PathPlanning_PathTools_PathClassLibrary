function fixed_length_string = fcn_Path_debugPrintStringToNCharacters(input_sequence,N)
% fcn_Path_debugPrintStringToNCharacters
% Given a string and an integer N representing the number of characters to
% keep or pad, creates a new string of exactly length N by cropping the
% string or padding it (to the right) with spaces.
%
% FORMAT:
%
%      fixed_length_string = ...
%      fcn_Path_debugPrintStringToNCharacters(input_sequence,N)
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
%      (none)
%
% EXAMPLES:
%
% See the script: script_test_fcn_Path_debugPrintStringToNCharacters
% for a full test suite.
%
% This function was written on 2021_01_24 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
%      2021_01_24:
%      -- first write of the code


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
    %     % Are there the right number of inputs?
    %     if nargin < 2 || nargin > 3
    %         error('Incorrect number of input arguments')
    %     end
    %
    %     % Check the traversal_1 input
    %     fcn_Path_checkInputsToFunctions(traversal_1, 'traversal');
    %
    %     % Check the traversal_2 input
    %     fcn_Path_checkInputsToFunctions(traversal_2, 'traversal');
              
end


% % Does user want to show the plots?
% if 3 == nargin
%     fig_num = varargin{1};
%     figure(fig_num);
%     flag_do_plot = 1;
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

