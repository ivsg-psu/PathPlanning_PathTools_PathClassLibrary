function fcn_DebugTools_debugPrintTableToNCharacters(...
    table_data, header_strings, formatter_strings, N_chars,varargin)
% fcn_DebugTools_debugPrintTableToNCharacters
% Given a matrix of data, prints the data in user-specified width to the
% workspace.
%
% FORMAT:
%
%      fcn_DebugTools_debugPrintTableToNCharacters(...
%         table_data, header_strings, formatter_strings, N_chars,varargin)
%
% INPUTS:
%
%      table_data: a matrix of N rows, M columns, containing data to be
%      printed
%
%      header_strings: a cell array, M long, containing characters to be
%      printed as headers
%
%      formatter_strings: a cell array, M long, containing the print
%      specification for each column
%
%      N_chars: an integeter saying how long the output string should be if
%      columns have constant spacing, or an array of M integers with each
%      integer corresponding to the print with of the respective column.
%
% OUTPUTS:
%
%      (none)
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
% See the script: script_test_fcn_DebugTools_debugPrintTableToNCharacters
% for a full test suite.
%
% This function was written on 2023_01_17 by S. Brennan.
% Questions or comments? sbrennan@psu.edu

% Revision history:
%      2023_01_17:
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
    
    % Are there the right number of inputs?
    narginchk(4,5);
    
    % Check the table_data input
    %fcn_DebugTools_checkInputsToFunctions(table_data, '_of_chars');

    % Check the header_strings input
    %fcn_DebugTools_checkInputsToFunctions(header_strings, '_of_chars');

    % Check the N_chars input   
    fcn_DebugTools_checkInputsToFunctions(N_chars, '_of_integers');

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

N_header_strings = length(header_strings);
N_data_rows = length(table_data(:,1));

% Print the header
fprintf(1,'\n\n');
for ith_header = 1:N_header_strings
    header_str = header_strings{ith_header};
    if length(N_chars)==N_header_strings
         fixed_header_str = fcn_DebugTools_debugPrintStringToNCharacters(header_str,N_chars(ith_header));
    else
        fixed_header_str = fcn_DebugTools_debugPrintStringToNCharacters(header_str,N_chars);
    end
    fprintf(1,'%s ', fixed_header_str);
end
fprintf(1,'\n');

% Print the results
if ~isempty(table_data)
    for ith_row =1:N_data_rows
        for jth_col = 1:N_header_strings
            data_str = sprintf(formatter_strings{jth_col},table_data(ith_row,jth_col));
            if length(N_chars)==N_header_strings
                fixed_data_str = fcn_DebugTools_debugPrintStringToNCharacters(data_str,N_chars(jth_col));
            else
                fixed_data_str = fcn_DebugTools_debugPrintStringToNCharacters(data_str,N_chars);
            end
            fprintf(1,'%s ',...
                fixed_data_str);
        end % ends looping through columns
        fprintf(1,'\n');
        
    end % Ends looping down the rows
end % Ends check to see if table isempty



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

