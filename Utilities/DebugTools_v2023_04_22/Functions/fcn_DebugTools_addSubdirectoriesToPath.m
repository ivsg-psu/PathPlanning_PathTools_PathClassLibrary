function fcn_DebugTools_addSubdirectoriesToPath(root_path, subdirectories)
% fcn_DebugTools_addSubdirectoriesToPath
% Given a path, and a cell array of subdirectories, adds the subdirectories
% to the path, if they exist under the given path. It throws an error if
% the subdirectories do not exist
%
% FORMAT:
%
%      fcn_DebugTools_addSubdirectoriesToPath(root_path,subdirectories)
%
% INPUTS:
%
%      root_path: a string of the path location
% 
%      subdirectories: a cell array of strings containing subdirectories,
%      for example {'Functions','Data','Utilities'}
%
% OUTPUTS:
%
%      (none)
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
% See the script: fcn_DebugTools_addSubdirectoriesToPath
% for a full test suite.
%
% This function was written on 2022_03_27 by S. Brennan 
% Questions or comments? sbrennan@psu.edu

% Revision history:
%      2022_03_27:
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


% add necessary directories
if isempty(subdirectories)
    if(exist(root_path,'dir'))
        addpath(root_path);
    else % Throw an error?
        error('There was an attempt to add directory: \n%s \nto the path but the subdirectory does not exist. Suggest checking the README file to ensure correct folders are included.',...
            root_path);
    end
else
    for ith_subdirectory = 1:length(subdirectories)
        subdirectory_name = subdirectories{ith_subdirectory};
        if(exist([root_path, filesep,  subdirectory_name],'dir'))
            addpath(genpath([root_path, filesep, subdirectory_name]))
        else % Throw an error?
            error('There was an attempt to add subdirectory: \n%s \nto the path: \n%s\nbut the subdirectory does not exist. Suggest checking the README file to ensure correct folders are included.',...
                subdirectory_name,root_path);
        end
    end
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

