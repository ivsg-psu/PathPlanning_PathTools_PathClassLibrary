% script_test_fcn_DebugTools_installDependencies.m
% Tests fcn_DebugTools_installDependencies
% Written in 2023_01_25 by S.Brennan

% Clear all prior global variable flags
clear global FLAG_*

flag_show_warnings = 0; % Set to 1 to see the warnings go by. Should keep off for people who are not familiar with code.

%% Basic test case
% NOTE: this installs under the current directory!
% Define the name of subfolder to be created in "Utilities" subfolder
dependency_name = 'DebugTools_v2023_01_25';

% Define sub-subfolders that are in the code package that also need to be
% added to the MATLAB path after install. Leave empty ({}) to only add
% the subfolder path without any sub-subfolder path additions.
dependency_subfolders = {'Functions','Data'};

% Define a universal resource locator (URL) pointing to the zip file to
% install. For example, here is the zip file location to the Debugtools
% package on GitHub:
dependency_url = 'https://github.com/ivsg-psu/Errata_Tutorials_DebugTools/blob/main/Releases/DebugTools_v2023_01_25.zip?raw=true';

% Call the function to do the install
fcn_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url)

disp('Library installed! Verify this now, as it will be deleted to complete the demo');
disp('Paused. Hit any key to continue...');
pause;

% Remove the folders from path, to avoid deletion warnings
temp_path = fullfile(pwd,'Utilities','DebugTools_v2023_01_25','Functions');
rmpath(temp_path);
temp_path = fullfile(pwd,'Utilities','DebugTools_v2023_01_25','Data');
rmpath(temp_path);

% Remove the example Utilities folder and all subfolders
[success_flag,error_message,message_ID] = rmdir('Utilities','s');

% Did it work?
if ~success_flag
    error('Unable to remove the example Utilities directory. Reason: %s with message ID: %s\n',error_message,message_ID);
elseif ~isempty(error_message)
    if flag_show_warnings
        warning('The Utilities directory was removed, but with a warning: %s\n and message ID: %s\n(continuing)\n',error_message, message_ID);
    end
end

%% Call the function again, to show that global flag blocks directory creation
% Call the function to do the install
fcn_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url)
disp('Library NOT installed! Verify this now, as it will be checked to complete the demo');
disp('Paused. Hit any key to continue...');
pause;

% Remove the example Utilities folder and all subfolders
[success_flag,error_message,message_ID] = rmdir('Utilities','s');

% Did it work?
if success_flag
    error('Was unexpectedly able to delete the Utilities directory. This should NOT have happened. Reason: %s with message ID: %s\n',error_message,message_ID);
elseif ~isempty(error_message)
    if flag_show_warnings
        warning('The Utilities directory was NOT removed, which is correct. We confirm this because MATLAB throws the warning: %s\n and message ID: %s\n(continuing)\n',error_message, message_ID);
    end
end

%% Call the function again, with override flag
% Call the function to do the install, with override
override = 1;
fcn_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url,override)

disp('Library installed! Verify this now, as it will be deleted to complete the demo');
disp('Paused. Hit any key to continue...');
pause;

% Remove the folders from path, to avoid deletion warnings
temp_path = fullfile(pwd,'Utilities','DebugTools_v2023_01_25','Functions');
rmpath(temp_path);
temp_path = fullfile(pwd,'Utilities','DebugTools_v2023_01_25','Data');
rmpath(temp_path);

% Remove the example Utilities folder and all subfolders
[success_flag,error_message,message_ID] = rmdir('Utilities','s');

% Did it work?
if ~success_flag
    error('Unable to remove the example Utilities directory. Reason: %s with message ID: %s\n',error_message,message_ID);
elseif ~isempty(error_message)
    if flag_show_warnings
        warning('The Utilities directory was removed, but with a warning: %s\n and message ID: %s\n(continuing)\n',error_message, message_ID);
    end
end


% Clear all prior global variable flags so that this function works in
% future
clear global FLAG_*

%% Fail cases
if 1==0
    % Incorrect arguments
    fcn_DebugTools_installDependencies(dependency_name);
end

