% script_test_fcn_DebugTools_addSubdirectoriesToPath
% Tests: fcn_DebugTools_addSubdirectoriesToPath.m

% 
% REVISION HISTORY:
%      2022_03_27:
%      -- first write of the code



% Show how to add user-defined folders to the path

%% Successful example
% Save current folder and move into the top folder, then move back into
% current folder
current_folder = pwd;
cd ..
top_folder = pwd;
cd(current_folder);

% Save the "which" command result
full_location = fullfile(top_folder,'Data','data_DebugToolsExample.m');

data_directory = fullfile(top_folder,'Data');
rmpath(data_directory)

% Show the file can't be found now
blank_response = which('data_DebugToolsExample');
assert(isempty(blank_response)); % Show that the file can't be found

% Add the path now
fcn_DebugTools_addSubdirectoriesToPath(top_folder,{'Functions','Data'});

% Show that the file can be found!
good_response = which('data_DebugToolsExample');
assert(isequal(good_response,full_location)); % Show that the file CAN be found

%% FAIL Case: try to add a bad folder
% Show that the file can be found!
fcn_DebugTools_addSubdirectoriesToPath(top_folder,{'BadName'});

