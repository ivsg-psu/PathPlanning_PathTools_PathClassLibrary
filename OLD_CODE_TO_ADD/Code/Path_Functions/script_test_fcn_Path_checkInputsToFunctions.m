% script_test_fcn_Path_checkInputsToFunctions.m
% Tests fcn_Path_checkInputsToFunctions
       
% Revision history:
%      2021_01_06
%      -- first write of the code
%      2021_01_07:
%      -- denoted sections with figlets
%      -- created 'path' and 'paths' checks

%% Station
%    _____ _        _   _             
%   / ____| |      | | (_)            
%  | (___ | |_ __ _| |_ _  ___  _ __  
%   \___ \| __/ _` | __| |/ _ \| '_ \ 
%   ____) | || (_| | |_| | (_) | | | |
%  |_____/ \__\__,_|\__|_|\___/|_| |_|
%                                     
                                    

%% Test the station type (success)
station_test = 4;
fcn_Path_checkInputsToFunctions(station_test, 'station');


%% Test the station type (fail since non-numeric)
station_test = 'junk';
fcn_Path_checkInputsToFunctions(station_test, 'station');


%% Test the station type (fail since not 1 column)
station_test = [4 0];
fcn_Path_checkInputsToFunctions(station_test, 'station');

%% Stations
% 
%    _____ _        _   _                 
%   / ____| |      | | (_)                
%  | (___ | |_ __ _| |_ _  ___  _ __  ___ 
%   \___ \| __/ _` | __| |/ _ \| '_ \/ __|
%   ____) | || (_| | |_| | (_) | | | \__ \
%  |_____/ \__\__,_|\__|_|\___/|_| |_|___/
%                                         
%                                         


%% Test the stations type (success)
station_test = [4; 2];
fcn_Path_checkInputsToFunctions(station_test, 'stations');


%% Test the stations type (fail since only one row)
station_test = 4;
fcn_Path_checkInputsToFunctions(station_test, 'stations');

%% Path
% 
%   _____      _   _     
%  |  __ \    | | | |    
%  | |__) |_ _| |_| |__  
%  |  ___/ _` | __| '_ \ 
%  | |  | (_| | |_| | | |
%  |_|   \__,_|\__|_| |_|
%                        
%                        


%% Test the path type (success)
clc;
path_test = [4 1; 2 1];
fcn_Path_checkInputsToFunctions(path_test, 'path');


%% Test the path type (fail since only one column)
clc
path_test = [4; 2];
fcn_Path_checkInputsToFunctions(path_test, 'path');

%% Test the path type (fail since only one row)
clc
path_test = [4 2];
fcn_Path_checkInputsToFunctions(path_test, 'path');

%% Paths
% 
%   _____      _   _         
%  |  __ \    | | | |        
%  | |__) |_ _| |_| |__  ___ 
%  |  ___/ _` | __| '_ \/ __|
%  | |  | (_| | |_| | | \__ \
%  |_|   \__,_|\__|_| |_|___/
%                            
%                            


%% Test the paths type (success)
clc;
paths_test = [4 1; 2 1; 3 2];
fcn_Path_checkInputsToFunctions(paths_test, 'paths');


%% Test the paths type (fail since only one column)
clc
paths_test = [4; 2];
fcn_Path_checkInputsToFunctions(paths_test, 'paths');

%% Test the paths type (fail since only two rows)
clc
paths_test = [4 2; 0 0];
fcn_Path_checkInputsToFunctions(paths_test, 'paths');

%% Traversal
% 
%   _______                                 _ 
%  |__   __|                               | |
%     | |_ __ __ ___   _____ _ __ ___  __ _| |
%     | | '__/ _` \ \ / / _ \ '__/ __|/ _` | |
%     | | | | (_| |\ V /  __/ |  \__ \ (_| | |
%     |_|_|  \__,_| \_/ \___|_|  |___/\__,_|_|
%                                             
%                                             


%% Test the traversal type (success)
% Fill in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
for i_Path = 1:length(paths)
    test_traversal = fcn_Path_convertPathToTraversalStructure(paths{i_Path});
    all_traversals.traversal{i_Path} = test_traversal;
end
fcn_Path_checkInputsToFunctions(test_traversal, 'traversal');


%% Test the traversal type (fail since field is missing)
clc
clear test_traversal_bad
test_traversal_bad = 3;
fcn_Path_checkInputsToFunctions(test_traversal_bad, 'traversal');

%% Test the traversal type (fail since field is not numeric)
clc;
clear test_traversal_bad
test_traversal_bad.X = 'junk';
test_traversal_bad.Y = 'junk';
test_traversal_bad.Station = 'junk';

fcn_Path_checkInputsToFunctions(test_traversal_bad, 'traversal');

%% Test the traversal type (fail since fields are not columns)
clc;
clear test_traversal_bad
test_traversal_bad.X = eye(2);
test_traversal_bad.Y = eye(2);
test_traversal_bad.Station = eye(2);

fcn_Path_checkInputsToFunctions(test_traversal_bad, 'traversal');

%% Test the traversal type (fail since fields have different lengths)
clc;
clear test_traversal_bad
test_traversal_bad.X = [1; 2; 3];
test_traversal_bad.Y = [1; 2; 3];
test_traversal_bad.Station = [1; 2];

fcn_Path_checkInputsToFunctions(test_traversal_bad, 'traversal');

%% Test the traversal type (fail since Station field is not strictly increasing)
clc;
clear test_traversal_bad
test_traversal_bad.X = [1; 2; 3];
test_traversal_bad.Y = [1; 2; 3];
test_traversal_bad.Station = [1; 2; 2];

fcn_Path_checkInputsToFunctions(test_traversal_bad, 'traversal');

%% Traversals
% 
%   _______                                 _     
%  |__   __|                               | |    
%     | |_ __ __ ___   _____ _ __ ___  __ _| |___ 
%     | | '__/ _` \ \ / / _ \ '__/ __|/ _` | / __|
%     | | | | (_| |\ V /  __/ |  \__ \ (_| | \__ \
%     |_|_|  \__,_| \_/ \___|_|  |___/\__,_|_|___/
%                                                 
%                                                 

%% Test the traversals type (success)
% Fill in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
for i_Path = 1:length(paths)
    test_traversal = fcn_Path_convertPathToTraversalStructure(paths{i_Path});
    test_traversals.traversal{i_Path} = test_traversal;
end
clc;
fcn_Path_checkInputsToFunctions(test_traversals, 'traversals');


%% Test the traversals type (fail since field is missing)
clc
clear test_traversals_bad
test_traversals_bad = 3;
fcn_Path_checkInputsToFunctions(test_traversals_bad, 'traversals');

%% Test the traversals type (fail since second index is bad)
clc
clear test_traversals_bad
test_traversals_bad.traversal{1} = test_traversal;

clear test_traversal_bad
test_traversal_bad.X = [1; 2; 3];
test_traversal_bad.Y = [1; 2; 3];
test_traversal_bad.Station = [1; 2; 2];

test_traversals_bad.traversal{2} = test_traversal_bad;
fcn_Path_checkInputsToFunctions(test_traversals_bad, 'traversals');


