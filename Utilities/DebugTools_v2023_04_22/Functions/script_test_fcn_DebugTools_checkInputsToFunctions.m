% script_test_fcn_DebugTools_checkInputsToFunctions
% Tests: fcn_DebugTools_checkInputsToFunctions

% 
% REVISION HISTORY:
% 
% 2021_06_20 by S. Brennan
% -- first write of script
% 2021_07_07 by S. Brennan
% -- modified to allow variable input types
% 2022_04_03 by S. Brennan
% -- added Path variable types
% 2023_01_16 by S. Brennan
% -- added narginchk
% -- fixed dbstack error to match error to source function
% -- added char type
% -- added string type
% -- added DoesFileExist type
% -- added DoesDirectoryExist type
% -- improved documentation
%%%%%%%%%%%%%%§


%% Echo options
%   ______     _            ____        _   _                 
%  |  ____|   | |          / __ \      | | (_)                
%  | |__   ___| |__   ___ | |  | |_ __ | |_ _  ___  _ __  ___ 
%  |  __| / __| '_ \ / _ \| |  | | '_ \| __| |/ _ \| '_ \/ __|
%  | |___| (__| | | | (_) | |__| | |_) | |_| | (_) | | | \__ \
%  |______\___|_| |_|\___/ \____/| .__/ \__|_|\___/|_| |_|___/
%                                | |                          
%                                |_|                                                                
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=EchoOptions
options = fcn_DebugTools_checkInputsToFunctions;
fprintf(1,'Here are a listing of all active input checking options: \n');
for ith_option = 1:length(options)
    fprintf('\t"%s"\n',options(ith_option).name)
    fprintf('\t\t%s\n',options(ith_option).description)    
end

%% 1column_of_numbers
% 
%   __           _                               __                       _                   
%  /_ |         | |                             / _|                     | |                  
%   | | ___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ __  _   _ _ __ ___ | |__   ___ _ __ ___ 
%   | |/ __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| '_ \| | | | '_ ` _ \| '_ \ / _ \ '__/ __|
%   | | (_| (_) | | |_| | | | | | | | | || (_) | | | | | | |_| | | | | | | |_) |  __/ |  \__ \
%   |_|\___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_| |_|\__,_|_| |_| |_|_.__/ \___|_|  |___/
%                                     ______   ______                                         
%                                    |______| |______|                                        
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs

%% Test the column_of_numbers type (success)
column_of_numbers_test = 4;
fcn_DebugTools_checkInputsToFunctions(column_of_numbers_test, '1column_of_numbers');

column_of_numbers_test = [4; 3; 2];
fcn_DebugTools_checkInputsToFunctions(column_of_numbers_test, '1column_of_numbers');

column_of_numbers_test = [4; 3; 2];
fcn_DebugTools_checkInputsToFunctions(column_of_numbers_test, '1column_of_numbers',3);

%% positive_column_of_numbers
% 
%                   _ _   _            __           _                               __                       _                   
%                  (_) | (_)          /_ |         | |                             / _|                     | |                  
%   _ __   ___  ___ _| |_ ___   _____  | | ___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ __  _   _ _ __ ___ | |__   ___ _ __ ___ 
%  | '_ \ / _ \/ __| | __| \ \ / / _ \ | |/ __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| '_ \| | | | '_ ` _ \| '_ \ / _ \ '__/ __|
%  | |_) | (_) \__ \ | |_| |\ V /  __/ | | (_| (_) | | |_| | | | | | | | | || (_) | | | | | | |_| | | | | | | |_) |  __/ |  \__ \
%  | .__/ \___/|___/_|\__|_| \_/ \___| |_|\___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_| |_|\__,_|_| |_| |_|_.__/ \___|_|  |___/
%  | |                             ______                                ______   ______                                         
%  |_|                            |______|                              |______| |______|                                        
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs

%% Test the column_of_numbers type (success)
column_of_numbers_test = 4;
fcn_DebugTools_checkInputsToFunctions(column_of_numbers_test, 'positive_1column_of_numbers');

column_of_numbers_test = [4; 3; 2];
fcn_DebugTools_checkInputsToFunctions(column_of_numbers_test, 'positive_1column_of_numbers');

column_of_numbers_test = [4; 3; 2];
fcn_DebugTools_checkInputsToFunctions(column_of_numbers_test, 'positive_1column_of_numbers',3);


%% 2column_of_numbers
% 
%   ___           _                               __                       _                   
%  |__ \         | |                             / _|                     | |                  
%     ) |___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ __  _   _ _ __ ___ | |__   ___ _ __ ___ 
%    / // __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| '_ \| | | | '_ ` _ \| '_ \ / _ \ '__/ __|
%   / /| (_| (_) | | |_| | | | | | | | | || (_) | | | | | | |_| | | | | | | |_) |  __/ |  \__ \
%  |____\___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_| |_|\__,_|_| |_| |_|_.__/ \___|_|  |___/
%                                      ______   ______                                         
%                                     |______| |______|                                        
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs

%% Test the column_of_numbers type (success)
Twocolumn_of_numbers_test = [4 2];
fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers');

Twocolumn_of_numbers_test = [4 1; 3 0; 2 5];
fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers');

Twocolumn_of_numbers_test = [4 1; 3 9; 2 7];
fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',3);

% Minimum length is 2 or greater
Twocolumn_of_numbers_test = [4 1; 3 9; 2 7];
fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',[2 3]);

% Maximum length is 5 or less
Twocolumn_of_numbers_test = [4 1; 3 9; 2 7];
fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',[5 4]);

% Maximum length is 3
Twocolumn_of_numbers_test = [4 1; 3 9; 2 7];
fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',[3 3]);


%% 2or3column_of_numbers
% 
%   ___           ____            _                               __                       _                   
%  |__ \         |___ \          | |                             / _|                     | |                  
%     ) |___  _ __ __) | ___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ __  _   _ _ __ ___ | |__   ___ _ __ ___ 
%    / // _ \| '__|__ < / __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| '_ \| | | | '_ ` _ \| '_ \ / _ \ '__/ __|
%   / /| (_) | |  ___) | (_| (_) | | |_| | | | | | | | | || (_) | | | | | | |_| | | | | | | |_) |  __/ |  \__ \
%  |____\___/|_| |____/ \___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_| |_|\__,_|_| |_| |_|_.__/ \___|_|  |___/
%                                                      ______   ______                                         
%                                                     |______| |______|                                        
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs


%% Test the column_of_numbers type (success)
% Test 1 by 2
TwoOrThreeColumn_of_numbers_test = [4 2];
fcn_DebugTools_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers');

% Test 1 by 3
TwoOrThreeColumn_of_numbers_test = [4 2 1];
fcn_DebugTools_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers');

% Test multiple points - 2 columns
TwoOrThreeColumn_of_numbers_test = [4 1; 3 0; 2 5];
fcn_DebugTools_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers');

% Test multiple points - 3 columns
TwoOrThreeColumn_of_numbers_test = [4 1 3; 3 0 5; 2 5 7];
fcn_DebugTools_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers');

% Test specified length - 2 columns
TwoOrThreeColumn_of_numbers_test = [4 1; 3 9; 2 7];
fcn_DebugTools_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',3);

% Test specified length - 3 columns
TwoOrThreeColumn_of_numbers_test = [4 1 5; 3 9 5; 2 7 5];
fcn_DebugTools_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',3);

% Minimum length is 2 or greater - 2 columns
TwoOrThreeColumn_of_numbers_test = [4 1; 3 9; 2 7];
fcn_DebugTools_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[2 3]);

% Minimum length is 2 or greater - 3 columns
TwoOrThreeColumn_of_numbers_test = [4 1 5; 3 9 5; 2 7 5];
fcn_DebugTools_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[2 3]);

% Maximum length is 5 or less - 2 column
TwoOrThreeColumn_of_numbers_test = [4 1; 3 9; 2 7];
fcn_DebugTools_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[5 4]);

% Maximum length is 5 or less - 3 column
TwoOrThreeColumn_of_numbers_test = [4 1 4; 3 9 4; 2 7 4];
fcn_DebugTools_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[5 4]);

% Length MUST be 3 - 2 column
TwoOrThreeColumn_of_numbers_test = [4 1; 3 9; 2 7];
fcn_DebugTools_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[3 3]);

% Length MUST be 3 - 3 column
TwoOrThreeColumn_of_numbers_test = [4 1 3; 3 9 3; 2 7 3];
fcn_DebugTools_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[3 3]);

%% 2or3column_of_numbers
%   _   _            __  __           _                               __                       _                   
%  | \ | |          |  \/  |         | |                             / _|                     | |                  
%  |  \| | ___  _ __| \  / | ___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ __  _   _ _ __ ___ | |__   ___ _ __ ___ 
%  | . ` |/ _ \| '__| |\/| |/ __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| '_ \| | | | '_ ` _ \| '_ \ / _ \ '__/ __|
%  | |\  | (_) | |  | |  | | (_| (_) | | |_| | | | | | | | | || (_) | | | | | | |_| | | | | | | |_) |  __/ |  \__ \
%  |_| \_|\___/|_|  |_|  |_|\___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_| |_|\__,_|_| |_| |_|_.__/ \___|_|  |___/
%                                                          ______   ______                                         
%                                                         |______| |______|                                        
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=NorMcolumn_of_numbers
%% Test the N to M column_of_numbers type (success)
% Test 1 by 2
NorMColumn_of_numbers_test = [4 2];
fcn_DebugTools_checkInputsToFunctions(NorMColumn_of_numbers_test, '2or4column_of_numbers');

% Test 1 by 3
NorMColumn_of_numbers_test = [4 2 1];
fcn_DebugTools_checkInputsToFunctions(NorMColumn_of_numbers_test, '2or4column_of_numbers');


%% 2column_of_integers
% 
% 
%   ___           _                               __ _       _                           
%  |__ \         | |                             / _(_)     | |                          
%     ) |___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ _ __ | |_ ___  __ _  ___ _ __ ___ 
%    / // __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| | '_ \| __/ _ \/ _` |/ _ \ '__/ __|
%   / /| (_| (_) | | |_| | | | | | | | | || (_) | | | | | | | ||  __/ (_| |  __/ |  \__ \
%  |____\___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_|_| |_|\__\___|\__, |\___|_|  |___/
%                                      ______   ______                __/ |              
%                                     |______| |______|              |___/               
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs

%% Test the 2column_of_integers type (success)
Twocolumn_of_integers_test = [4 2];
fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers');

Twocolumn_of_integers_test = [4 1; 3 0; 2 5];
fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers');

Twocolumn_of_integers_test = [4 1; 3 9; 2 7];
fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers',3);

% Minimum length is 2 or greater
Twocolumn_of_integers_test = [4 1; 3 9; 2 7];
fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers',[2 3]);

% Maximum length is 5 or less
Twocolumn_of_integers_test = [4 1; 3 9; 2 7];
fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers',[5 4]);

% Maximum length is 3
Twocolumn_of_integers_test = [4 1; 3 9; 2 7];
fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers',[3 3]);



%% Char Types
%    _____ _                    
%   / ____| |                   
%  | |    | |__   __ _ _ __ ___ 
%  | |    | '_ \ / _` | '__/ __|
%  | |____| | | | (_| | |  \__ \
%   \_____|_| |_|\__,_|_|  |___/
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Chars                              

% Basic character test
Chars_test = 'abcdefg';
fcn_DebugTools_checkInputsToFunctions(Chars_test, '_of_chars');
                              

%% String Types
%    _____ _        _                 
%   / ____| |      (_)                
%  | (___ | |_ _ __ _ _ __   __ _ ___ 
%   \___ \| __| '__| | '_ \ / _` / __|
%   ____) | |_| |  | | | | | (_| \__ \
%  |_____/ \__|_|  |_|_| |_|\__, |___/
%                            __/ |    
%                           |___/     
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Strings

% Basic string test
String_test = "abcdefg";
fcn_DebugTools_checkInputsToFunctions(String_test, '_of_strings');

%% DoesFileExist
% 
%   ______ _ _      ______      _     _   
%  |  ____(_) |    |  ____|    (_)   | |  
%  | |__   _| | ___| |__  __  ___ ___| |_ 
%  |  __| | | |/ _ \  __| \ \/ / / __| __|
%  | |    | | |  __/ |____ >  <| \__ \ |_ 
%  |_|    |_|_|\___|______/_/\_\_|___/\__|
%                                         
%                                         
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=FileExist

% Basic file test
FileExist_test = 'fcn_DebugTools_checkInputsToFunctions.m';
fcn_DebugTools_checkInputsToFunctions(FileExist_test, 'DoesFileExist');


%% DoesDirectoryExist
%
%   _____  _               _                   ______      _     _
%  |  __ \(_)             | |                 |  ____|    (_)   | |
%  | |  | |_ _ __ ___  ___| |_ ___  _ __ _   _| |__  __  ___ ___| |_
%  | |  | | | '__/ _ \/ __| __/ _ \| '__| | | |  __| \ \/ / / __| __|
%  | |__| | | | |  __/ (__| || (_) | |  | |_| | |____ >  <| \__ \ |_
%  |_____/|_|_|  \___|\___|\__\___/|_|   \__, |______/_/\_\_|___/\__|
%                                         __/ |
%                                        |___/
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=DirectoryExist

% Basic directory test
DirectoryExist_test = pwd;
fcn_DebugTools_checkInputsToFunctions(DirectoryExist_test, 'DoesDirectoryExist');

%% Structure Types
% 
%    _____ _                   _                       
%   / ____| |                 | |                      
%  | (___ | |_ _ __ _   _  ___| |_ _   _ _ __ ___  ___ 
%   \___ \| __| '__| | | |/ __| __| | | | '__/ _ \/ __|
%   ____) | |_| |  | |_| | (__| |_| |_| | | |  __/\__ \
%  |_____/ \__|_|   \__,_|\___|\__|\__,_|_|  \___||___/
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=PathLibraryTypes

% Define a test structure and compare
template_structure = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});
structure_to_test = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});
fcn_DebugTools_checkInputsToFunctions(structure_to_test, 'likestructure',template_structure);

%% Test the structure type (success since more fields)
template_structure = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});
structure_to_test  = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{},'temp',{});
fcn_DebugTools_checkInputsToFunctions(structure_to_test, 'likestructure',template_structure);

%% Path Library Types
%   _____      _   _     _      _ _                       _______                    
%  |  __ \    | | | |   | |    (_) |                     |__   __|                   
%  | |__) |_ _| |_| |__ | |     _| |__  _ __ __ _ _ __ _   _| |_   _ _ __   ___  ___ 
%  |  ___/ _` | __| '_ \| |    | | '_ \| '__/ _` | '__| | | | | | | | '_ \ / _ \/ __|
%  | |  | (_| | |_| | | | |____| | |_) | | | (_| | |  | |_| | | |_| | |_) |  __/\__ \
%  |_|   \__,_|\__|_| |_|______|_|_.__/|_|  \__,_|_|   \__, |_|\__, | .__/ \___||___/
%                                                       __/ |   __/ | |              
%                                                      |___/   |___/|_|     
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=PathLibraryTypes
            

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
fcn_DebugTools_checkInputsToFunctions(station_test, 'station');



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
fcn_DebugTools_checkInputsToFunctions(station_test, 'stations');



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
fcn_DebugTools_checkInputsToFunctions(path_test, 'path');

%% path2or3D
% 
%               _   _     ___           ____  _____  
%              | | | |   |__ \         |___ \|  __ \ 
%   _ __   __ _| |_| |__    ) |___  _ __ __) | |  | |
%  | '_ \ / _` | __| '_ \  / // _ \| '__|__ <| |  | |
%  | |_) | (_| | |_| | | |/ /| (_) | |  ___) | |__| |
%  | .__/ \__,_|\__|_| |_|____\___/|_| |____/|_____/ 
%  | |                                               
%  |_|                                               

%% Test the path2or3D type (success)
clc;
path_test = [4 1; 2 1];
fcn_DebugTools_checkInputsToFunctions(path_test, 'path2or3D');

path_test = [4 1 0; 2 1 4];
fcn_DebugTools_checkInputsToFunctions(path_test, 'path2or3D');


%% Elevated Path
%   ______ _                 _           _   _____      _   _     
%  |  ____| |               | |         | | |  __ \    | | | |    
%  | |__  | | _____   ____ _| |_ ___  __| | | |__) |_ _| |_| |__  
%  |  __| | |/ _ \ \ / / _` | __/ _ \/ _` | |  ___/ _` | __| '_ \ 
%  | |____| |  __/\ V / (_| | ||  __/ (_| | | |  | (_| | |_| | | |
%  |______|_|\___| \_/ \__,_|\__\___|\__,_| |_|   \__,_|\__|_| |_|
%                                                                 
%                                                                 
%% Test the elevated path type (success)
clc;
elevated_path_test = [4 1 0.1; 2 1 0.2];
fcn_DebugTools_checkInputsToFunctions(elevated_path_test, 'elevated_path');

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
fcn_DebugTools_checkInputsToFunctions(paths_test, 'paths');



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
fcn_DebugTools_checkInputsToFunctions(test_traversal, 'traversal');



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
fcn_DebugTools_checkInputsToFunctions(test_traversals, 'traversals');





%% Fail conditions
if 1==0
    %% Bad string
    clc;
    column_of_numbers_test = [4 1];
    fcn_DebugTools_checkInputsToFunctions(column_of_numbers_test, 'dumb_text');
    
    %% 1column_of_numbers
    %
    %   __           _                               __                       _
    %  /_ |         | |                             / _|                     | |
    %   | | ___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ __  _   _ _ __ ___ | |__   ___ _ __ ___
    %   | |/ __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| '_ \| | | | '_ ` _ \| '_ \ / _ \ '__/ __|
    %   | | (_| (_) | | |_| | | | | | | | | || (_) | | | | | | |_| | | | | | | |_) |  __/ |  \__ \
    %   |_|\___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_| |_|\__,_|_| |_| |_|_.__/ \___|_|  |___/
    %                                     ______   ______
    %                                    |______| |______|
    % See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
    
    
    %% Test the column_of_numbers type (FAILURE because 1 x 2)
    clc;
    column_of_numbers_test = [4 1];
    fcn_DebugTools_checkInputsToFunctions(column_of_numbers_test, '1column_of_numbers');
    
    %% Test the column_of_numbers type (FAILURE because 3 long, not 2)
    clc;
    column_of_numbers_test = [4; 3; 2];
    fcn_DebugTools_checkInputsToFunctions(column_of_numbers_test, '1column_of_numbers',2);
    
    %% Test the column_of_numbers type (FAILURE because NaN)
    clc;
    column_of_numbers_test = [4; nan; 2];
    fcn_DebugTools_checkInputsToFunctions(column_of_numbers_test, '1column_of_numbers',3);
    
    
    %% positive_column_of_numbers
    %
    %                   _ _   _            __           _                               __                       _
    %                  (_) | (_)          /_ |         | |                             / _|                     | |
    %   _ __   ___  ___ _| |_ ___   _____  | | ___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ __  _   _ _ __ ___ | |__   ___ _ __ ___
    %  | '_ \ / _ \/ __| | __| \ \ / / _ \ | |/ __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| '_ \| | | | '_ ` _ \| '_ \ / _ \ '__/ __|
    %  | |_) | (_) \__ \ | |_| |\ V /  __/ | | (_| (_) | | |_| | | | | | | | | || (_) | | | | | | |_| | | | | | | |_) |  __/ |  \__ \
    %  | .__/ \___/|___/_|\__|_| \_/ \___| |_|\___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_| |_|\__,_|_| |_| |_|_.__/ \___|_|  |___/
    %  | |                             ______                                ______   ______
    %  |_|                            |______|                              |______| |______|
    % See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
    
    %% Test the positive_column_of_numbers type (FAILURE because 1 x 2)
    clc;
    positive_column_of_numbers_test = [4 1];
    fcn_DebugTools_checkInputsToFunctions(positive_column_of_numbers_test, 'positive_1column_of_numbers');
    
    %% Test the positive_column_of_numbers type (FAILURE because 3 long, not 2)
    clc;
    positive_column_of_numbers_test = [4; 3; 2];
    fcn_DebugTools_checkInputsToFunctions(positive_column_of_numbers_test, 'positive_1column_of_numbers',2);
    
    %% Test the positive_column_of_numbers type (FAILURE because NaN)
    clc;
    positive_column_of_numbers_test = [4; nan; 2];
    fcn_DebugTools_checkInputsToFunctions(positive_column_of_numbers_test, 'positive_1column_of_numbers',3);

    %% Test the positive_column_of_numbers type (FAILURE because negative value)
    clc;
    positive_column_of_numbers_test = [4; -1; 2];
    fcn_DebugTools_checkInputsToFunctions(positive_column_of_numbers_test, 'positive_1column_of_numbers',3);

    %% Test the positive_column_of_numbers type (FAILURE because zero value)
    clc;
    positive_column_of_numbers_test = [4; 1; 0];
    fcn_DebugTools_checkInputsToFunctions(positive_column_of_numbers_test, 'positive_1column_of_numbers',3);


    %% 2column_of_numbers
    %
    %   ___           _                               __                       _
    %  |__ \         | |                             / _|                     | |
    %     ) |___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ __  _   _ _ __ ___ | |__   ___ _ __ ___
    %    / // __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| '_ \| | | | '_ ` _ \| '_ \ / _ \ '__/ __|
    %   / /| (_| (_) | | |_| | | | | | | | | || (_) | | | | | | |_| | | | | | | |_) |  __/ |  \__ \
    %  |____\___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_| |_|\__,_|_| |_| |_|_.__/ \___|_|  |___/
    %                                      ______   ______
    %                                     |______| |______|
    %
    % See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
    
    %% Test the column_of_numbers type (FAILURE because 1 x 1)
    clc;
    Twocolumn_of_numbers_test = 4;
    fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers');
    
    %% Test the column_of_numbers type (FAILURE because 3 long, not 2)
    clc;
    Twocolumn_of_numbers_test = [4 1; 3 1; 2 2];
    fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',2);
        
    %% Test the column_of_numbers type (FAILURE because 4 columns)
    clc;
    Twocolumn_of_numbers_test = [4 1 1 4];
    fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers');
       
    %% Test the column_of_numbers type (FAILURE because NaN)
    clc;
    Twocolumn_of_numbers_test = [4 1; nan 1; 2 0];
    fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',3);
    
    
    %% Minimum length is 4 or greater
    clc;
    Twocolumn_of_numbers_test = [4 1; 3 9; 2 7];
    fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',[4 5]);
    
    %% Maximum length is 2 or less
    clc;
    Twocolumn_of_numbers_test = [4 1; 3 9; 2 7];
    fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',[2 1]);
    
    %% Maximum length is 2
    clc;
    Twocolumn_of_numbers_test = [4 1; 3 9; 2 7];
    fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',[2 2]);
    
    %% 2or3column_of_numbers
    %
    %   ___           ____            _                               __                       _
    %  |__ \         |___ \          | |                             / _|                     | |
    %     ) |___  _ __ __) | ___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ __  _   _ _ __ ___ | |__   ___ _ __ ___
    %    / // _ \| '__|__ < / __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| '_ \| | | | '_ ` _ \| '_ \ / _ \ '__/ __|
    %   / /| (_) | |  ___) | (_| (_) | | |_| | | | | | | | | || (_) | | | | | | |_| | | | | | | |_) |  __/ |  \__ \
    %  |____\___/|_| |____/ \___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_| |_|\__,_|_| |_| |_|_.__/ \___|_|  |___/
    %                                                      ______   ______
    %                                                     |______| |______|
    % See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
    
    

    %% Test the column_of_numbers type (FAILURE because 1 x 1)
    clc;
    TwoOrThreeColumn_of_numbers_test = 4;
    fcn_DebugTools_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers');

    %% Test the column_of_numbers type (FAILURE because 1 x 4)
    clc;
    TwoOrThreeColumn_of_numbers_test = [4 1 1 1];
    fcn_DebugTools_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers');

    %% Test the column_of_numbers type (FAILURE because 3 long, not 2)
    clc;
    TwoOrThreeColumn_of_numbers_test = [4 1; 3 1; 2 2];
    fcn_DebugTools_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',2);
       
    %% Test the column_of_numbers type (FAILURE because NaN)
    clc;
    TwoOrThreeColumn_of_numbers_test = [4 1; nan 1; 2 0];
    fcn_DebugTools_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',3);
    
    
    %% Minimum length is 4 or greater
    clc;
    TwoOrThreeColumn_of_numbers_test = [4 1; 3 9; 2 7];
    fcn_DebugTools_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[4 5]);
    
    %% Maximum length is 2 or less
    clc;
    TwoOrThreeColumn_of_numbers_test = [4 1; 3 9; 2 7];
    fcn_DebugTools_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[2 1]);
    
    %% Maximum length is 2
    clc;
    TwoOrThreeColumn_of_numbers_test = [4 1; 3 9; 2 7];
    fcn_DebugTools_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[2 2]);


    %% 2or3column_of_numbers
    %   _   _            __  __           _                               __                       _
    %  | \ | |          |  \/  |         | |                             / _|                     | |
    %  |  \| | ___  _ __| \  / | ___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ __  _   _ _ __ ___ | |__   ___ _ __ ___
    %  | . ` |/ _ \| '__| |\/| |/ __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| '_ \| | | | '_ ` _ \| '_ \ / _ \ '__/ __|
    %  | |\  | (_) | |  | |  | | (_| (_) | | |_| | | | | | | | | || (_) | | | | | | |_| | | | | | | |_) |  __/ |  \__ \
    %  |_| \_|\___/|_|  |_|  |_|\___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_| |_|\__,_|_| |_| |_|_.__/ \___|_|  |___/
    %                                                          ______   ______
    %                                                         |______| |______|
    % See: http://patorjk.com/software/taag/#p=display&f=Big&t=NorMcolumn_of_numbers
    %% Test the N to M column_of_numbers type (fail as only 1 column)
    clc;
    NorMColumn_of_numbers_test = 4;
    fcn_DebugTools_checkInputsToFunctions(NorMColumn_of_numbers_test, '2or4column_of_numbers');

    %% Test the N to M column_of_numbers type (fail as has 5 columns)
    clc;
    NorMColumn_of_numbers_test = [4 4 2 5 1];    
    fcn_DebugTools_checkInputsToFunctions(NorMColumn_of_numbers_test, '2or4column_of_numbers');



    %% 2column_of_integers
    %
    %
    %   ___           _                               __ _       _
    %  |__ \         | |                             / _(_)     | |
    %     ) |___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ _ __ | |_ ___  __ _  ___ _ __ ___
    %    / // __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| | '_ \| __/ _ \/ _` |/ _ \ '__/ __|
    %   / /| (_| (_) | | |_| | | | | | | | | || (_) | | | | | | | ||  __/ (_| |  __/ |  \__ \
    %  |____\___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_|_| |_|\__\___|\__, |\___|_|  |___/
    %                                      ______   ______                __/ |
    %                                     |______| |______|              |___/
    % See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
    %% Test the 2column_of_integers type (FAILURE because not integers)
    clc;
    Twocolumn_of_integers_test = [4 3.2];
    fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers');
    
    %% Test the 2column_of_integers type (FAILURE because 1 x 1)
    clc;
    Twocolumn_of_integers_test = 4;
    fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers');
    
    %% Test the 2column_of_integers type (FAILURE because 3 long, not 2)
    clc;
    Twocolumn_of_integers_test = [4 1; 3 1; 2 2];
    fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers',2);
       
    %% Test the 2column_of_integers type (FAILURE because NaN)
    clc;
    Twocolumn_of_integers_test = [4 1; nan 1; 2 0];
    fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers',3);
    
    
    %% Minimum length is 4 or greater
    clc;
    Twocolumn_of_integers_test = [4 1; 3 9; 2 7];
    fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers',[4 5]);
    
    %% Maximum length is 2 or less
    clc;
    Twocolumn_of_integers_test = [4 1; 3 9; 2 7];
    fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers',[2 1]);
    
    %% Maximum length is 2
    clc;
    Twocolumn_of_integers_test = [4 1; 3 9; 2 7];
    fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers',[2 2]);

    %% Char Types
    %    _____ _
    %   / ____| |
    %  | |    | |__   __ _ _ __ ___
    %  | |    | '_ \ / _` | '__/ __|
    %  | |____| | | | (_| | |  \__ \
    %   \_____|_| |_|\__,_|_|  |___/
    % See: http://patorjk.com/software/taag/#p=display&f=Big&t=Chars

    %% Basic character test (fail - not character)
    clc
    Chars_test = 1;
    fcn_DebugTools_checkInputsToFunctions(Chars_test, '_of_chars');

    %% Basic character test (fail - not character)
    clc
    Chars_test = "abcdefg";
    fcn_DebugTools_checkInputsToFunctions(Chars_test, '_of_chars');


    %% String Types
    %    _____ _        _
    %   / ____| |      (_)
    %  | (___ | |_ _ __ _ _ __   __ _ ___
    %   \___ \| __| '__| | '_ \ / _` / __|
    %   ____) | |_| |  | | | | | (_| \__ \
    %  |_____/ \__|_|  |_|_| |_|\__, |___/
    %                            __/ |
    %                           |___/
    % See: http://patorjk.com/software/taag/#p=display&f=Big&t=Strings

    %% Basic string test (pass)
    clc
    String_test = "abcdefg";
    fcn_DebugTools_checkInputsToFunctions(String_test, '_of_strings');

    %% Basic string test (fail - not a string)
    clc
    String_test = 'abcdefg';
    fcn_DebugTools_checkInputsToFunctions(String_test, '_of_strings');

    %% DoesFileExist
    %
    %   ______ _ _      ______      _     _
    %  |  ____(_) |    |  ____|    (_)   | |
    %  | |__   _| | ___| |__  __  ___ ___| |_
    %  |  __| | | |/ _ \  __| \ \/ / / __| __|
    %  | |    | | |  __/ |____ >  <| \__ \ |_
    %  |_|    |_|_|\___|______/_/\_\_|___/\__|
    %
    %
    % See: http://patorjk.com/software/taag/#p=display&f=Big&t=FileExist

    % Basic string test
    FileExist_test = 'junkjunkjunk.m';
    fcn_DebugTools_checkInputsToFunctions(FileExist_test, 'DoesFileExist');

    %% DoesDirectoryExist
    %
    %   _____  _               _                   ______      _     _
    %  |  __ \(_)             | |                 |  ____|    (_)   | |
    %  | |  | |_ _ __ ___  ___| |_ ___  _ __ _   _| |__  __  ___ ___| |_
    %  | |  | | | '__/ _ \/ __| __/ _ \| '__| | | |  __| \ \/ / / __| __|
    %  | |__| | | | |  __/ (__| || (_) | |  | |_| | |____ >  <| \__ \ |_
    %  |_____/|_|_|  \___|\___|\__\___/|_|   \__, |______/_/\_\_|___/\__|
    %                                         __/ |
    %                                        |___/
    %
    % See: http://patorjk.com/software/taag/#p=display&f=Big&t=FileExist

    % Basic string test
    DirectoryExist_test = 'junkjunkjunk';
    fcn_DebugTools_checkInputsToFunctions(DirectoryExist_test, 'DoesDirectoryExist');
    
    %% Station
    %    _____ _        _   _
    %   / ____| |      | | (_)
    %  | (___ | |_ __ _| |_ _  ___  _ __
    %   \___ \| __/ _` | __| |/ _ \| '_ \
    %   ____) | || (_| | |_| | (_) | | | |
    %  |_____/ \__\__,_|\__|_|\___/|_| |_|
    %
    
    
    %% Test the station type (success)
    clc;
    station_test = 4;
    fcn_DebugTools_checkInputsToFunctions(station_test, 'station');
    
    
    %% Test the station type (fail since non-numeric)
    clc;
    station_test = 'junk';
    fcn_DebugTools_checkInputsToFunctions(station_test, 'station');
    
    
    %% Test the station type (fail since not 1 column)
    clc;
    station_test = [4 0];
    fcn_DebugTools_checkInputsToFunctions(station_test, 'station');
    
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
    clc;
    station_test = [4; 2];
    fcn_DebugTools_checkInputsToFunctions(station_test, 'stations');
    
    
    %% Test the stations type (fail since only one row)
    clc;
    station_test = 4;
    fcn_DebugTools_checkInputsToFunctions(station_test, 'stations');
    
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
    fcn_DebugTools_checkInputsToFunctions(path_test, 'path');
    
    
    %% Test the path type (fail since only one column)
    clc
    path_test = [4; 2];
    fcn_DebugTools_checkInputsToFunctions(path_test, 'path');
    
    %% Test the path type (fail since only one row)
    clc
    path_test = [4 2];
    fcn_DebugTools_checkInputsToFunctions(path_test, 'path');
    
    %% path2or3D
    %
    %               _   _     ___           ____  _____
    %              | | | |   |__ \         |___ \|  __ \
    %   _ __   __ _| |_| |__    ) |___  _ __ __) | |  | |
    %  | '_ \ / _` | __| '_ \  / // _ \| '__|__ <| |  | |
    %  | |_) | (_| | |_| | | |/ /| (_) | |  ___) | |__| |
    %  | .__/ \__,_|\__|_| |_|____\___/|_| |____/|_____/
    %  | |
    %  |_|
    
    %% Test the path2or3D type (success)
    clc;
    path_test = [4 1; 2 1];
    fcn_DebugTools_checkInputsToFunctions(path_test, 'path2or3D');
    
    path_test = [4 1 0; 2 1 4];
    fcn_DebugTools_checkInputsToFunctions(path_test, 'path2or3D');
    
    %% Test the path type (fail since only one column)
    clc
    path_test = [4; 2];
    fcn_DebugTools_checkInputsToFunctions(path_test, 'path2or3D');
    
    %% Test the path type (fail since only one row)
    clc
    path_test = [4 2 1];
    fcn_DebugTools_checkInputsToFunctions(path_test, 'path2or3D');
    
    
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
    fcn_DebugTools_checkInputsToFunctions(paths_test, 'paths');
    
    
    %% Test the paths type (fail since only one column)
    clc
    paths_test = [4; 2];
    fcn_DebugTools_checkInputsToFunctions(paths_test, 'paths');
    
    %% Test the paths type (fail since only two rows)
    clc
    paths_test = [4 2; 0 0];
    fcn_DebugTools_checkInputsToFunctions(paths_test, 'paths');
    
    %% Structure Types
    %
    %    _____ _                   _
    %   / ____| |                 | |
    %  | (___ | |_ _ __ _   _  ___| |_ _   _ _ __ ___  ___
    %   \___ \| __| '__| | | |/ __| __| | | | '__/ _ \/ __|
    %   ____) | |_| |  | |_| | (__| |_| |_| | | |  __/\__ \
    %  |_____/ \__|_|   \__,_|\___|\__|\__,_|_|  \___||___/
    %
    % See: http://patorjk.com/software/taag/#p=display&f=Big&t=PathLibraryTypes
    
    %% Test the structure type (success)
    clc;
    template_structure = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});
    structure_to_test = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});
    fcn_DebugTools_checkInputsToFunctions(structure_to_test, 'likestructure',template_structure);

    
    %% Test the structure type (success since more fields)
    clc;
    template_structure = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});
    structure_to_test  = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{},'temp',{});
    fcn_DebugTools_checkInputsToFunctions(structure_to_test, 'likestructure',template_structure);

    %% Test the structure type (fail since not a structure)
    clc;
    template_structure = [1 2 3];
    structure_to_test = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});
    fcn_DebugTools_checkInputsToFunctions(structure_to_test, 'likestructure',template_structure);

    %% Test the structure type (fail since fields are wrong)
    clc;
    template_structure = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});
    structure_to_test = struct('A',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});
    fcn_DebugTools_checkInputsToFunctions(structure_to_test, 'likestructure',template_structure);

    %% Test the structure type (fail since less fields)
    clc;
    template_structure = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});
    structure_to_test  = struct('color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});
    fcn_DebugTools_checkInputsToFunctions(structure_to_test, 'likestructure',template_structure);

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
    clc;
    % Fill in sample paths (as a starter)
    paths = fcn_Path_fillSamplePaths;
    
    % Convert paths to traversal structures
    for i_Path = 1:length(paths)
        test_traversal = fcn_Path_convertPathToTraversalStructure(paths{i_Path});
        all_traversals.traversal{i_Path} = test_traversal;
    end
    fcn_DebugTools_checkInputsToFunctions(test_traversal, 'traversal');
    
    
    %% Test the traversal type (fail since field is missing)
    clc
    clear test_traversal_bad
    test_traversal_bad = 3;
    fcn_DebugTools_checkInputsToFunctions(test_traversal_bad, 'traversal');
    
    %% Test the traversal type (fail since Z field is missing)
    clc;
    clear test_traversal_bad
    test_traversal_bad.X = 'junk';
    test_traversal_bad.Y = 'junk';
    test_traversal_bad.Station = 'junk';
    
    fcn_DebugTools_checkInputsToFunctions(test_traversal_bad, 'traversal');
    
    %% Test the traversal type (fail since field is not numeric)
    clc;
    clear test_traversal_bad
    test_traversal_bad.X = 'junk';
    test_traversal_bad.Y = 'junk';
    test_traversal_bad.Z = 'junk';
    test_traversal_bad.Station = 'junk';
    
    fcn_DebugTools_checkInputsToFunctions(test_traversal_bad, 'traversal');
    
    %% Test the traversal type (fail since fields are not columns)
    clc;
    clear test_traversal_bad
    test_traversal_bad.X = eye(2);
    test_traversal_bad.Y = eye(2);
    test_traversal_bad.Z = eye(2);
    test_traversal_bad.Station = eye(2);
    
    fcn_DebugTools_checkInputsToFunctions(test_traversal_bad, 'traversal');
    
    %% Test the traversal type (fail since fields have different lengths)
    clc;
    clear test_traversal_bad
    test_traversal_bad.X = [1; 2; 3];
    test_traversal_bad.Y = [1; 2; 3];
    test_traversal_bad.Z = [1; 2; 3];
    test_traversal_bad.Station = [1; 2];
    
    fcn_DebugTools_checkInputsToFunctions(test_traversal_bad, 'traversal');
    
    %% Test the traversal type (fail since Station field is not strictly increasing)
    clc;
    clear test_traversal_bad
    test_traversal_bad.X = [1; 2; 3];
    test_traversal_bad.Y = [1; 2; 3];
    test_traversal_bad.Z = [1; 2; 3];
    test_traversal_bad.Station = [1; 2; 2];
    
    fcn_DebugTools_checkInputsToFunctions(test_traversal_bad, 'traversal');
    
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
    fcn_DebugTools_checkInputsToFunctions(test_traversals, 'traversals');
    
    
    %% Test the traversals type (fail since field is missing)
    clc
    clear test_traversals_bad
    test_traversals_bad = 3;
    fcn_DebugTools_checkInputsToFunctions(test_traversals_bad, 'traversals');
    
    %% Test the traversals type (fail since second index is bad)
    clc
    clear test_traversals_bad
    test_traversals_bad.traversal{1} = test_traversal;
    
    clear test_traversal_bad
    test_traversal_bad.X = [1; 2; 3];
    test_traversal_bad.Y = [1; 2; 3];
    test_traversal_bad.Z = [1; 2; 3];
    test_traversal_bad.Station = [1; 2; 2];
    
    test_traversals_bad.traversal{2} = test_traversal_bad;
    fcn_DebugTools_checkInputsToFunctions(test_traversals_bad, 'traversals');


end

