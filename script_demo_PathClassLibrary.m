% This is a demonstration script to show the primary functionality of the
% path class library.

%% Set up workspace
if ~exist('flag_paths_were_added_already','var')
    
    clc
    close all
    
    % add necessary directories for functions recursively
    addpath(genpath([pwd '/Functions']))
    
    % % add necessary directories for data
    % addpath([pwd '/Data'])  % This is where sample maps are stored
    
    % add necessary directories for Utilities to the path
    if(exist([pwd '/Utilities'],'dir'))
        addpath(genpath([pwd '/Utilities']))  % This is where GPS utilities are stored
    else
        error('No Utilities directory exists to be added to the path. Please create one (see README.md) and run again.');
    end
    
    % set a flag so we do not have to do this again
    flag_paths_were_added_already = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See http://patorjk.com/software/taag/#p=display&f=Big&t=Getting%20Started
% 
%    _____      _   _   _                _____ _             _           _ 
%   / ____|    | | | | (_)              / ____| |           | |         | |
%  | |  __  ___| |_| |_ _ _ __   __ _  | (___ | |_ __ _ _ __| |_ ___  __| |
%  | | |_ |/ _ \ __| __| | '_ \ / _` |  \___ \| __/ _` | '__| __/ _ \/ _` |
%  | |__| |  __/ |_| |_| | | | | (_| |  ____) | || (_| | |  | ||  __/ (_| |
%   \_____|\___|\__|\__|_|_| |_|\__, | |_____/ \__\__,_|_|   \__\___|\__,_|
%                                __/ |                                     
%                               |___/                                      
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Show how input arguments are checked, fcn_Path_checkInputsToFunctions
% TO-DO - move debug tools into utilities and remove the checkinputs capability out of Path class
path_test = [4 1; 2 1];
fcn_Path_checkInputsToFunctions(path_test, 'path');

%% Show how to load sample paths, fcn_Path_fillSamplePaths  
% Call the function to fill in an array of "path" type
paths_array = fcn_Path_fillSamplePaths;

% We can even save one of these as a single "path"
single_path = paths_array{1};

%% Show how to calculate the relative angle changes between path segments, fcn_Path_calcDiffAnglesBetweenPathSegments
% Avoids the use of atan function since this breaks near -180 degree point

% Pick first path as reference_traversal structure
paths_array = fcn_Path_fillSamplePaths;
paths_to_check = paths_array{1};
fig_num = 11111; 
diff_angles = fcn_Path_calcDiffAnglesBetweenPathSegments(paths_to_check,fig_num);


%% Show how to calculate the yaw angles along a path, fcn_Path_calcYawFromPathSegments

% Basic call with one path
fig_num = 22222;
yaw_angles = fcn_Path_calcYawFromPathSegments(path_to_check,fig_num); %#ok<*NASGU>


% Multiple paths
for i_path = 1:length(paths_array)
    % Pick first path as reference_traversal structure
    path_to_check = paths_array{i_path};
    
    % Pick first path as reference_traversal structure
    traversal_to_check = fcn_Path_convertPathToTraversalStructure(paths_array{i_path});
    all_traversals.traversal{i_path} = traversal_to_check;
    
    fig_num = 33333;
    yaw_angles = fcn_Path_calcYawFromPathSegments(path_to_check,fig_num);
end

% Plot the results? (Note: they are plotted below as well)
if 1==1
    fig_num = 333331;
    fcn_Path_plotTraversalsYaw(all_traversals,fig_num);
    fig_num = 333332;
    fcn_Path_plotTraversalsXY(all_traversals,fig_num);
end

%% Show how to convert paths to traversals, fcn_Path_convertPathToTraversalStructure
fig_num = 44444;
traversal = fcn_Path_convertPathToTraversalStructure(single_path,fig_num);
simple_example.traversal{1} = traversal;

fcn_Path_plotTraversalsXY(simple_example,fig_num);
xlabel('X [m]');
ylabel('Y [m]');

%% Show how to plot a traversal, fcn_Path_convertPathToTraversalStructure
fig_num = 55555;

fcn_Path_plotTraversalsXY(simple_example,fig_num);
xlabel('X [m]');
ylabel('Y [m]');

%% Show how to plot the yaw angle of traversal, fcn_Path_plotTraversalsYaw

% Fill in some dummy data
paths_array = fcn_Path_fillSamplePaths;
 

% Convert paths into traversals
for i_traveral = 1:length(paths_array)
    traversal = fcn_Path_convertPathToTraversalStructure(paths_array{i_traveral});
    data.traversal{i_traveral} = traversal;
end

fig_num = 666661;
fcn_Path_plotTraversalsXY(data,fig_num);


% Next, specify the figure number to show that it will NOT auto-label the
% axes if figure is already given and it puts the plots into this figure.
fig_num = 666662;
fcn_Path_plotTraversalsYaw(data,fig_num);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See http://patorjk.com/software/taag/#p=display&f=Big&t=Basic%20Path%20Operations
%   ____            _        _____      _   _        ____                       _   _                 
%  |  _ \          (_)      |  __ \    | | | |      / __ \                     | | (_)                
%  | |_) | __ _ ___ _  ___  | |__) |_ _| |_| |__   | |  | |_ __   ___ _ __ __ _| |_ _  ___  _ __  ___ 
%  |  _ < / _` / __| |/ __| |  ___/ _` | __| '_ \  | |  | | '_ \ / _ \ '__/ _` | __| |/ _ \| '_ \/ __|
%  | |_) | (_| \__ \ | (__  | |  | (_| | |_| | | | | |__| | |_) |  __/ | | (_| | |_| | (_) | | | \__ \
%  |____/ \__,_|___/_|\___| |_|   \__,_|\__|_| |_|  \____/| .__/ \___|_|  \__,_|\__|_|\___/|_| |_|___/
%                                                         | |                                         
%                                                         |_|                                         
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Show how to find intersection of a path onto another path, fcn_Path_findProjectionHitOntoPath
% This is the core function of the entire library

% Simple test 3 - multiple intersections
fprintf(1,'Multiple intersections result: \n');
path = [0 10; 10 10; 0 6; 10 6; 0 2];
sensor_vector_start = [0 0]; 
sensor_vector_end   = [5 12];
fig_debugging = 2343;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);

URHERE


