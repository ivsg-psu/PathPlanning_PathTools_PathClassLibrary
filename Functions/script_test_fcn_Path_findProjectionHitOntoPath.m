% script_test_fcn_Path_findProjectionHitOntoPath
% This is a script to exercise the function: fcn_Path_findProjectionHitOntoPath.m
% This function was written on 2020_11_10 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Modification history:
% 2020_12_31 - S. Brennan
% -- updated for new argument list
%
% 2021_01_08 - S. Brennan
% -- updated comments
%
% 2020_01_09 - S. Brennan
% -- added more comments during clean-up
%
% 2021_01_23 - S. Brennan
% -- Added flag_search_type = 2 option, to allow multiple cross points to
% be returned
% -- Fixed bug with partially overlapping vectors not returning a result
% -- Added path segment output so that we can ID which segment was hit
%
% 2021_01_24 - S. Brennan
% -- Fixed bug with overlapping colinear where two path segments identified
% when there is only one
%
% 2021_01_24 - S. Brennan
% -- Added assertions to force checking
% -- Added test cases for flag=3 flag=4 options
%
% 2024_05_15 - Aneesh Batchu
% -- Added a BUG case
%
% 2024_05_15 - S. Brennan
% -- Organized sections
% -- Added fast mode tests
% -- Added full assertion tests on output variables
% -- Added figure open/close assertions
% -- Added automated assertion testing

close all

%%%%%%%%%%%%
% FIGURE NUMBERING:
% FSIXXX
%
% F is First figure number, starting with:
% 1: demonstration cases
% 2: single point intersection cases
% 3: non intersection cases
% 4: infinite intersection cases
% 5: multi-hit cases
% 6: multi-hit overlapping cases
% 7: multi-hit multi-path cases
% 8: fast mode cases
% 9: known bug cases
%
% S is second figure number, flag_search_type:
% 0: first intersection if there is any overlap
% 1: first intersection of sensor vector to the path, for ANY directional extension of the sensor.
% 2: same as flag_search_type 0, except all intersections if there is overlap, and only first and last contact points if infinte
% 3: first intersection of sensor vector to path, for ANY directional extension of the path. Note: this is opposite of flag 1
% 4: first intersection of sensor vector to path, for ANY directional extension of the path OR the sensor. Note: this is the combinations of flags 1 and 3
%
% XXX: 2nd to 5th number: a counter that counts up through the cases in this
% section.
%
% Example:
% 24006 plots a single point intersection, flag type 4, 6th test case


%% Demonstration cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% _____                                 _             _   _                _____
% |  __ \                               | |           | | (_)              / ____|
% | |  | | ___ _ __ ___   ___  _ __  ___| |_ _ __ __ _| |_ _  ___  _ __   | |     __ _ ___  ___  ___
% | |  | |/ _ \ '_ ` _ \ / _ \| '_ \/ __| __| '__/ _` | __| |/ _ \| '_ \  | |    / _` / __|/ _ \/ __|
% | |__| |  __/ | | | | | (_) | | | \__ \ |_| | | (_| | |_| | (_) | | | | | |___| (_| \__ \  __/\__ \
% |_____/ \___|_| |_| |_|\___/|_| |_|___/\__|_|  \__,_|\__|_|\___/|_| |_|  \_____\__,_|___/\___||___/
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Demonstration%20Cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All demonstration case figures start with the number 1

%% Demo of standard format checking
fig_num = 10001;
fprintf(1,'Figure: %.0f :Demo of short format checking\n',fig_num);
figure(fig_num); clf;
plotting.FigureExpected = 1;

clear inputs
inputs.fig_num = fig_num;
inputs.path = [0 10; 10 10];
inputs.sensor_vector_start = [5 0];
inputs.sensor_vector_end   = [5 15];
inputs.flag_search_type = 0;

% SHORT format checking
clear expected
expected.distance = 10;
expected.location = [5 10];
expected.path_segment = 1;
expected.t = 0.5;
expected.u = 0.6667;

actual = struct;
[actual.distance, actual.location, actual.path_segment, actual.t, actual.u] = ...
    fcn_Path_findProjectionHitOntoPath(...
    inputs.path,inputs.sensor_vector_start,inputs.sensor_vector_end,...
    inputs.flag_search_type,inputs.fig_num);

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)
%fcn_INTERNAL_printResults(actual.distance,actual.location);

%% Demo of long format checking
fig_num = 10002;
fprintf(1,'Figure: %.0f :Demo of long format checking\n',fig_num);
figure(fig_num); clf;
plotting.FigureExpected = 1;

Nsolutions = 1;

clear inputs
inputs.fig_num = fig_num;
inputs.path = [0 10; 10 10];
inputs.sensor_vector_start = [5 0];
inputs.sensor_vector_end   = [5 15];
inputs.flag_search_type = 0;

NpathSegments = length(inputs.path(:,1))-1;

% LONG format checking
clear expected
expected.distance.TypeString      = 'numeric'; % See 'help isa' for a listing
expected.distance.Size            = [Nsolutions 1];
expected.distance.Value           = 10;
expected.distance.ValueDigits     = 4; % How many digits must match

expected.location.TypeString      = 'numeric'; % See 'help isa' for a listing
expected.location.Size            = [Nsolutions 2];
expected.location.Value           = [5 10];
expected.location.ValueDigits     = 4; % How many digits must match

expected.path_segment.TypeString  = 'numeric'; % See 'help isa' for a listing
expected.path_segment.Size        = [Nsolutions 1];
expected.path_segment.Value       = 1;
expected.path_segment.ValueDigits = 0; % How many digits must match

expected.t.TypeString  = 'numeric'; % See 'help isa' for a listing
expected.t.Size        = [NpathSegments 1]; % Expecting 1 solution
expected.t.Value       = 0.5;
expected.t.ValueDigits = 4; % How many digits must match

expected.u.TypeString  = 'numeric'; % See 'help isa' for a listing
expected.u.Size        = [NpathSegments 1]; % Expecting 1 solution
expected.u.Value       = 0.6667;
expected.u.ValueDigits = 4; % How many digits must match

actual = struct;
[actual.distance, actual.location, actual.path_segment, actual.t, actual.u] = ...
    fcn_Path_findProjectionHitOntoPath(...
    inputs.path,inputs.sensor_vector_start,inputs.sensor_vector_end,...
    inputs.flag_search_type,inputs.fig_num);


fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)
%fcn_INTERNAL_printResults(actual.distance,actual.location);


%% Single point intersection cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____ _             _        _____      _       _     _____       _                          _   _
%  / ____(_)           | |      |  __ \    (_)     | |   |_   _|     | |                        | | (_)
% | (___  _ _ __   __ _| | ___  | |__) |__  _ _ __ | |_    | |  _ __ | |_ ___ _ __ ___  ___  ___| |_ _  ___  _ __  ___
%  \___ \| | '_ \ / _` | |/ _ \ |  ___/ _ \| | '_ \| __|   | | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __| |/ _ \| '_ \/ __|
%  ____) | | | | | (_| | |  __/ | |  | (_) | | | | | |_   _| |_| | | | ||  __/ |  \__ \  __/ (__| |_| | (_) | | | \__ \
% |_____/|_|_| |_|\__, |_|\___| |_|   \___/|_|_| |_|\__| |_____|_| |_|\__\___|_|  |___/\___|\___|\__|_|\___/|_| |_|___/
%                  __/ |
%                 |___/
%
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Single%20Point%20Intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All single point intersection figures start with the number 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____ _             _                    _       _     _       _                          _   _
%  / ____(_)           | |                  (_)     | |   (_)     | |                        | | (_)
% | (___  _ _ __   __ _| | ___   _ __   ___  _ _ __ | |_   _ _ __ | |_ ___ _ __ ___  ___  ___| |_ _  ___  _ __
%  \___ \| | '_ \ / _` | |/ _ \ | '_ \ / _ \| | '_ \| __| | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __| |/ _ \| '_ \
%  ____) | | | | | (_| | |  __/ | |_) | (_) | | | | | |_  | | | | | ||  __/ |  \__ \  __/ (__| |_| | (_) | | | |
% |_____/|_|_| |_|\__, |_|\___| | .__/ \___/|_|_| |_|\__| |_|_| |_|\__\___|_|  |___/\___|\___|\__|_|\___/|_| |_|
%                  __/ |        | |
%                 |___/         |_|
%   __ _                         ___
%  / _| |              ______   / _ \
% | |_| | __ _  __ _  |______| | | | |
% |  _| |/ _` |/ _` |  ______  | | | |
% | | | | (_| | (_| | |______| | |_| |
% |_| |_|\__,_|\__, |           \___/
%               __/ |
%              |___/
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Single%20point%20intersection%20%0Aflag%20%3D%200
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;

%% Single point intersection (2), flag (0), test 1 - a simple intersection

fig_num = 20001;
fprintf(1,'Figure: %.0f :Single point intersection (2), flag (0) - a simple intersection  \n',fig_num);
figure(fig_num); clf;
plotting.FigureExpected = 1;

clear inputs
inputs.fig_num = fig_num;
inputs.path = [0 10; 10 10];
inputs.sensor_vector_start = [5 0];
inputs.sensor_vector_end   = [5 15];
inputs.flag_search_type = 0;

% SHORT format checking
clear expected
expected.distance = 10;
expected.location = [5 10];
expected.path_segment = 1;
expected.t = 0.5;
expected.u = 0.6667;

actual = struct;
[actual.distance, actual.location, actual.path_segment, actual.t, actual.u] = ...
    fcn_Path_findProjectionHitOntoPath(...
    inputs.path,inputs.sensor_vector_start,inputs.sensor_vector_end,...
    inputs.flag_search_type,inputs.fig_num);

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)
%fcn_INTERNAL_printResults(actual.distance,actual.location);

%% Single point intersection (2), flag (0), test 2 - intersection through an acute vertex

fig_num = 20002;
fprintf(1,'Figure: %.0f : Single point intersection (2), flag (0) - intersection through an acute vertex  \n',fig_num);
figure(fig_num); clf;
plotting.FigureExpected = 1;

clear inputs
inputs.fig_num = fig_num;
inputs.path = [0 10; 5 10; 15 15];
inputs.sensor_vector_start = [5 0];
inputs.sensor_vector_end   = [5 15];
inputs.flag_search_type = 0;

% SHORT format checking
clear expected
expected.distance = 10;
expected.location = [5 10];
expected.path_segment = 1;
expected.t = 1;
expected.u = 0.6667;

actual = struct;
[actual.distance, actual.location, actual.path_segment, actual.t, actual.u] = ...
    fcn_Path_findProjectionHitOntoPath(...
    inputs.path,inputs.sensor_vector_start,inputs.sensor_vector_end,...
    inputs.flag_search_type,inputs.fig_num);

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)
%fcn_INTERNAL_printResults(actual.distance,actual.location);

%% Single point intersection (2), flag (0), test 3 - intersection through an obtuse vertex

fig_num = 20003;
fprintf(1,'Figure: %.0f : Single point intersection (2), flag (0) - intersection through an obtuse vertex  \n',fig_num);
figure(fig_num); clf;
plotting.FigureExpected = 1;

clear inputs
inputs.fig_num = fig_num;
inputs.path = [0 10; 5 10; 4 15];
inputs.sensor_vector_start = [5 0];
inputs.sensor_vector_end   = [5 15];
inputs.flag_search_type = 0;

% SHORT format checking
clear expected
expected.distance = 10;
expected.location = [5 10];
expected.path_segment = 1;
expected.t = 1;
expected.u = 0.6667;

actual = struct;
[actual.distance, actual.location, actual.path_segment, actual.t, actual.u] = ...
    fcn_Path_findProjectionHitOntoPath(...
    inputs.path,inputs.sensor_vector_start,inputs.sensor_vector_end,...
    inputs.flag_search_type,inputs.fig_num);

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)
%fcn_INTERNAL_printResults(actual.distance,actual.location);

%% Single point intersection (2), flag (0), test 4 -  intersection at start of sensor

fig_num = 20004;
fprintf(1,'Figure: %.0f : Single point intersection (2), flag (0) - intersection at start of sensor  \n',fig_num);
figure(fig_num); clf;
plotting.FigureExpected = 1;

clear inputs
inputs.fig_num = fig_num;
inputs.path = [0 10; 10 10];
inputs.sensor_vector_start = [5 10];
inputs.sensor_vector_end   = [5 15];
inputs.flag_search_type = 0;

% SHORT format checking
clear expected
expected.distance = 0;
expected.location = [5 10];
expected.path_segment = 1;
expected.t = 0.5;
expected.u = 0;

actual = struct;
[actual.distance, actual.location, actual.path_segment, actual.t, actual.u] = ...
    fcn_Path_findProjectionHitOntoPath(...
    inputs.path,inputs.sensor_vector_start,inputs.sensor_vector_end,...
    inputs.flag_search_type,inputs.fig_num);

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)
%fcn_INTERNAL_printResults(actual.distance,actual.location);

%% Single point intersection (2), flag (0), test 5 -  intersection at end of sensor

fig_num = 20005;
fprintf(1,'Figure: %.0f : Single point intersection (2), flag (0) - intersection at end of sensor  \n',fig_num);
figure(fig_num); clf;
plotting.FigureExpected = 1;

clear inputs
inputs.fig_num = fig_num;
inputs.path = [0 10; 10 10];
inputs.sensor_vector_start = [5 5];
inputs.sensor_vector_end   = [5 10];
inputs.flag_search_type = 0;

% SHORT format checking
clear expected
expected.distance = 5;
expected.location = [5 10];
expected.path_segment = 1;
expected.t = 0.5;
expected.u = 1;

actual = struct;
[actual.distance, actual.location, actual.path_segment, actual.t, actual.u] = ...
    fcn_Path_findProjectionHitOntoPath(...
    inputs.path,inputs.sensor_vector_start,inputs.sensor_vector_end,...
    inputs.flag_search_type,inputs.fig_num);

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)
%fcn_INTERNAL_printResults(actual.distance,actual.location);


%% Single point intersection (2), flag (0), test 6 -  intersection at start of path

fig_num = 20006;
fprintf(1,'Figure: %.0f : Single point intersection (2), flag (0) - intersection at start of path  \n',fig_num);
figure(fig_num); clf;
plotting.FigureExpected = 1;

clear inputs
inputs.fig_num = fig_num;
inputs.path = [0 10; 10 10];
inputs.sensor_vector_start = [0 0];
inputs.sensor_vector_end   = [0 15];
inputs.flag_search_type = 0;

% SHORT format checking
clear expected
expected.distance = 10;
expected.location = [0 10];
expected.path_segment = 1;
expected.t = 0;
expected.u = 0.6667;

actual = struct;
[actual.distance, actual.location, actual.path_segment, actual.t, actual.u] = ...
    fcn_Path_findProjectionHitOntoPath(...
    inputs.path,inputs.sensor_vector_start,inputs.sensor_vector_end,...
    inputs.flag_search_type,inputs.fig_num);

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)
%fcn_INTERNAL_printResults(actual.distance,actual.location);

%% Single point intersection (2), flag (0), test 7 -  intersection at end of path

fig_num = 20007;
fprintf(1,'Figure: %.0f : Single point intersection (2), flag (0), - intersection at end of path  \n',fig_num);
figure(fig_num); clf;
plotting.FigureExpected = 1;

clear inputs
inputs.fig_num = fig_num;
inputs.path = [0 10; 10 10];
inputs.sensor_vector_start = [10 0];
inputs.sensor_vector_end   = [10 15];
inputs.flag_search_type = 0;

% SHORT format checking
clear expected
expected.distance = 10;
expected.location = [10 10];
expected.path_segment = 1;
expected.t = 1;
expected.u = 0.6667;

actual = struct;
[actual.distance, actual.location, actual.path_segment, actual.t, actual.u] = ...
    fcn_Path_findProjectionHitOntoPath(...
    inputs.path,inputs.sensor_vector_start,inputs.sensor_vector_end,...
    inputs.flag_search_type,inputs.fig_num);

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)
%fcn_INTERNAL_printResults(actual.distance,actual.location);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____ _             _                    _       _     _       _                          _   _
%  / ____(_)           | |                  (_)     | |   (_)     | |                        | | (_)
% | (___  _ _ __   __ _| | ___   _ __   ___  _ _ __ | |_   _ _ __ | |_ ___ _ __ ___  ___  ___| |_ _  ___  _ __
%  \___ \| | '_ \ / _` | |/ _ \ | '_ \ / _ \| | '_ \| __| | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __| |/ _ \| '_ \
%  ____) | | | | | (_| | |  __/ | |_) | (_) | | | | | |_  | | | | | ||  __/ |  \__ \  __/ (__| |_| | (_) | | | |
% |_____/|_|_| |_|\__, |_|\___| | .__/ \___/|_|_| |_|\__| |_|_| |_|\__\___|_|  |___/\___|\___|\__|_|\___/|_| |_|
%                  __/ |        | |
%                 |___/         |_|
%   __ _                        __
%  / _| |              ______  /_ |
% | |_| | __ _  __ _  |______|  | |
% |  _| |/ _` |/ _` |  ______   | |
% | | | | (_| | (_| | |______|  | |
% |_| |_|\__,_|\__, |           |_|
%               __/ |
%              |___/
% http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Single%20point%20intersection%20%0Aflag%20%3D%201
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;


%% Single point intersection (2), flag (1) -  using flag to find intersection beyond sensor range, one path

fig_num = 21001;
fprintf(1,'Figure: %.0f : Single point intersection (2), flag (1), - using flag to find intersection beyond sensor range, one path  \n',fig_num);
figure(fig_num); clf;
plotting.FigureExpected = 1;

clear inputs
inputs.fig_num = fig_num;
inputs.path = [0 10; 10 10];
inputs.sensor_vector_start = [5 0];
inputs.sensor_vector_end   = [5 5];
inputs.flag_search_type = 1;

% SHORT format checking
clear expected
expected.distance = 10;
expected.location = [5 10];
expected.path_segment = 1;
expected.t = 0.5;
expected.u = 2;

actual = struct;
[actual.distance, actual.location, actual.path_segment, actual.t, actual.u] = ...
    fcn_Path_findProjectionHitOntoPath(...
    inputs.path,inputs.sensor_vector_start,inputs.sensor_vector_end,...
    inputs.flag_search_type,inputs.fig_num);

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)
%fcn_INTERNAL_printResults(actual.distance,actual.location);


%% Single point intersection (2), flag (1) -  using flag to find intersection before sensor range, one path

fig_num = 21002;
fprintf(1,'Figure: %.0f : Single point intersection (2), flag (1), - using flag to find intersection before sensor range, one path  \n',fig_num);
figure(fig_num); clf;
plotting.FigureExpected = 1;

clear inputs
inputs.fig_num = fig_num;
inputs.path = [0 10; 10 10];
inputs.sensor_vector_start = [5 15];
inputs.sensor_vector_end   = [5 20];
inputs.flag_search_type = 1;

% SHORT format checking
clear expected
expected.distance = -5;
expected.location = [5 10];
expected.path_segment = 1;
expected.t = 0.5;
expected.u = -1;

actual = struct;
[actual.distance, actual.location, actual.path_segment, actual.t, actual.u] = ...
    fcn_Path_findProjectionHitOntoPath(...
    inputs.path,inputs.sensor_vector_start,inputs.sensor_vector_end,...
    inputs.flag_search_type,inputs.fig_num);

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)
%fcn_INTERNAL_printResults(actual.distance,actual.location);


%% Single point intersection (2), flag (1) -  using flag to find intersection beyond sensor range, multi path

fig_num = 21003;
fprintf(1,'Figure: %.0f : Single point intersection (2), flag (1), - using flag to find intersection beyond sensor range, multi path  \n',fig_num);
figure(fig_num); clf;
plotting.FigureExpected = 1;

clear inputs
inputs.fig_num = fig_num;
inputs.path = [0 10; 10 10; 0 15];
inputs.sensor_vector_start = [5 0];
inputs.sensor_vector_end   = [5 5];
inputs.flag_search_type = 1;

% SHORT format checking
clear expected
expected.distance = 10;
expected.location = [5 10];
expected.path_segment = 1;
expected.t = 0.5;
expected.u = 2;

actual = struct;
[actual.distance, actual.location, actual.path_segment, actual.t, actual.u] = ...
    fcn_Path_findProjectionHitOntoPath(...
    inputs.path,inputs.sensor_vector_start,inputs.sensor_vector_end,...
    inputs.flag_search_type,inputs.fig_num);

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)
%fcn_INTERNAL_printResults(actual.distance,actual.location);


%% Single point intersection (2), flag (1) -  using flag to find intersection before sensor range, multi path

fig_num = 21004;
fprintf(1,'Figure: %.0f : Single point intersection (2), flag (1), - using flag to find intersection before sensor range, multi path  \n',fig_num);
figure(fig_num); clf;
plotting.FigureExpected = 1;

clear inputs
inputs.fig_num = fig_num;
inputs.path = [0 10; 10 10; 0 5];
inputs.sensor_vector_start = [5 15];
inputs.sensor_vector_end   = [5 20];
inputs.flag_search_type = 1;

% SHORT format checking
clear expected
expected.distance = -5;
expected.location = [5 10];
expected.path_segment = 1;
expected.t = 0.5;
expected.u = -1;

actual = struct;
[actual.distance, actual.location, actual.path_segment, actual.t, actual.u] = ...
    fcn_Path_findProjectionHitOntoPath(...
    inputs.path,inputs.sensor_vector_start,inputs.sensor_vector_end,...
    inputs.flag_search_type,inputs.fig_num);

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)
%fcn_INTERNAL_printResults(actual.distance,actual.location);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____ _             _                    _       _     _       _                          _   _
%  / ____(_)           | |                  (_)     | |   (_)     | |                        | | (_)
% | (___  _ _ __   __ _| | ___   _ __   ___  _ _ __ | |_   _ _ __ | |_ ___ _ __ ___  ___  ___| |_ _  ___  _ __
%  \___ \| | '_ \ / _` | |/ _ \ | '_ \ / _ \| | '_ \| __| | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __| |/ _ \| '_ \
%  ____) | | | | | (_| | |  __/ | |_) | (_) | | | | | |_  | | | | | ||  __/ |  \__ \  __/ (__| |_| | (_) | | | |
% |_____/|_|_| |_|\__, |_|\___| | .__/ \___/|_|_| |_|\__| |_|_| |_|\__\___|_|  |___/\___|\___|\__|_|\___/|_| |_|
%                  __/ |        | |
%                 |___/         |_|
%   __ _                        ___
%  / _| |              ______  |__ \
% | |_| | __ _  __ _  |______|    ) |
% |  _| |/ _` |/ _` |  ______    / /
% | | | | (_| | (_| | |______|  / /_
% |_| |_|\__,_|\__, |          |____|
%               __/ |
%              |___/
% http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Single%20point%20intersection%20%0Aflag%20%3D%202
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;

% 2: same as flag_search_type 0, except all intersections if there is overlap, and only first and last contact points if infinte


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____ _             _                    _       _     _       _                          _   _
%  / ____(_)           | |                  (_)     | |   (_)     | |                        | | (_)
% | (___  _ _ __   __ _| | ___   _ __   ___  _ _ __ | |_   _ _ __ | |_ ___ _ __ ___  ___  ___| |_ _  ___  _ __
%  \___ \| | '_ \ / _` | |/ _ \ | '_ \ / _ \| | '_ \| __| | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __| |/ _ \| '_ \
%  ____) | | | | | (_| | |  __/ | |_) | (_) | | | | | |_  | | | | | ||  __/ |  \__ \  __/ (__| |_| | (_) | | | |
% |_____/|_|_| |_|\__, |_|\___| | .__/ \___/|_|_| |_|\__| |_|_| |_|\__\___|_|  |___/\___|\___|\__|_|\___/|_| |_|
%                  __/ |        | |
%                 |___/         |_|
%   __ _                        ____
%  / _| |              ______  |___ \
% | |_| | __ _  __ _  |______|   __) |
% |  _| |/ _` |/ _` |  ______   |__ <
% | | | | (_| | (_| | |______|  ___) |
% |_| |_|\__,_|\__, |          |____/
%               __/ |
%              |___/
% http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Single%20point%20intersection%20%0Aflag%20%3D%203
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;

%% Single point intersection (2), flag (3) - intersection by extending path

fig_num = 23004;
fprintf(1,'Figure: %.0f : Single point intersection (2), flag (3), - intersection by extending single path  \n',fig_num);
figure(fig_num); clf;
plotting.FigureExpected = 1;

clear inputs
inputs.fig_num = fig_num;
inputs.path = [0 10; 10 10];
inputs.sensor_vector_start = [15 0];
inputs.sensor_vector_end   = [15 15];
inputs.flag_search_type = 3;

% SHORT format checking
clear expected
expected.distance = 10;
expected.location = [15 10];
expected.path_segment = 1;
expected.t = 1.5;
expected.u = 0.6667;

actual = struct;
[actual.distance, actual.location, actual.path_segment, actual.t, actual.u] = ...
    fcn_Path_findProjectionHitOntoPath(...
    inputs.path,inputs.sensor_vector_start,inputs.sensor_vector_end,...
    inputs.flag_search_type,inputs.fig_num);

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)
%fcn_INTERNAL_printResults(actual.distance,actual.location);



%% Single point intersection (2), flag (3) - intersection by extending 2 paths

fig_num = 23004;
fprintf(1,'Figure: %.0f : Single point intersection (2), flag (3), - intersection by extending 2 paths  \n',fig_num);
figure(fig_num); clf;
plotting.FigureExpected = 1;

clear inputs
inputs.fig_num = fig_num;
inputs.path = [0 10; 10 10; 10 15; 5 20];
inputs.sensor_vector_start = [15 0];
inputs.sensor_vector_end   = [15 15];
inputs.flag_search_type = 3;

% SHORT format checking
clear expected
expected.distance = 10;
expected.location = [15 10];
expected.path_segment = 1;
expected.t = 1.5;
expected.u = 0.6667;

actual = struct;
[actual.distance, actual.location, actual.path_segment, actual.t, actual.u] = ...
    fcn_Path_findProjectionHitOntoPath(...
    inputs.path,inputs.sensor_vector_start,inputs.sensor_vector_end,...
    inputs.flag_search_type,inputs.fig_num);

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)
%fcn_INTERNAL_printResults(actual.distance,actual.location);




%% Advanced test 1 - intersection beyond a sensor's range with flag (for flag_search_type = 4)
fprintf(1,'Intersection beyond sensor range  \n');
path = [0 5; 4 5; 8 2];

sensor_vector_start = [4 0];
sensor_vector_end   = [4 2];
fig_debugging = 34345;

flag_search_type =0;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
title('Intersection beyond sensor range result, Flag = 0')
fcn_INTERNAL_printResults(distance,location);

assert(isnan(distance));
assert(all(isnan(location)));
assert(isnan(path_segments));

fig_debugging = 34346;

flag_search_type = 4;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
title('Intersection beyond sensor range result, Flag = 4')
fcn_INTERNAL_printResults(distance,location);

assert(isequal(round(distance,4),5));
assert(isequal(round(location,4),[4 5]));
assert(isequal(path_segments,1));

% Test the negative condition
sensor_vector_start = [4 6];
sensor_vector_end   = [4 8];
fig_debugging = 34347;

flag_search_type = 0;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
title('Intersection beyond sensor range result, Flag = 0')
fcn_INTERNAL_printResults(distance,location);

assert(isnan(distance));
assert(all(isnan(location)));
assert(isnan(path_segments));

fig_debugging = 34348;
flag_search_type = 4;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
title('Intersection beyond sensor range result, Flag = 4')
fcn_INTERNAL_printResults(distance,location);

assert(isequal(round(distance,4),-1));
assert(isequal(round(location,4),[4 5]));
assert(isequal(path_segments,1));


%% Simple test - Extend both the vectors (for flag_search_type = 4)
fig_debugging = 34343;
fprintf(1,'Simple test result, flag = 4, long sensor: \n');
path = [-4 10; 2 10];
sensor_vector_start = [0 0];
sensor_vector_end   = 2*[4 6];

flag_search_type =4;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
fcn_INTERNAL_printResults(distance,location);

assert(isequal(round(distance,4),12.0185));
assert(isequal(round(location,4),[6.6667   10.0000]));
assert(isequal(path_segments,1));

%% Simple test - Extend both the vectors (for flag_search_type = 4)
fig_debugging = 34344;
fprintf(1,'Simple test result, flag = 4, short sensor: \n');
path = [-4 10; 2 10];
sensor_vector_start = [0 0];
sensor_vector_end   = [4 6];

flag_search_type =4;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
fcn_INTERNAL_printResults(distance,location);

assert(isequal(round(distance,4),12.0185));
assert(isequal(round(location,4),[6.6667   10.0000]));
assert(isequal(path_segments,1));

%% Multi-hit test - flag_search_type = 3
fig_debugging = 343433;
fprintf(1,'No intersection  \n');
path = [-4 0; -2 0; 0 -2; 2 0; 4 0];
sensor_vector_start = [0 0];
sensor_vector_end   = [0 4];

flag_search_type = 3;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
fcn_INTERNAL_printResults(distance,location);

assert(isequal(round(distance,4),0));
assert(isequal(round(location,4),[0   0]));
assert(isequal(path_segments,1));

%% Multi-hit test - flag_search_type = 4
fig_debugging = 343434;
fprintf(1,'No intersection  \n');
path = [-4 0; -2 0; 0 -2; 2 0; 4 0];
sensor_vector_start = [0 0];
sensor_vector_end   = [0 4];

flag_search_type = 4;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
fcn_INTERNAL_printResults(distance,location);

assert(isequal(round(distance,4),0));
assert(isequal(round(location,4),[0   0]));
assert(isequal(path_segments,1));

%% Multi-hit test - flag_search_type = 4
fig_debugging = 343438;
fprintf(1,'No intersection  \n');
path = [-4 0; -2 0; 0 -2; 2 0; 4 0];
sensor_vector_start = [0 1];
sensor_vector_end   = [0 4];

flag_search_type = 4;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
fcn_INTERNAL_printResults(distance,location);

assert(isequal(round(distance,4),-1));
assert(isequal(round(location,4),[0   0]));
assert(isequal(path_segments,1));


% %% Non intersection cases
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  _   _               _____       _                          _   _
% % | \ | |             |_   _|     | |                        | | (_)
% % |  \| | ___  _ __     | |  _ __ | |_ ___ _ __ ___  ___  ___| |_ _  ___  _ __  ___
% % | . ` |/ _ \| '_ \    | | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __| |/ _ \| '_ \/ __|
% % | |\  | (_) | | | |  _| |_| | | | ||  __/ |  \__ \  __/ (__| |_| | (_) | | | \__ \
% % |_| \_|\___/|_| |_| |_____|_| |_|\__\___|_|  |___/\___|\___|\__|_|\___/|_| |_|___/
% %
% % See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Non%20Intersections
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % All non intersection figures start with the number 3
%
% close all;
%
% %% Non intersection (3), flag (0), test 1
% fig_num = 30001;
% figure(fig_num); clf;
%
% fprintf(1,'\nNon intersection (3), flag (0), test 1  \n');
%
% path = [-4 10; 2 10];
% sensor_vector_start = [0 0];
% sensor_vector_end   = [5 12];
% flag_search_type = 0;
%
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_num);
% fcn_INTERNAL_printResults(distance,location);
%
%
% % Check variable types
% assert(isnumeric(distance));
% assert(isnumeric(location));
%
% % Check variable sizes
% assert(isequal(size(distance),[1 1]));
% assert(isequal(size(location),[1 2]));
%
% % Check variable values
% assert(isnan(distance));
% assert(all(isnan(location)));
%
% % Make sure plot opened up
% assert(isequal(get(gcf,'Number'),fig_num));
%
%
% %% Simple test - No intersection (for flag_search_type = 3)
% fig_debugging = 34342;
% fprintf(1,'Simple test result, flag = 3, short sensor: \n');
% path = [-4 10; 2 10];
% sensor_vector_start = [0 0];
% sensor_vector_end   = [4 6];
%
% flag_search_type =3;
% [distance,location,path_segments] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isnan(distance));
% assert(all(isnan(location)));
% assert(isnan(path_segments));
%
% %% Parallel test for flag_search_type = 3
% fig_debugging = 3434666;
% fprintf(1,'Parallel test for flag=3: \n');
% path = [0 0; 10 0];
% sensor_vector_start = [0 1];
% sensor_vector_end   = [10 1];
%
% flag_search_type =3;
% [distance,location,path_segments] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isnan(distance));
% assert(all(isnan(location)));
% assert(isnan(path_segments));
%
% %% Parallel test for flag_search_type = 4
% fig_debugging = 3434777;
% fprintf(1,'Parallel test for flag=4: \n');
% path = [0 0; 10 0];
% sensor_vector_start = [0 1];
% sensor_vector_end   = [10 1];
%
% flag_search_type =4;
% [distance,location,path_segments] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isnan(distance));
% assert(all(isnan(location)));
% assert(isnan(path_segments));
%
%
% %% Multi-hit test - flag_search_type = 3
% fig_debugging = 343437;
% fprintf(1,'No intersection  \n');
% path = [-4 0; -2 0; 0 -2; 2 0; 4 0];
% sensor_vector_start = [0 1];
% sensor_vector_end   = [0 4];
%
% flag_search_type = 3;
% [distance,location,path_segments] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isnan(distance));
% assert(all(isnan(location)));
% assert(isnan(path_segments));
%
%
%
% %% Infinite intersection cases
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  _____        __ _       _ _         _____       _                          _   _
% % |_   _|      / _(_)     (_) |       |_   _|     | |                        | | (_)
% %   | |  _ __ | |_ _ _ __  _| |_ ___    | |  _ __ | |_ ___ _ __ ___  ___  ___| |_ _  ___  _ __  ___
% %   | | | '_ \|  _| | '_ \| | __/ _ \   | | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __| |/ _ \| '_ \/ __|
% %  _| |_| | | | | | | | | | | ||  __/  _| |_| | | | ||  __/ |  \__ \  __/ (__| |_| | (_) | | | \__ \
% % |_____|_| |_|_| |_|_| |_|_|\__\___| |_____|_| |_|\__\___|_|  |___/\___|\___|\__|_|\___/|_| |_|___/
% %
% % See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Infinite%20Intersections
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % All infinte intersection figures start with the number 4
%
% close all;
%
% %% Simple test 3 - multiple intersections, returns only the first one
% fprintf(1,'Multiple intersections  \n');
% path = [0 10; 10 10; 0 6; 10 6; 0 2];
% sensor_vector_start = [0 0];
% sensor_vector_end   = [5 12];
% fig_debugging = 2343;
% flag_search_type = 0;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),2.6));
% assert(isequal(round(location,4),[1.0000    2.4000]));
%
%
%
% %% Simple test 7 - identically overlapping colinear
% fprintf(1,'identically overlapping colinear  \n');
% path = [0 10; 10 10];
% sensor_vector_start = [0 10];
% sensor_vector_end   = [10 10];
% fig_debugging = 2343;
% flag_search_type = 0;
% [distance,location,path_segments] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
% fcn_INTERNAL_printMoreResults(distance,location,path_segments);
%
% assert(isequal(round(distance,4),0));
% assert(isequal(round(location,4),[0 10]));
%
%
% %% Simple test 8 - partially overlapping colinear 1
% fprintf(1,'Partially overlapping colinear  \n');
% path = [0 10; 10 10];
% sensor_vector_start = [-2 10];
% sensor_vector_end   = [10 10];
% fig_debugging = 2343;
% flag_search_type = 0;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
%
% assert(isequal(round(distance,4),2));
% assert(isequal(round(location,4),[0 10]));
%
% %% Simple test 9 - partially overlapping colinear 1
% fprintf(1,'Partially overlapping colinear  \n');
% path = [0 10; 10 10];
% sensor_vector_start = [-2 10];
% sensor_vector_end   = [5 10];
% fig_debugging = 2343;
% flag_search_type = 0;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),2));
% assert(isequal(round(location,4),[0 10]));
%
% %% Simple test 10 - partially overlapping colinear 1
% fprintf(1,'Partially overlapping colinear  \n');
% path = [0 10; 10 10];
% sensor_vector_start = [3 10];
% sensor_vector_end   = [5 10];
% fig_debugging = 2343;
% flag_search_type = 0;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
%
% assert(isequal(round(distance,4),0));
% assert(isequal(round(location,4),[3 10]));
%
% %% Simple test 11 - partially overlapping colinear 1
% fprintf(1,'Partially overlapping colinear  \n');
% path = [0 10; 10 10];
% sensor_vector_start = [3 10];
% sensor_vector_end   = [15 10];
% fig_debugging = 2343;
% flag_search_type = 0;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),0));
% assert(isequal(round(location,4),[3 10]));
%
% %% Simple test 12 - super overlapping colinear 1
% fprintf(1,'Super overlapping colinear  \n');
% path = [0 10; 10 10];
% sensor_vector_start = [-3 10];
% sensor_vector_end   = [15 10];
% fig_debugging = 2343;
% flag_search_type = 0;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
%
% assert(isequal(round(distance,4),3));
% assert(isequal(round(location,4),[0 10]));
%
% %% Simple test 13 - end overlapping colinear 1
% fprintf(1,'End overlapping colinear  \n');
% path = [0 10; 10 10];
% sensor_vector_start = [-3 10];
% sensor_vector_end   = [0 10];
% fig_debugging = 2343;
% flag_search_type = 0;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),3));
% assert(isequal(round(location,4),[0 10]));
%
% %% Simple test 14 - end overlapping colinear 2
% fprintf(1,'End overlapping colinear  \n');
% path = [0 10; 10 10];
% sensor_vector_start = [10 10];
% sensor_vector_end   = [15 10];
% fig_debugging = 2343;
% flag_search_type = 0;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),0));
% assert(isequal(round(location,4),[10 10]));
%
%
% %% Simple test 27 - identically overlapping colinear BACKWARDS
% fprintf(1,'identically overlapping colinear BACKWARDS  \n');
% path = [0 10; 10 10];
% sensor_vector_end = [0 10];
% sensor_vector_start   = [10 10];
% fig_debugging = 2343;
% flag_search_type = 0;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),0));
% assert(isequal(round(location,4),[10 10]));
%
% %% Simple test 28 - partially overlapping colinear 1 BACKWARDS
% fprintf(1,'Partially overlapping colinear BACKWARDS  \n');
% path = [0 10; 10 10];
% sensor_vector_end = [-2 10];
% sensor_vector_start   = [10 10];
% fig_debugging = 2343;
% flag_search_type = 0;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),0));
% assert(isequal(round(location,4),[10 10]));
%
% %% Simple test 29 - partially overlapping colinear 1 BACKWARDS
% fprintf(1,'Partially overlapping colinear BACKWARDS  \n');
% path = [0 10; 10 10];
% sensor_vector_end = [-2 10];
% sensor_vector_start   = [5 10];
% fig_debugging = 2343;
% flag_search_type = 0;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),0));
% assert(isequal(round(location,4),[5 10]));
%
% %% Simple test 30 - partially overlapping colinear 1 BACKWARDS
% fprintf(1,'Partially overlapping colinear BACKWARDS  \n');
% path = [0 10; 10 10];
% sensor_vector_end = [3 10];
% sensor_vector_start   = [5 10];
% fig_debugging = 2343;
% flag_search_type = 0;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),0));
% assert(isequal(round(location,4),[5 10]));
%
% %% Simple test 31 - partially overlapping colinear 1 BACKWARDS
% fprintf(1,'Partially overlapping colinear BACKWARDS  \n');
% path = [0 10; 10 10];
% sensor_vector_end = [3 10];
% sensor_vector_start   = [15 10];
% fig_debugging = 2343;
% flag_search_type = 0;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),5));
% assert(isequal(round(location,4),[10 10]));
%
% %% Simple test 32 - super overlapping colinear 1 BACKWARDS
% fprintf(1,'Super overlapping colinear BACKWARDS  \n');
% path = [0 10; 10 10];
% sensor_vector_end = [-3 10];
% sensor_vector_start   = [15 10];
% fig_debugging = 2343;
% flag_search_type = 0;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),5));
% assert(isequal(round(location,4),[10 10]));
%
%
% %% Simple test 33 - end overlapping colinear 1 BACKWARDS
% fprintf(1,'End overlapping colinear  \n');
% path = [0 10; 10 10];
% sensor_vector_end = [-3 10];
% sensor_vector_start   = [0 10];
% fig_debugging = 2343;
% flag_search_type = 0;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),0));
% assert(isequal(round(location,4),[0 10]));
%
% %% Simple test 34 - end overlapping colinear 2 BACKWARDS
% fprintf(1,'End overlapping colinear BACKWARDS  \n');
% path = [0 10; 10 10];
% sensor_vector_end = [10 10];
% sensor_vector_start   = [15 10];
% fig_debugging = 2343;
% flag_search_type = 0;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
%
% assert(isequal(round(distance,4),5));
% assert(isequal(round(location,4),[10 10]));
%
%
% %% Simple test 15 - non overlapping colinear 1
% fprintf(1,'Non overlapping colinear  \n');
% path = [0 10; 10 10];
% sensor_vector_start = [-3 10];
% sensor_vector_end   = [-1 10];
% fig_debugging = 2343;
% flag_search_type = 0;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isnan(distance));
% assert(all(isnan(location)));
%
%
% %% Simple test 15 - non overlapping colinear 2
% fprintf(1,'Non overlapping colinear  \n');
% path = [0 10; 10 10];
% sensor_vector_start = [13 10];
% sensor_vector_end   = [15 10];
% fig_debugging = 2343;
% flag_search_type = 0;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isnan(distance));
% assert(all(isnan(location)));
%
%
%
%
% %% Multi-hit tests
% %   __  __       _ _   _ _    _ _ _
% %  |  \/  |     | | | (_) |  | (_) |
% %  | \  / |_   _| | |_ _| |__| |_| |_
% %  | |\/| | | | | | __| |  __  | | __|
% %  | |  | | |_| | | |_| | |  | | | |_
% %  |_|  |_|\__,_|_|\__|_|_|  |_|_|\__|
% %
% %
%
% close all;
%
% %% Advanced test 2 - multiple intersections
% fprintf(1,'Single intersections reporting only first  \n');
% path = [0 10; 10 10; 0 6; 10 6; 0 2];
% sensor_vector_start = [0 0];
% sensor_vector_end   = [5 12];
% fig_debugging = 23487;
% flag_search_type = 0;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),2.6));
% assert(isequal(round(location,4),[1 2.4]));
%
% fprintf(1,'Multiple intersections reporting all results: \n');
% path = [0 10; 10 10; 0 6; 10 6; 0 2];
% sensor_vector_start = [0 0];
% sensor_vector_end   = [5 12];
% fig_debugging = 23488;
% flag_search_type = 2;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance',4),[10.8333    7.8000    6.5000    2.6000]));
% assert(isequal(round(location,4),[    4.1667   10.0000
%     3.0000    7.2000
%     2.5000    6.0000
%     1.0000    2.4000]));
%
% %% Advanced test 3 - multiple intersections possible, but no hits
% fprintf(1,'Multiple intersections possible but no hits, reporting all results: \n');
% path = [0 10; 10 10; 0 6; 10 6; 0 2];
% sensor_vector_start = [0 0];
% sensor_vector_end   = [0.5 1.2];
% fig_debugging = 23499;
% flag_search_type = 2;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isempty(distance));
% assert(isempty(location));
%
% %% Advanced test 4 - multiple intersections possible, but few hits
% fprintf(1,'Multiple intersections possible but few hits, reporting all results: \n');
% path = [0 10; 10 10; 0 6; 10 6; 0 2];
% sensor_vector_start = [0 0];
% sensor_vector_end   = [2.5 6];
% fig_debugging = 1010;
% flag_search_type = 2;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),[6.5; 2.6]));
% assert(isequal(round(location,4),[2.5000    6.0000
%     1.0000    2.4000]));
%
%
% %% Multi-hit overlapping tests
% %   __  __       _ _   _ _    _ _ _    ____                 _                   _
% %  |  \/  |     | | | (_) |  | (_) |  / __ \               | |                 (_)
% %  | \  / |_   _| | |_ _| |__| |_| |_| |  | |_   _____ _ __| | __ _ _ __  _ __  _ _ __   __ _
% %  | |\/| | | | | | __| |  __  | | __| |  | \ \ / / _ \ '__| |/ _` | '_ \| '_ \| | '_ \ / _` |
% %  | |  | | |_| | | |_| | |  | | | |_| |__| |\ V /  __/ |  | | (_| | |_) | |_) | | | | | (_| |
% %  |_|  |_|\__,_|_|\__|_|_|  |_|_|\__|\____/  \_/ \___|_|  |_|\__,_| .__/| .__/|_|_| |_|\__, |
% %                                                                  | |   | |             __/ |
% %                                                                  |_|   |_|            |___/
%
% close all;
%
% %% Advanced Multihit Overlapping test - identically overlapping colinear
% fprintf(1,'identically overlapping colinear  \n');
% path = [0 10; 10 10];
% sensor_vector_start = [0 10];
% sensor_vector_end   = [10 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segments] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
% fcn_INTERNAL_printMoreResults(distance,location,path_segments);
%
% assert(isequal(round(distance,4),[0; 10]));
% assert(isequal(round(location,4),[     0    10
%     10    10]));
%
%
% %% Advanced Multihit Overlapping  test 8 - partially overlapping colinear 1
% fprintf(1,'Partially overlapping colinear  \n');
% path = [0 10; 10 10];
% sensor_vector_start = [-2 10];
% sensor_vector_end   = [10 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),[2; 12]));
% assert(isequal(round(location,4),[     0    10
%     10    10]));
%
% %% Advanced Multihit Overlapping  test 9 - partially overlapping colinear 1
% fprintf(1,'Partially overlapping colinear  \n');
% path = [0 10; 10 10];
% sensor_vector_start = [-2 10];
% sensor_vector_end   = [5 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),[2; 7]));
% assert(isequal(round(location,4),[0   10.0000
%     5.0000   10.0000]));
%
% %% Advanced Multihit Overlapping  test 10 - partially overlapping colinear 1
% fprintf(1,'Partially overlapping colinear  \n');
% path = [0 10; 10 10];
% sensor_vector_start = [3 10];
% sensor_vector_end   = [5 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),[0; 2]));
% assert(isequal(round(location,4),[     3    10
%      5    10]));
%
% %% Advanced Multihit Overlapping  test 11 - partially overlapping colinear 1
% fprintf(1,'Partially overlapping colinear  \n');
% path = [0 10; 10 10];
% sensor_vector_start = [3 10];
% sensor_vector_end   = [15 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),[0; 7]));
% assert(isequal(round(location,4),[     3    10
%     10    10]));
%
% %% Advanced Multihit Overlapping  test 12 - super overlapping colinear 1
% fprintf(1,'Super overlapping colinear  \n');
% path = [0 10; 10 10];
% sensor_vector_start = [-3 10];
% sensor_vector_end   = [15 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),[3; 13]));
% assert(isequal(round(location,4),[0    10
%     10    10]));
%
% %% Advanced Multihit Overlapping  test 13 - end overlapping colinear 1
% fprintf(1,'End overlapping colinear  \n');
% path = [0 10; 10 10];
% sensor_vector_start = [-3 10];
% sensor_vector_end   = [0 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),3));
% assert(isequal(round(location,4),[ 0    10]));
%
% %% Advanced Multihit Overlapping  test 14 - end overlapping colinear 2
% fprintf(1,'End overlapping colinear  \n');
% path = [0 10; 10 10];
% sensor_vector_start = [10 10];
% sensor_vector_end   = [15 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),0));
% assert(isequal(round(location,4),[10 10]));
%
%
% %% Advanced Multihit Overlapping  test 27 - identically overlapping colinear BACKWARDS
% fprintf(1,'identically overlapping colinear BACKWARDS  \n');
% path = [0 10; 10 10];
% sensor_vector_end = [0 10];
% sensor_vector_start   = [10 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),[0; 10]));
% assert(isequal(round(location,4),[    10    10
%      0    10]));
%
% %% Advanced Multihit Overlapping  test 28 - partially overlapping colinear 1 BACKWARDS
% fprintf(1,'Partially overlapping colinear BACKWARDS  \n');
% path = [0 10; 10 10];
% sensor_vector_end = [-2 10];
% sensor_vector_start   = [10 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),[0; 10]));
% assert(isequal(round(location,4),[    10    10
%      0    10]));
%
% %% Advanced Multihit Overlapping  test 29 - partially overlapping colinear 1 BACKWARDS
% fprintf(1,'Partially overlapping colinear BACKWARDS  \n');
% path = [0 10; 10 10];
% sensor_vector_end = [-2 10];
% sensor_vector_start   = [5 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),[0; 5]));
% assert(isequal(round(location,4),[    5    10
%      0    10]));
%
% %% Advanced Multihit Overlapping  test 30 - partially overlapping colinear 1 BACKWARDS
% fprintf(1,'Partially overlapping colinear BACKWARDS  \n');
% path = [0 10; 10 10];
% sensor_vector_end = [3 10];
% sensor_vector_start   = [5 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),[0; 2]));
% assert(isequal(round(location,4),[    5    10
%      3    10]));
%
% %% Advanced Multihit Overlapping  test 31 - partially overlapping colinear 1 BACKWARDS
% fprintf(1,'Partially overlapping colinear BACKWARDS  \n');
% path = [0 10; 10 10];
% sensor_vector_end = [3 10];
% sensor_vector_start   = [15 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),[5; 12]));
% assert(isequal(round(location,4),[    10    10
%      3    10]));
%
% %% Advanced Multihit Overlapping  test 32 - super overlapping colinear 1 BACKWARDS
% fprintf(1,'Super overlapping colinear BACKWARDS  \n');
% path = [0 10; 10 10];
% sensor_vector_end = [-3 10];
% sensor_vector_start   = [15 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),[5; 15]));
% assert(isequal(round(location,4),[    10    10
%      0    10]));
%
% %% Advanced Multihit Overlapping  test 33 - end overlapping colinear 1 BACKWARDS
% fprintf(1,'End overlapping colinear  \n');
% path = [0 10; 10 10];
% sensor_vector_end = [-3 10];
% sensor_vector_start   = [0 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),0));
% assert(isequal(round(location,4),[    0    10]));
%
% %% Advanced Multihit Overlapping  test 34 - end overlapping colinear 2 BACKWARDS
% fprintf(1,'End overlapping colinear BACKWARDS  \n');
% path = [0 10; 10 10];
% sensor_vector_end = [10 10];
% sensor_vector_start   = [15 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isequal(round(distance,4),5));
% assert(isequal(round(location,4),[10    10]));
%
%
% %% Advanced Multihit Overlapping  test 15 - non overlapping colinear 1
% fprintf(1,'Non overlapping colinear  \n');
% path = [0 10; 10 10];
% sensor_vector_start = [-3 10];
% sensor_vector_end   = [-1 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isempty(distance));
% assert(isempty(location));
%
% %% Advanced Multihit Overlapping  test 15 - non overlapping colinear 2
% fprintf(1,'Non overlapping colinear  \n');
% path = [0 10; 10 10];
% sensor_vector_start = [13 10];
% sensor_vector_end   = [15 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printResults(distance,location);
%
% assert(isempty(distance));
% assert(isempty(location));
%
%
% %% multiple hits with multiple paths
% %   __  __       _ _   _ _    _ _ _   __  __       _ _ _   _____      _   _
% %  |  \/  |     | | | (_) |  | (_) | |  \/  |     | (_) | |  __ \    | | | |
% %  | \  / |_   _| | |_ _| |__| |_| |_| \  / |_   _| |_| |_| |__) |_ _| |_| |__
% %  | |\/| | | | | | __| |  __  | | __| |\/| | | | | | | __|  ___/ _` | __| '_ \
% %  | |  | | |_| | | |_| | |  | | | |_| |  | | |_| | | | |_| |  | (_| | |_| | | |
% %  |_|  |_|\__,_|_|\__|_|_|  |_|_|\__|_|  |_|\__,_|_|_|\__|_|   \__,_|\__|_| |_|
% %
% %
%
% close all;
%
% %% Advanced Multihit Overlapping test - identically overlapping colinear
% fprintf(1,'identically overlapping colinear  \n');
% path = [0 10; 10 10; 12 8; 14 10; 15 10];
% sensor_vector_start = [0 10];
% sensor_vector_end   = [10 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segments] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printMoreResults(distance,location,path_segments);
%
% assert(isequal(round(distance,4),[     0
%     10
%     10]));
% assert(isequal(round(location,4),[     0    10
%     10    10
%     10    10]));
%
% %% Advanced Multihit Overlapping  test 8 - partially overlapping colinear 1
% fprintf(1,'Partially overlapping colinear  \n');
% path = [0 10; 10 10; 12 8; 14 10; 15 10];
% sensor_vector_start = [-2 10];
% sensor_vector_end   = [10 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segments] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printMoreResults(distance,location,path_segments);
%
% assert(isequal(round(distance,4),[     2
%     12
%     12]));
% assert(isequal(round(location,4),[     0    10
%     10    10
%     10    10]));
%
% %% Advanced Multihit Overlapping  test 9 - partially overlapping colinear 1
% fprintf(1,'Partially overlapping colinear  \n');
% path = [0 10; 10 10; 12 8; 14 10; 15 10];
% sensor_vector_start = [-2 10];
% sensor_vector_end   = [5 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segments] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printMoreResults(distance,location,path_segments);
%
% assert(isequal(round(distance,4),[    2.0000
%     7.0000]));
% assert(isequal(round(location,4),[         0   10.0000
%     5.0000   10.0000]));
%
% %% Advanced Multihit Overlapping  test 10 - partially overlapping colinear 1
% fprintf(1,'Partially overlapping colinear  \n');
% path = [0 10; 10 10; 12 8; 14 10; 15 10];
% sensor_vector_start = [3 10];
% sensor_vector_end   = [5 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segments] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printMoreResults(distance,location,path_segments);
%
% assert(isequal(round(distance,4),[0; 2]));
% assert(isequal(round(location,4),[     3    10
%      5    10]));
%
% %% Advanced Multihit Overlapping  test 11 - partially overlapping colinear 1
% fprintf(1,'Partially overlapping colinear  \n');
% path = [0 10; 10 10; 12 8; 14 10; 15 10];
% sensor_vector_start = [3 10];
% sensor_vector_end   = [15 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segments] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printMoreResults(distance,location,path_segments);
%
% assert(isequal(round(distance,4),[     0
%      7
%     11
%     11
%      7
%     12]));
% assert(isequal(round(location,4),[     3    10
%     10    10
%     14    10
%     14    10
%     10    10
%     15    10]));
%
% %% Advanced Multihit Overlapping  test 12 - super overlapping colinear 1
% fprintf(1,'Super overlapping colinear  \n');
% path = [0 10; 10 10; 12 8; 14 10; 15 10];
% sensor_vector_start = [-3 10];
% sensor_vector_end   = [15 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segments] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printMoreResults(distance,location,path_segments);
%
% assert(isequal(round(distance,4),[     3
%     13
%     17
%     17
%     13
%     18]));
% assert(isequal(round(location,4),[     0    10
%     10    10
%     14    10
%     14    10
%     10    10
%     15    10]));
%
% %% Advanced Multihit Overlapping  test 13 - end overlapping colinear 1
% fprintf(1,'End overlapping colinear  \n');
% path = [0 10; 10 10; 12 8; 14 10; 15 10];
% sensor_vector_start = [-3 10];
% sensor_vector_end   = [0 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segments] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printMoreResults(distance,location,path_segments);
%
% assert(isequal(round(distance,4),3));
% assert(isequal(round(location,4),[0 10]));
%
% %% Advanced Multihit Overlapping  test 14 - end overlapping colinear 2
% fprintf(1,'End overlapping colinear  \n');
% path = [0 10; 10 10; 12 8; 14 10; 15 10];
% sensor_vector_start = [10 10];
% sensor_vector_end   = [15 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segments] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printMoreResults(distance,location,path_segments);
%
% assert(isequal(round(distance,4),[     0
%      0
%      4
%      4
%      5]));
% assert(isequal(round(location,4),[    10    10
%     10    10
%     14    10
%     14    10
%     15    10]));
%
% %% Advanced Multihit Overlapping  test 27 - identically overlapping colinear BACKWARDS
% fprintf(1,'identically overlapping colinear BACKWARDS  \n');
% path = [0 10; 10 10; 12 8; 14 10; 15 10];
% sensor_vector_end = [0 10];
% sensor_vector_start   = [10 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segments] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printMoreResults(distance,location,path_segments);
%
% assert(isequal(round(distance,4),[          0
%      0
%     10]));
% assert(isequal(round(location,4),[    10    10
%     10    10
%      0    10]));
%
% %% Advanced Multihit Overlapping  test 28 - partially overlapping colinear 1 BACKWARDS
% fprintf(1,'Partially overlapping colinear BACKWARDS  \n');
% path = [0 10; 10 10; 12 8; 14 10; 15 10];
% sensor_vector_end = [-2 10];
% sensor_vector_start   = [10 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segments] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printMoreResults(distance,location,path_segments);
%
% assert(isequal(round(distance,4),[     0
%      0
%     10]));
% assert(isequal(round(location,4),[    10    10
%     10    10
%      0    10]));
%
% %% Advanced Multihit Overlapping  test 29 - partially overlapping colinear 1 BACKWARDS
% fprintf(1,'Partially overlapping colinear BACKWARDS  \n');
% path = [0 10; 10 10; 12 8; 14 10; 15 10];
% sensor_vector_end = [-2 10];
% sensor_vector_start   = [5 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segments] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printMoreResults(distance,location,path_segments);
%
% assert(isequal(round(distance,4),[     0
%      5]));
% assert(isequal(round(location,4),[
%      5    10
%      0    10]));
%
% %% Advanced Multihit Overlapping  test 30 - partially overlapping colinear 1 BACKWARDS
% fprintf(1,'Partially overlapping colinear BACKWARDS  \n');
% path = [0 10; 10 10; 12 8; 14 10; 15 10];
% sensor_vector_end = [3 10];
% sensor_vector_start   = [5 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segments] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printMoreResults(distance,location,path_segments);
%
% assert(isequal(round(distance,4),[0; 2]));
% assert(isequal(round(location,4),[     5    10
%      3    10]));
%
% %% Advanced Multihit Overlapping  test 31 - partially overlapping colinear 1 BACKWARDS
% fprintf(1,'Partially overlapping colinear BACKWARDS  \n');
% path = [0 10; 10 10; 12 8; 14 10; 15 10];
% sensor_vector_end = [3 10];
% sensor_vector_start   = [15 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segments] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printMoreResults(distance,location,path_segments);
%
% assert(isequal(round(distance,4),[     5
%      5
%      1
%      0
%     12
%      1]));
% assert(isequal(round(location,4),[   10.0000   10.0000
%    10.0000   10.0000
%    14.0000   10.0000
%    15.0000   10.0000
%     3.0000   10.0000
%    14.0000   10.0000]));
%
% %% Advanced Multihit Overlapping  test 32 - super overlapping colinear 1 BACKWARDS
% fprintf(1,'Super overlapping colinear BACKWARDS  \n');
% path = [0 10; 10 10; 12 8; 14 10; 15 10];
% sensor_vector_end = [-3 10];
% sensor_vector_start   = [15 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segments] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printMoreResults(distance,location,path_segments);
%
% assert(isequal(round(distance,4),[     5
%      5
%      1
%      0
%     15
%      1]));
% assert(isequal(round(location,4),[    10    10
%     10    10
%     14    10
%     15    10
%      0    10
%     14    10]));
%
%
% %% Advanced Multihit Overlapping  test 33 - end overlapping colinear 1 BACKWARDS
% fprintf(1,'End overlapping colinear  \n');
% path = [0 10; 10 10; 12 8; 14 10; 15 10];
% sensor_vector_end = [-3 10];
% sensor_vector_start   = [0 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segments] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printMoreResults(distance,location,path_segments);
%
% assert(isequal(round(distance,4),0));
% assert(isequal(round(location,4),[0 10]));
%
% %% Advanced Multihit Overlapping  test 34 - end overlapping colinear 2 BACKWARDS
% fprintf(1,'End overlapping colinear BACKWARDS  \n');
% path = [0 10; 10 10; 12 8; 14 10; 15 10];
% sensor_vector_end = [10 10];
% sensor_vector_start   = [15 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segments] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printMoreResults(distance,location,path_segments);
%
% assert(isequal(round(distance,4),[     5
%      5
%      1
%      0
%      1]));
% assert(isequal(round(location,4),[    10    10
%     10    10
%     14    10
%     15    10
%     14    10]));
%
% %% Advanced Multihit Overlapping  test 15 - non overlapping colinear 1
% fprintf(1,'Non overlapping colinear  \n');
% path = [0 10; 10 10; 12 8; 14 10; 15 10];
% sensor_vector_start = [-3 10];
% sensor_vector_end   = [-1 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segments] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printMoreResults(distance,location,path_segments);
%
% assert(isempty(distance));
% assert(isempty(location));
%
% %% Advanced Multihit Overlapping  test 15 - partially non-overlapping colinear 2
% fprintf(1,'Non overlapping colinear  \n');
% path = [0 10; 10 10; 12 8; 14 10; 15 10];
% sensor_vector_start = [13 10];
% sensor_vector_end   = [15 10];
% fig_debugging = 2343;
% flag_search_type = 2;
% [distance,location,path_segments] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_debugging);
% fcn_INTERNAL_printMoreResults(distance,location,path_segments);
%
% assert(isequal(round(distance,4),[     1
%      1
%      2]));
% assert(isequal(round(location,4),[    14    10
%     14    10
%     15    10]));
%
%
% %% BUG cases
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  ____  _    _  _____
% % |  _ \| |  | |/ ____|
% % | |_) | |  | | |  __    ___ __ _ ___  ___  ___
% % |  _ <| |  | | | |_ |  / __/ _` / __|/ _ \/ __|
% % | |_) | |__| | |__| | | (_| (_| \__ \  __/\__ \
% % |____/ \____/ \_____|  \___\__,_|___/\___||___/
% %
% % See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=BUG%20cases
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % All bug case figures start with the number 9
%
% close all;
%
% %% BUG - infinite point intersection (4), flag (1), test 1 - should be infinite intersections
%
% fig_num = 91001;
% figure(fig_num); clf;
%
% fprintf(1,'\n BUG - infinite point intersection (4), flag (1), test 1 - should be infinite intersections  \n');
%
% path = [0 10; 10 10];
% sensor_vector_start = [13 10];
% sensor_vector_end   = [11 10];
% flag_search_type = 1;
%
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,fig_num);
% fcn_INTERNAL_printResults(distance,location);
%
%
% % Check variable types
% assert(isnumeric(distance));
% assert(isnumeric(location));
%
% % Check variable sizes
% assert(isequal(size(distance),[1 1]));
% assert(isequal(size(location),[1 2]));
%
% % Check variable values
% % assert(isequal(round(distance,4),9.2043));
% % assert(isequal(round(location,4),[3.9286,10.0000]));
%
% % Make sure plot opened up
% assert(isequal(get(gcf,'Number'),fig_num));
%
%
%
%
%
% %% Fast Mode Tests
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  ______        _     __  __           _        _______        _
% % |  ____|      | |   |  \/  |         | |      |__   __|      | |
% % | |__ __ _ ___| |_  | \  / | ___   __| | ___     | | ___  ___| |_ ___
% % |  __/ _` / __| __| | |\/| |/ _ \ / _` |/ _ \    | |/ _ \/ __| __/ __|
% % | | | (_| \__ \ |_  | |  | | (_) | (_| |  __/    | |  __/\__ \ |_\__ \
% % |_|  \__,_|___/\__| |_|  |_|\___/ \__,_|\___|    |_|\___||___/\__|___/
% %
% %
% % See: http://patorjk.com/software/taag/#p=display&f=Big&t=Fast%20Mode%20Tests
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Figures start with 8
%
% close all;
%
% %% Basic example - NO FIGURE
%
% fig_num = 80001;
% figure(fig_num);
% close(fig_num);
%
% fprintf(1,'\nSingle point intersection (2), flag (0), test 1  \n');
%
% path = [0 10; 10 10];
% sensor_vector_start = [2 1];
% sensor_vector_end   = [5 15];
% flag_search_type = 0;
%
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,[]);
% fcn_INTERNAL_printResults(distance,location);
%
%
% % Check variable types
% assert(isnumeric(distance));
% assert(isnumeric(location));
%
% % Check variable sizes
% assert(isequal(size(distance),[1 1]));
% assert(isequal(size(location),[1 2]));
%
% % Check variable values
% assert(isequal(round(distance,4),9.2043));
% assert(isequal(round(location,4),[3.9286,10.0000]));
%
% % Make sure plot did NOT open up
% figHandles = get(groot, 'Children');
% assert(~any(figHandles==fig_num));
%
% %% Basic fast mode - NO FIGURE, FAST MODE
% fig_num = 80002;
% figure(fig_num);
% close(fig_num);
%
% fprintf(1,'\nSingle point intersection (2), flag (0), test 1  \n');
%
% path = [0 10; 10 10];
% sensor_vector_start = [2 1];
% sensor_vector_end   = [5 15];
% flag_search_type = 0;
%
% [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,-1);
% fcn_INTERNAL_printResults(distance,location);
%
%
% % Check variable types
% assert(isnumeric(distance));
% assert(isnumeric(location));
%
% % Check variable sizes
% assert(isequal(size(distance),[1 1]));
% assert(isequal(size(location),[1 2]));
%
% % Check variable values
% assert(isequal(round(distance,4),9.2043));
% assert(isequal(round(location,4),[3.9286,10.0000]));
%
% % Make sure plot did NOT open up
% figHandles = get(groot, 'Children');
% assert(~any(figHandles==fig_num));
%
%
% %% Compare speeds of pre-calculation versus post-calculation versus a fast variant
% fig_num = 80003;
% figure(fig_num);
% close(fig_num);
%
% fprintf(1,'\nSingle point intersection (2), flag (0), test 1  \n');
%
% path = [0 10; 10 10];
% sensor_vector_start = [2 1];
% sensor_vector_end   = [5 15];
% flag_search_type = 0;
%
% Niterations = 100;
%
% % Do calculation without pre-calculation
% tic;
% for ith_test = 1:Niterations
%     % Call the function
%     [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,[]);
% end
% slow_method = toc;
%
% % Do calculation with pre-calculation, FAST_MODE on
% tic;
% for ith_test = 1:Niterations
%     % Call the function
%     [distance,location,path_segment, t, u] = ...
%     fcn_Path_findProjectionHitOntoPath(...
%     path,sensor_vector_start,sensor_vector_end,...
%     flag_search_type,-1);
%
% end
% fast_method = toc;
%
% % Plot results as bar chart
% figure(373737);
% clf;
% X = categorical({'Normal mode','Fast mode'});
% X = reordercats(X,{'Normal mode','Fast mode'}); % Forces bars to appear in this exact order, not alphabetized
% Y = [slow_method fast_method ]*1000/Niterations;
% bar(X,Y)
% ylabel('Execution time (Milliseconds)')
%
%
% % Make sure plot did NOT open up
% figHandles = get(groot, 'Children');
% assert(~any(figHandles==fig_num));

%% fcn_INTERNAL_printResults
function fcn_INTERNAL_printResults(distance,location)
fprintf(1,'Distance \t Location X \t Location Y \n');
if ~isempty(distance)
    for i_result = 1:length(distance(:,1))
        fprintf(1,'%.3f \t\t %.3f \t\t\t %.3f\n',distance(i_result),location(i_result,1),location(i_result,2));
    end
end
end % Ends fcn_INTERNAL_printResults

%% fcn_INTERNAL_printMoreResults
function fcn_INTERNAL_printMoreResults(distance,location,path_segments) %#ok<DEFNU>
fprintf(1,'\nDistance \t Location X \t Location Y \t PathSegment \n');
if ~isempty(distance)
    for i_result = 1:length(distance(:,1))
        fprintf(1,'%.3f \t\t %.3f \t\t\t %.3f \t\t %.0d\n',distance(i_result),location(i_result,1),location(i_result,2),path_segments(i_result));
    end
end
end % Ends fcn_INTERNAL_printMoreResults

%% fcn_INTERNAL_checkTestCases
function fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)

% Get all the fields
expectedFields = fieldnames(expected);
Nfields = length(expectedFields);

% Loop through expected fields, checking if they are in actual
for ith_field = 1:Nfields
    thisField = expectedFields{ith_field};

    if ~isstruct(expected.(thisField))
        % Inheret values
        expectedValue  = expected.(thisField);
        expectedSize   = size(expectedValue);
        expectedValueDigits = 4;

        % Find inhereted TypeString
        typesToTest = {'numeric','logical','char','string', 'struct','handle'};
        expectedType = 'nan'; % Default is nan
        for ithTypeToCheck = 1:length(typesToTest)
            if isa(expectedValue,typesToTest{ithTypeToCheck})
                expectedType = typesToTest{ithTypeToCheck};
            end
        end

    else
        expectedType  = expected.(thisField).TypeString;
        expectedSize  = expected.(thisField).Size;
        expectedValue = expected.(thisField).Value;
        expectedValueDigits = expected.(thisField).ValueDigits;
    end

    thisVariable = actual.(thisField);
    % Check variable type
    if ~isa(thisVariable,expectedType)|| ~isequal(size(thisVariable),expectedSize) | ~isequal(round(thisVariable,expectedValueDigits),expectedValue)
        fprintf(1,'Error will be thrown for variable: %s:\n',thisField);
        fprintf(1,'Expected type of: %s, with size: [%.0d %.0d]\n',expectedType, expectedSize(1,1), expectedSize(1,2));
        fprintf(1,'Expected value:\n');
        disp(round(expectedValue,expectedValueDigits));
        fprintf(1,'Actual value:\n');
        disp(round(thisVariable,expectedValueDigits));


    end
    assert(isa(thisVariable,expectedType));
    assert(isequal(size(thisVariable),expectedSize));
    assert(isequal(round(thisVariable,expectedValueDigits),expectedValue));

end

% Make sure plot opened up
if isfield(plotting,'FigureExpected') && plotting.FigureExpected==1
    assert(isequal(get(gcf,'Number'),inputs.fig_num));
end

end