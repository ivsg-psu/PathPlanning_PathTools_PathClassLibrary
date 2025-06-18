% script_test_fcn_Path_findProjectionHitOntoPathImproved
% This is a script to exercise the function: fcn_Path_findProjectionHitOntoPathImproved.m
% This function was written on 2025_06_15 by S. Brennan, from original
% version (started 2020_12_31)
% Questions or comments? sbrennan@psu.edu

% Modification history:
% 2020_12_31 - S. Brennan
% -- Updated for new release
% -- Organized sections
% -- Added fast mode tests
% -- Added full assertion tests on output variables
% -- Added figure open/close assertions
% -- Added automated assertion testing

close all

%%%%%%%%%%%%
% FIGURE NUMBERING:
% FSTXXX
%
% F is First figure number, starting with:
% 1: demonstration cases
% 2: single point intersection cases with one wall
% 3: infinite intersection cases with one wall
% 4: multi-hit cases
% 6: multi-hit overlapping cases
% 7: multi-hit multi-wall cases
% 8: fast mode cases
% 9: known bug cases
%
% S is second figure number, flag_search_return_type:
% 0: first intersection if there is any overlap
% 1: all intersections, if there is overlap
%
% T is third figure number, flag_search_range_type:
% 0: (default) the GIVEN sensor and GIVEN wall used.
% 1: ANY projection of the sensor used with the GIVEN wall
% 2: ANY projection of the wall used with the GIVEN sensor
% 3: ANY projection of BOTH the wall and sensor
%
% t: is tolerance
% 0: uses defaults
% 1: uses positive tolerances
% 2: uses negative tolerances
%
% XX: 4th to 6th number: a counter that counts up through the cases in this
% section.
%
% Example:
% 213206 plots a single point intersection, flag_search_return_type type 4, flag_search_range_type of 3, negative tolerances, 6th test case


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
inputs.wall_start = [0 10];
inputs.wall_end   = [10 10];
inputs.sensor_vector_start = [5 0];
inputs.sensor_vector_end   = [5 15];
inputs.flag_search_return_type = 0;
inputs.flag_search_range_type = 0;
inputs.tolerance = [];

% SHORT format checking
clear expected
expected.distance = 10;
expected.location = [5 10];
expected.path_segment = 1;
expected.t = 0.5;
expected.u = 0.6667;

actual = struct;
[actual.distance, actual.location, actual.path_segment, actual.t, actual.u] = ...
    fcn_Path_findProjectionHitOntoPathImproved(...
    inputs.wall_start, inputs.wall_end,...
    inputs.sensor_vector_start,inputs.sensor_vector_end,...
    (inputs.flag_search_return_type), (inputs.flag_search_range_type), ...
    (inputs.tolerance), (inputs.fig_num));

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
inputs.wall_start = [0 10];
inputs.wall_end   = [10 10];

inputs.sensor_vector_start = [5 0];
inputs.sensor_vector_end   = [5 15];
inputs.flag_search_return_type = 0;
inputs.flag_search_range_type = 0;
inputs.tolerance = [];

NwallSegments = length(inputs.wall_start(:,1));

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

expected.wall_segment.TypeString  = 'numeric'; % See 'help isa' for a listing
expected.wall_segment.Size        = [Nsolutions 1];
expected.wall_segment.Value       = 1;
expected.wall_segment.ValueDigits = 0; % How many digits must match

expected.t.TypeString  = 'numeric'; % See 'help isa' for a listing
expected.t.Size        = [NwallSegments 1]; % Expecting 1 solution
expected.t.Value       = 0.5;
expected.t.ValueDigits = 4; % How many digits must match

expected.u.TypeString  = 'numeric'; % See 'help isa' for a listing
expected.u.Size        = [NwallSegments 1]; % Expecting 1 solution
expected.u.Value       = 0.6667;
expected.u.ValueDigits = 4; % How many digits must match

actual = struct;
[actual.distance, actual.location, actual.wall_segment, actual.t, actual.u] = ...
    fcn_Path_findProjectionHitOntoPathImproved(...
    inputs.wall_start, inputs.wall_end,...
    inputs.sensor_vector_start,inputs.sensor_vector_end,...
    (inputs.flag_search_return_type), (inputs.flag_search_range_type), ...
    (inputs.tolerance), (inputs.fig_num));

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

close all;

%% Single point intersection 2XXX1 cases
% flag_search_return_type 0: first intersection if there is any overlap
% flag_search_range_type  0: (default) the GIVEN sensor and GIVEN wall used.
% tolerance 0: uses defaults

intersectionTestType = 2;
return_flags = [0 1];
range_flags  = [0 1 2 3];
tolerance_flags = [0 1 2];
thisCase = 1;

for ith_return = 1:length(return_flags)
    thisReturn = return_flags(ith_return);
    for ith_range = 1:length(range_flags)
        thisRange = range_flags(ith_range);
        for ith_tolerance = 1:length(tolerance_flags)
            thisTolerance = tolerance_flags(ith_tolerance);

            % Build the figure number
            fig_string = sprintf('%.0f%.0f%.0f%.0f%.0f',intersectionTestType, thisReturn,thisRange,thisTolerance,thisCase);
            fig_num = str2double(fig_string);
            plotting.FigureExpected = 1;

            % Build the test cases
            testCases = fcn_INTERNAL_fillTestCasesVerticalArrowSensors(fig_num);

            for ith_testCase = 1:length(testCases)

                % Plot the sensor vector in red, to highlight which one we're on
                quiver(testCases(ith_testCase).sensor_vector_start(1,1),testCases(ith_testCase).sensor_vector_start(1,2),0,1,'r','Linewidth',1,'MaxHeadSize',1);


                clear inputs
                inputs = testCases(ith_testCase);

                % SHORT format checking
                clear expected
                expected.distance = testCases(ith_testCase).expected.Distance;
                expected.location = testCases(ith_testCase).expected.Intersection;
                expected.wall_segment = testCases(ith_testCase).expected.wall_segment;
                expected.t = testCases(ith_testCase).expected.t;
                expected.u = testCases(ith_testCase).expected.u;

                actual = struct;
                [actual.distance, actual.location, actual.wall_segment, actual.t, actual.u] = ...
                    fcn_Path_findProjectionHitOntoPathImproved(...
                    inputs.wall_start, inputs.wall_end,...
                    inputs.sensor_vector_start,inputs.sensor_vector_end,...
                    (inputs.flag_search_return_type), (inputs.flag_search_range_type), ...
                    (inputs.tolerance), (inputs.fig_num));

                fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)
                %fcn_INTERNAL_printResults(actual.distance,actual.location);

            end % Ends loop through cases

            close all;

        end % Ends loop through tolerance flags
    end % Ends loop through range flags
end % Ends loop throough return flags


%% Infinite intersection cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  _____        __ _       _ _         _____       _                          _   _
% |_   _|      / _(_)     (_) |       |_   _|     | |                        | | (_)
%   | |  _ __ | |_ _ _ __  _| |_ ___    | |  _ __ | |_ ___ _ __ ___  ___  ___| |_ _  ___  _ __  ___
%   | | | '_ \|  _| | '_ \| | __/ _ \   | | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __| |/ _ \| '_ \/ __|
%  _| |_| | | | | | | | | | | ||  __/  _| |_| | | | ||  __/ |  \__ \  __/ (__| |_| | (_) | | | \__ \
% |_____|_| |_|_| |_|_| |_|_|\__\___| |_____|_| |_|\__\___|_|  |___/\___|\___|\__|_|\___/|_| |_|___/
%
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Infinite%20Intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All infinte intersection figures start with the number 3

close all;

%% Infinite intersection 3XXX1 cases
% flag_search_return_type 0: first intersection if there is any overlap
% flag_search_range_type  0: (default) the GIVEN sensor and GIVEN wall used.
% tolerance 0: uses defaults

% URHERE

intersectionTestType = 3;
return_flags = 0; % [0 1];
range_flags  = 1;  % [0 1 2 3];
tolerance_flags = 0; %[0 1 2];
thisCase = 1;

for ith_return = 1:length(return_flags)
    thisReturn = return_flags(ith_return);
    for ith_range = 1:length(range_flags)
        thisRange = range_flags(ith_range);
        for ith_tolerance = 1:length(tolerance_flags)
            thisTolerance = tolerance_flags(ith_tolerance);

            % Build the figure number
            fig_string = sprintf('%.0f%.0f%.0f%.0f%.0f',intersectionTestType, thisReturn,thisRange,thisTolerance,thisCase);
            fig_num = str2double(fig_string);
            plotting.FigureExpected = 1;

            % Build the test cases
            testCases = fcn_INTERNAL_fillTestCasesVerticalArrowSensors(fig_num);

            for ith_testCase = 5:length(testCases)

                % Plot the sensor vector in red, to highlight which one we're on
                vectorX = testCases(ith_testCase).sensor_vector_end(1,1) - testCases(ith_testCase).sensor_vector_start(1,1);
                vectorY = testCases(ith_testCase).sensor_vector_end(1,2) - testCases(ith_testCase).sensor_vector_start(1,2);
                quiver(testCases(ith_testCase).sensor_vector_start(1,1),testCases(ith_testCase).sensor_vector_start(1,2),vectorX,vectorY,'r','Linewidth',1,'MaxHeadSize',1);


                clear inputs
                inputs = testCases(ith_testCase);

                % SHORT format checking
                clear expected
                expected.distance = testCases(ith_testCase).expected.Distance;
                expected.location = testCases(ith_testCase).expected.Intersection;
                expected.wall_segment = testCases(ith_testCase).expected.wall_segment;
                expected.t = testCases(ith_testCase).expected.t;
                expected.u = testCases(ith_testCase).expected.u;

                actual = struct;
                [actual.distance, actual.location, actual.wall_segment, actual.t, actual.u] = ...
                    fcn_Path_findProjectionHitOntoPathImproved(...
                    inputs.wall_start, inputs.wall_end,...
                    inputs.sensor_vector_start,inputs.sensor_vector_end,...
                    (inputs.flag_search_return_type), (inputs.flag_search_range_type), ...
                    (inputs.tolerance), (inputs.fig_num));

                fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)
                %fcn_INTERNAL_printResults(actual.distance,actual.location);

            end % Ends loop through cases

            close all;

        end % Ends loop through tolerance flags
    end % Ends loop through range flags
end % Ends loop throough return flags

%%


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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
%     fcn_Path_findProjectionHitOntoPathImproved(...
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
function fcn_INTERNAL_printResults(distance,location) %#ok<DEFNU>
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
    flag_errorWillBeThrown = 0;
    if ~isa(thisVariable,expectedType)|| ~isequal(size(thisVariable),expectedSize)
        flag_errorWillBeThrown = 1;
    end
    if ~any(isnan(thisVariable),'all') && ~isequal(round(thisVariable,expectedValueDigits),expectedValue)
        flag_errorWillBeThrown = 1;
    end
    try
        if ~any(isnan(expectedValue),'all') && ~isequal(round(thisVariable,expectedValueDigits),expectedValue)
            flag_errorWillBeThrown = 1;
        end
    catch
        disp('stop here');
    end
    if all(isnan(thisVariable)) && ~all(isnan(thisVariable))
        flag_errorWillBeThrown = 1;
    end
    if 1==flag_errorWillBeThrown
        fprintf(1,'Error will be thrown for variable: %s:\n',thisField);
        fprintf(1,'Expected type of: %s, with size: [%.0d %.0d]\n',expectedType, expectedSize(1,1), expectedSize(1,2));
        fprintf(1,'Expected value:\n');
        disp(round(expectedValue,expectedValueDigits));
        fprintf(1,'Actual value:\n');
        disp(round(thisVariable,expectedValueDigits));
    end
    assert(isa(thisVariable,expectedType));
    assert(isequal(size(thisVariable),expectedSize));
    if all(isnan(expectedValue))
        assert(all(isnan(thisVariable)));
    else
        assert(isequal(round(thisVariable,expectedValueDigits),expectedValue));
    end

end

% Make sure plot opened up
if isfield(plotting,'FigureExpected') && plotting.FigureExpected==1
    assert(isequal(get(gcf,'Number'),inputs.fig_num));
end

end

%% fcn_INTERNAL_plotTestCases
function fcn_INTERNAL_plotTestCases(fig_num, wall_starts, wall_ends, testCases)
% Plot all the test case results
figure(fig_num);
clf; hold on;
axis equal;
grid on;
grid minor;

% Plot the walls in black
Nwalls = length(wall_starts(:,1));
allWallsX = [wall_starts(:,1) wall_ends(:,1) nan(Nwalls,1)];
allWallsX = reshape(allWallsX',1,[]);
allWallsY = [wall_starts(:,2) wall_ends(:,2) nan(Nwalls,1)];
allWallsY = reshape(allWallsY',1,[]);
allWalls = [allWallsX' allWallsY'];

plot(allWalls(:,1),allWalls(:,2),'k.-','Linewidth',2);

for ith_testCase = 1:length(testCases)
    handle_text = text(allWalls(1,1),allWalls(1,2),'Walls');
    set(handle_text,'Color',[0 0 0]);

    % Plot the sensor vector
    vectorX = testCases(ith_testCase).sensor_vector_end(1,1) - testCases(ith_testCase).sensor_vector_start(1,1);
    vectorY = testCases(ith_testCase).sensor_vector_end(1,2) - testCases(ith_testCase).sensor_vector_start(1,2);
    quiver(testCases(ith_testCase).sensor_vector_start(1,1),testCases(ith_testCase).sensor_vector_start(1,2),vectorX,vectorY,'b','Linewidth',1,'MaxHeadSize',1);

    % Plot the predicted solution
    solution = testCases(ith_testCase).expected.Intersection;
    if ~all(isnan(solution))
        plot(solution(:,1),solution(:,2),'r.','MarkerSize',10);
    end
end
end % Ends fcn_INTERNAL_plotTestCases

function trimmed = fcn_INTERNAL_trimEnds(vector)

if length(vector)>2
    trimmed = vector(2:end-1);
else
    trimmed = vector;
end
end

%% fcn_INTERNAL_fillTestCases
function testCases = fcn_INTERNAL_fillTestCasesVerticalArrowSensors(fig_num)

temp = sprintf('%.0d',fig_num);


%%%%%%%%%%%%
% FIGURE NUMBERING:
% FSTXXX
firstFigureNumber = str2double(temp(1));
secondFigureNumber = str2double(temp(2));
thirdFigureNumber = str2double(temp(3));
fourthFigureNumber = str2double(temp(4));


% F is First figure number, starting with:
% 1: demonstration cases
% 2: single point intersection cases with one wall
% 3: infinite intersection cases with one wall


% S is second figure number, flag_search_return_type:
% 0: first intersection if there is any overlap
% 1: all intersections, if there is overlap
flag_search_return_type = secondFigureNumber;

% T is third figure number, flag_search_range_type:
% 0: (default) the GIVEN sensor and GIVEN path used.
% 1: ANY projection of the sensor used with the GIVEN path
% 2: ANY projection of the path used with the GIVEN sensor
% 3: ANY projection of BOTH the path and sensor
% The yStarts represents the y value that each sensor vector will start at
% The xStarts represents the x value that each sensor vector will start at
%

%%%%%%
% Set all locations based on first number

switch firstFigureNumber
    case {2} % Single intersect cases
        % All are tested with single wall
        wall_starts          = [0 0];
        wall_ends            = [1 0];
        sensorLengthX = 0;
        sensorLengthY = 1;
        totalRangeX = [-1 0 0.5 1 2];
        totalRangeY = [-2 -1 -0.5 0 1];
        flag_testTolerance = 1;

        Nintersections = 1;

        % If ANY projections are used, the end points of the range will have
        % intersections. But if GIVEN projections are used, only the values from
        % 2:end-1 will have intersections. The trimEnds function gives values from
        % 2:end-1 for vectors that have length greater than 2.
        flag_search_range_type = thirdFigureNumber;
        switch thirdFigureNumber
            case {0}
                xStartsWithIntersections = fcn_INTERNAL_trimEnds(totalRangeX);
                yStartsWithIntersections = fcn_INTERNAL_trimEnds(totalRangeY);
            case {1}
                xStartsWithIntersections = fcn_INTERNAL_trimEnds(totalRangeX);
                yStartsWithIntersections = totalRangeY;
            case {2}
                xStartsWithIntersections = totalRangeX;
                yStartsWithIntersections = fcn_INTERNAL_trimEnds(totalRangeY);
            case {3}
                xStartsWithIntersections = totalRangeX;
                yStartsWithIntersections = totalRangeY;
            otherwise
                warning('on','backtrace');
                warning('Expecting a third figure number with integer values of 0 to 3, but found: %.3f',thirdFigureNumber);
                error('Third figure number not recognized');
        end

    case {3} % Infinite intersect cases
        % All are tested with single wall
        wall_starts          = [0 0];
        wall_ends            = [1 0];
        sensorLengthX = 1;
        sensorLengthY = 0;
        totalRangeX = [-1.5 -0.5 0 0.5 1.5];
        totalRangeY = 0;
        flag_testTolerance = 0;

        % If ANY projections are used, the end points of the range will have
        % intersections. But if GIVEN projections are used, only the values from
        % 2:end-1 will have intersections. The trimEnds function gives values from
        % 2:end-1 for vectors that have length greater than 2.
        flag_search_range_type = thirdFigureNumber;
        switch thirdFigureNumber
            case {0}
                xStartsWithIntersections = fcn_INTERNAL_trimEnds(totalRangeX);
                yStartsWithIntersections = fcn_INTERNAL_trimEnds(totalRangeY);
            case {1,2,3}
                xStartsWithIntersections = totalRangeX;
                yStartsWithIntersections = totalRangeY;
            otherwise
                warning('on','backtrace');
                warning('Expecting a third figure number with integer values of 0 to 3, but found: %.3f',thirdFigureNumber);
                error('Third figure number not recognized');
        end

    case {5,6,7} % All multi-hit cases
        % All are tested with single wall
        wall_starts          = [0 0; 1 0;];
        wall_ends            = [1 0; 1 1];
        switch secondFigureNumber
            case {0} % First intersection only
                % Nintersections = 1;
            case {1} % All intersections
                % Nintersections = 1;
            otherwise
                warning('on','backtrace');
                warning('Expecting a second figure number with integer values of 0 to 1, but found: %.3f',secondFigureNumber);
                error('Second figure number not recognized');

        end
    otherwise
        warning('on','backtrace');
        warning('Expecting a first figure number with integer values of 2 to 7, but found: %.3f',firstFigureNumber);
        error('Third figure number not recognized');
end





% t: is tolerance
% 0: uses defaults
% 1: uses positive tolerances
% 2: uses negative tolerances
switch fourthFigureNumber
    case {0}
        % do nothing
        tolerance = [];
    case {1}
        % Positive tolerance
        tolerance = 0.001;
    case {2}
        % Negative tolerance
        tolerance = -0.001;
        % xStarts are trimmed?
        if flag_testTolerance == 1
            if isequal(xStartsWithIntersections, [0 0.5 1])
                xStartsWithIntersections = 0.5;
            end
            if isequal(yStartsWithIntersections,[-1 -0.5 0])
                yStartsWithIntersections = -0.5;
            end
        end
    otherwise
        warning('on','backtrace');
        warning('Expecting a fourth figure number with integer values of 0 to 2, but found: %.3f',fourthFigureNumber);
        error('Fourth figure number not recognized');
end

ith_case = 0;
testCases = struct;

all_yStarts = totalRangeY;
all_xStarts = totalRangeX;


% Fill in all the test cases
for jth_yStart = 1:length(all_yStarts)
    thisY = all_yStarts(jth_yStart);
    for ith_xStart = 1:length(all_xStarts)
        thisX = all_xStarts(ith_xStart);
        ith_case = ith_case+1;
        % FOR DEBUGGING
        % if ith_case == 5
        %     disp('stop here');
        % end
        testCases(ith_case).wall_start          = wall_starts;
        testCases(ith_case).wall_end            = wall_ends;
        testCases(ith_case).sensor_vector_start = [thisX thisY];
        testCases(ith_case).sensor_vector_end   = [thisX+sensorLengthX thisY+sensorLengthY];
        testCases(ith_case).flag_search_return_type = flag_search_return_type;
        testCases(ith_case).flag_search_range_type = flag_search_range_type;
        testCases(ith_case).tolerance = tolerance;
        testCases(ith_case).fig_num = fig_num;

        % Fill in expected values
        switch firstFigureNumber
            case {2} % Single intersect cases
                if ismember(thisY,yStartsWithIntersections) && ismember(thisX, xStartsWithIntersections)
                    testCases(ith_case).expected.Intersection    = [thisX 0];
                    testCases(ith_case).expected.Distance    = 0 - thisY;
                    testCases(ith_case).expected.wall_segment = 1;
                else
                    testCases(ith_case).expected.Intersection    = [nan nan];
                    testCases(ith_case).expected.Distance    = nan;
                    testCases(ith_case).expected.wall_segment = nan;
                end
                testCases(ith_case).expected.t     = thisX;
                testCases(ith_case).expected.u     = -thisY;

            case {3} % Infinite intersect cases
                % If tolerance is negative, the lines NEVER intersect
                % (fourthFigureNumber == 2)
                if ismember(thisX, xStartsWithIntersections) && fourthFigureNumber~=2
                    maxXstart = max(thisX,0);
                    minXend   = min(thisX+1,1);
                    if thisX<0
                        distanceStart = -thisX;
                    else
                        distanceStart = 0;
                    end
                    if thirdFigureNumber==0 || secondFigureNumber==0
                        testCases(ith_case).expected.Intersection    = [maxXstart 0];
                    else
                        testCases(ith_case).expected.Intersection    = [maxXstart 0; minXend 0];
                    end
                    testCases(ith_case).expected.Distance        = distanceStart;
                    testCases(ith_case).expected.wall_segment = 1;
                    testCases(ith_case).expected.t     = [maxXstart; minXend];
                    testCases(ith_case).expected.u     = [distanceStart; minXend-thisX];

                else
                    testCases(ith_case).expected.Intersection    = [nan nan];
                    testCases(ith_case).expected.Distance    = nan;
                    testCases(ith_case).expected.wall_segment = nan;
                    testCases(ith_case).expected.t     = nan;
                    testCases(ith_case).expected.u     = nan;
                end
            otherwise
                warning('on','backtrace');
                warning('Expecting a first figure number with integer values of 2 to 7, but found: %.3f',firstFigureNumber);
                error('Third figure number not recognized');
        end

    end
end

fcn_INTERNAL_plotTestCases(fig_num, wall_starts, wall_ends, testCases);

end % Ends fcn_INTERNAL_fillTestCases