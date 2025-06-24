% script_test_fcn_Path_findSensorHitOnWall
% This is a script to exercise the function: fcn_Path_findSensorHitOnWall.m
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
% FSTtXX
%
% F is First figure number, starting with:
% 1: demonstration cases
% 2: single point intersection cases with one wall
% 3: infinite intersection cases with one wall
% 4: multi-hit cases
% 5: wall as single point cases
% 6: sensor as single point cases
% 7: (unused)
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
% Fourth number, t: is tolerance
% 0: uses defaults
% 1: uses positive tolerances
% 2: uses negative tolerances
%
% XX: 5th to 6th number: a counter that counts up through the cases in this
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
    fcn_Path_findSensorHitOnWall(...
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
    fcn_Path_findSensorHitOnWall(...
    inputs.wall_start, inputs.wall_end,...
    inputs.sensor_vector_start,inputs.sensor_vector_end,...
    (inputs.flag_search_return_type), (inputs.flag_search_range_type), ...
    (inputs.tolerance), (inputs.fig_num));

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)
%fcn_INTERNAL_printResults(actual.distance,actual.location);


%% Single hit intersection cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____ _             _        _    _ _ _     _____       _                          _   _
%  / ____(_)           | |      | |  | (_) |   |_   _|     | |                        | | (_)
% | (___  _ _ __   __ _| | ___  | |__| |_| |_    | |  _ __ | |_ ___ _ __ ___  ___  ___| |_ _  ___  _ __  ___
%  \___ \| | '_ \ / _` | |/ _ \ |  __  | | __|   | | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __| |/ _ \| '_ \/ __|
%  ____) | | | | | (_| | |  __/ | |  | | | |_   _| |_| | | | ||  __/ |  \__ \  __/ (__| |_| | (_) | | | \__ \
% |_____/|_|_| |_|\__, |_|\___| |_|  |_|_|\__| |_____|_| |_|\__\___|_|  |___/\___|\___|\__|_|\___/|_| |_|___/
%                  __/ |
%                 |___/
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Single%20Hit%20Intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All single point intersection figures start with the number 2

close all;
fprintf(1,'Figure: 2XXXXXX: Demo of single intersection cases\n');

%% Single hit intersection 2XXXX cases
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
                    fcn_Path_findSensorHitOnWall(...
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
fprintf(1,'Figure: 3XXXXXX: Demo of infinite intersection cases\n');

%% Infinite intersection 3XXX1 cases
% flag_search_return_type 0: first intersection if there is any overlap
% flag_search_range_type  0: (default) the GIVEN sensor and GIVEN wall used.
% tolerance 0: uses defaults

intersectionTestType = 3;
return_flags = [0 1];
range_flags  = [0 1 2 3];
tolerance_flags = 2; % [0 1 2];
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
                    fcn_Path_findSensorHitOnWall(...
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


%% Multi-hit tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  __  __       _ _   _   _    _ _ _     _____       _                          _   _
% |  \/  |     | | | (_) | |  | (_) |   |_   _|     | |                        | | (_)
% | \  / |_   _| | |_ _  | |__| |_| |_    | |  _ __ | |_ ___ _ __ ___  ___  ___| |_ _  ___  _ __  ___
% | |\/| | | | | | __| | |  __  | | __|   | | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __| |/ _ \| '_ \/ __|
% | |  | | |_| | | |_| | | |  | | | |_   _| |_| | | | ||  __/ |  \__ \  __/ (__| |_| | (_) | | | \__ \
% |_|  |_|\__,_|_|\__|_| |_|  |_|_|\__| |_____|_| |_|\__\___|_|  |___/\___|\___|\__|_|\___/|_| |_|___/
%
%
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Multi%20Hit%20Intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All multi-hit intersection figures start with the number 4

close all;
fprintf(1,'Figure: 4XXXXXX: Demo of multi-hit intersection cases\n');

%% Multi hit intersection 4XXXXX cases

intersectionTestType = 4;
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
                    fcn_Path_findSensorHitOnWall(...
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

%% Wall as a single point tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% __          __   _ _                        _____ _             _        _____      _       _
% \ \        / /  | | |                      / ____(_)           | |      |  __ \    (_)     | |
%  \ \  /\  / /_ _| | |   __ _ ___    __ _  | (___  _ _ __   __ _| | ___  | |__) |__  _ _ __ | |_
%   \ \/  \/ / _` | | |  / _` / __|  / _` |  \___ \| | '_ \ / _` | |/ _ \ |  ___/ _ \| | '_ \| __|
%    \  /\  / (_| | | | | (_| \__ \ | (_| |  ____) | | | | | (_| | |  __/ | |  | (_) | | | | | |_
%     \/  \/ \__,_|_|_|  \__,_|___/  \__,_| |_____/|_|_| |_|\__, |_|\___| |_|   \___/|_|_| |_|\__|
%                                                            __/ |
%                                                           |___/
%   _____       _                          _   _
%  |_   _|     | |                        | | (_)
%    | |  _ __ | |_ ___ _ __ ___  ___  ___| |_ _  ___  _ __  ___
%    | | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __| |/ _ \| '_ \/ __|
%   _| |_| | | | ||  __/ |  \__ \  __/ (__| |_| | (_) | | | \__ \
%  |_____|_| |_|\__\___|_|  |___/\___|\___|\__|_|\___/|_| |_|___/
%
%
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Wall%20as%20a%20Single%20Point%0A%20Intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All wall as a single point wall intersection figures start with the number 5

close all;
fprintf(1,'Figure: 3XXXXXX: Demo of single point wall intersection cases\n');

%% Wall as a single point intersection 5XXXX cases

intersectionTestType = 5;
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
                    fcn_Path_findSensorHitOnWall(...
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


%% Sensor as a single point tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____                                                  _____ _             _        _____      _       _
%  / ____|                                                / ____(_)           | |      |  __ \    (_)     | |
% | (___   ___ _ __  ___  ___  _ __    __ _ ___    __ _  | (___  _ _ __   __ _| | ___  | |__) |__  _ _ __ | |_
%  \___ \ / _ \ '_ \/ __|/ _ \| '__|  / _` / __|  / _` |  \___ \| | '_ \ / _` | |/ _ \ |  ___/ _ \| | '_ \| __|
%  ____) |  __/ | | \__ \ (_) | |    | (_| \__ \ | (_| |  ____) | | | | | (_| | |  __/ | |  | (_) | | | | | |_
% |_____/ \___|_| |_|___/\___/|_|     \__,_|___/  \__,_| |_____/|_|_| |_|\__, |_|\___| |_|   \___/|_|_| |_|\__|
%                                                                         __/ |
%                                                                        |___/
%   _____       _                          _   _
%  |_   _|     | |                        | | (_)
%    | |  _ __ | |_ ___ _ __ ___  ___  ___| |_ _  ___  _ __  ___
%    | | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __| |/ _ \| '_ \/ __|
%   _| |_| | | | ||  __/ |  \__ \  __/ (__| |_| | (_) | | | \__ \
%  |_____|_| |_|\__\___|_|  |___/\___|\___|\__|_|\___/|_| |_|___/
%
%
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Sensor%20as%20a%20Single%20Point%0A%20Intersections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All wall as a single point sensor intersection figures start with the number 6

close all;
fprintf(1,'Figure: 3XXXXXX: Demo of single point sensor intersection cases\n');

%% Wall as a single point intersection 6XXXX cases

intersectionTestType = 6;
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
                    fcn_Path_findSensorHitOnWall(...
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

%% Fast Mode Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ______        _     __  __           _        _______        _
% |  ____|      | |   |  \/  |         | |      |__   __|      | |
% | |__ __ _ ___| |_  | \  / | ___   __| | ___     | | ___  ___| |_ ___
% |  __/ _` / __| __| | |\/| |/ _ \ / _` |/ _ \    | |/ _ \/ __| __/ __|
% | | | (_| \__ \ |_  | |  | | (_) | (_| |  __/    | |  __/\__ \ |_\__ \
% |_|  \__,_|___/\__| |_|  |_|\___/ \__,_|\___|    |_|\___||___/\__|___/
%
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Fast%20Mode%20Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures start with 8

close all;
fprintf(1,'Figure: 8XXXXXX: Demo of fast mode cases\n');

%% Basic example - NO FIGURE
fig_num = 80001;
fprintf(1,'Figure: %.0f :Demo of fast mode, empty fig_num\n',fig_num);
figure(fig_num); close(fig_num);
plotting.FigureExpected = 0;

clear inputs
inputs.fig_num = [];
inputs.wall_start = [0 0; -3  3; 15 15; 7 7; 11 11;  7 3; 0 10; 10 16; 18 18; 20 0; 5.9 6; 13 14];
inputs.wall_end   = [5 0;  3 -3; 15 10; 9 9; 11 11; -1 7; 5 15; 14 16; 20 20; 20 0; 4 7;   13 10];
inputs.sensor_vector_start = [0 0];
inputs.sensor_vector_end   = [15 15];
inputs.flag_search_return_type = 1;
inputs.flag_search_range_type = 0;
inputs.tolerance = [];

% SHORT format checking
clear expected
expected.distance = [0         0   21.2132    9.8995   15.5563    6.1283   18.3848   12.7279]';
expected.location = [...
    0         0   15.0000    7.0000   11.0000    4.3333   13.0000    9.0000
    0         0   15.0000    7.0000   11.0000    4.3333   13.0000    9.0000]';
expected.path_segment = [1     2     3     4     5     6    12     4]';
expected.t = [0    0.5000         0         0         0    0.3333    0.2500    1.0000]';
expected.u = [0         0    1.0000    0.4667    0.7333    0.2889    0.8667    0.6000]';

actual = struct;
[actual.distance, actual.location, actual.path_segment, actual.t, actual.u] = ...
    fcn_Path_findSensorHitOnWall(...
    inputs.wall_start, inputs.wall_end,...
    inputs.sensor_vector_start,inputs.sensor_vector_end,...
    (inputs.flag_search_return_type), (inputs.flag_search_range_type), ...
    (inputs.tolerance), (inputs.fig_num));

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f :Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);
plotting.FigureExpected = 0;

clear inputs
inputs.fig_num = -1;
inputs.wall_start = [0 0; -3  3; 15 15; 7 7; 11 11;  7 3; 0 10; 10 16; 18 18; 20 0; 5.9 6; 13 14];
inputs.wall_end   = [5 0;  3 -3; 15 10; 9 9; 11 11; -1 7; 5 15; 14 16; 20 20; 20 0; 4 7;   13 10];
inputs.sensor_vector_start = [0 0];
inputs.sensor_vector_end   = [15 15];
inputs.flag_search_return_type = 1;
inputs.flag_search_range_type = 0;
inputs.tolerance = [];

% SHORT format checking
clear expected
expected.distance = [0         0   21.2132    9.8995   15.5563    6.1283   18.3848   12.7279]';
expected.location = [...
    0         0   15.0000    7.0000   11.0000    4.3333   13.0000    9.0000
    0         0   15.0000    7.0000   11.0000    4.3333   13.0000    9.0000]';
expected.path_segment = [1     2     3     4     5     6    12     4]';
expected.t = [0    0.5000         0         0         0    0.3333    0.2500    1.0000]';
expected.u = [0         0    1.0000    0.4667    0.7333    0.2889    0.8667    0.6000]';

actual = struct;
[actual.distance, actual.location, actual.path_segment, actual.t, actual.u] = ...
    fcn_Path_findSensorHitOnWall(...
    inputs.wall_start, inputs.wall_end,...
    inputs.sensor_vector_start,inputs.sensor_vector_end,...
    (inputs.flag_search_return_type), (inputs.flag_search_range_type), ...
    (inputs.tolerance), (inputs.fig_num));

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

inputs.wall_start = [0 0; -3  3; 15 15; 7 7; 11 11;  7 3; 0 10; 10 16; 18 18; 20 0; 5.9 6; 13 14];
inputs.wall_end   = [5 0;  3 -3; 15 10; 9 9; 11 11; -1 7; 5 15; 14 16; 20 20; 20 0; 4 7;   13 10];
inputs.sensor_vector_start = [0 0];
inputs.sensor_vector_end   = [15 15];
inputs.flag_search_return_type = 1;
inputs.flag_search_range_type = 0;
inputs.tolerance = [];

Niterations = 100;

% Do calculation without pre-calculation
inputs.fig_num = [];
tic;
for ith_test = 1:Niterations
    % Call the function
    actual = struct;
    [actual.distance, actual.location, actual.path_segment, actual.t, actual.u] = ...
        fcn_Path_findSensorHitOnWall(...
        inputs.wall_start, inputs.wall_end,...
        inputs.sensor_vector_start,inputs.sensor_vector_end,...
        (inputs.flag_search_return_type), (inputs.flag_search_range_type), ...
        (inputs.tolerance), (inputs.fig_num)); %#ok<STRNU>
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
inputs.fig_num = -1;
tic;
for ith_test = 1:Niterations
    % Call the function
    actual = struct;
    [actual.distance, actual.location, actual.path_segment, actual.t, actual.u] = ...
        fcn_Path_findSensorHitOnWall(...
        inputs.wall_start, inputs.wall_end,...
        inputs.sensor_vector_start,inputs.sensor_vector_end,...
        (inputs.flag_search_return_type), (inputs.flag_search_range_type), ...
        (inputs.tolerance), (inputs.fig_num)); %#ok<STRNU>

end
fast_method = toc;

% Plot results as bar chart
figure(373737);
clf;
X = categorical({'Normal mode','Fast mode'});
X = reordercats(X,{'Normal mode','Fast mode'}); % Forces bars to appear in this exact order, not alphabetized
Y = [slow_method fast_method ]*1000/Niterations;
bar(X,Y)
ylabel('Execution time (Milliseconds)')


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% BUG cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ____  _    _  _____
% |  _ \| |  | |/ ____|
% | |_) | |  | | |  __    ___ __ _ ___  ___  ___
% |  _ <| |  | | | |_ |  / __/ _` / __|/ _ \/ __|
% | |_) | |__| | |__| | | (_| (_| \__ \  __/\__ \
% |____/ \____/ \_____|  \___\__,_|___/\___||___/
%
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=BUG%20cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All bug case figures start with the number 9

close all;

%% BUG - infinite point intersection (4), flag (1), test 1 - should be infinite intersections

fig_num = 90001;
fprintf(1,'Figure: %.0f : known BUG (solved?)\n',fig_num);
figure(fig_num); clf;
plotting.FigureExpected = 1;

clear inputs
inputs.fig_num = fig_num;
inputs.wall_start = [0 10];
inputs.wall_end   = [10 10];
inputs.sensor_vector_start = [13 10];
inputs.sensor_vector_end   = [11 10];
inputs.flag_search_return_type = 0;
inputs.flag_search_range_type = 1;
inputs.tolerance = [];

% SHORT format checking
clear expected
expected.distance = 3;
expected.location = [10 10];
expected.path_segment = 1;
expected.t = 1;
expected.u = 1.5;

actual = struct;
[actual.distance, actual.location, actual.path_segment, actual.t, actual.u] = ...
    fcn_Path_findSensorHitOnWall(...
    inputs.wall_start, inputs.wall_end,...
    inputs.sensor_vector_start,inputs.sensor_vector_end,...
    (inputs.flag_search_return_type), (inputs.flag_search_range_type), ...
    (inputs.tolerance), (inputs.fig_num));

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)



%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§

% %% fcn_INTERNAL_printResults
% function fcn_INTERNAL_printResults(distance,location) 
% fprintf(1,'Distance \t Location X \t Location Y \n');
% if ~isempty(distance)
%     for i_result = 1:length(distance(:,1))
%         fprintf(1,'%.3f \t\t %.3f \t\t\t %.3f\n',distance(i_result),location(i_result,1),location(i_result,2));
%     end
% end
% end % Ends fcn_INTERNAL_printResults

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
    if ~any(isnan(thisVariable),'all') && ~isequal(round(thisVariable,expectedValueDigits),round(expectedValue,expectedValueDigits))
        flag_errorWillBeThrown = 1;
    end
    if ~any(isnan(expectedValue),'all') && ~isequal(round(thisVariable,expectedValueDigits),round(expectedValue,expectedValueDigits))
        flag_errorWillBeThrown = 1;
    end

    if all(isnan(thisVariable),'all') && ~all(isnan(thisVariable),'all')
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
        assert(isequal(round(thisVariable,expectedValueDigits),round(expectedValue,expectedValueDigits)));
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

% Tolerance is the fourth number
flag_testTolerance = fourthFigureNumber;

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

        % Nintersections = 1;

        % If ANY projections are used, the end points of the range will have
        % intersections. But if GIVEN projections are used, only the values from
        % 2:end-1 will have intersections. The trimEnds function gives values from
        % 2:end-1 for vectors that have length greater than 2.
        flag_search_range_type = thirdFigureNumber;
        switch flag_search_range_type
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
        totalRangeX = [-1.5 -0.4 0 0.4 1.5];
        totalRangeY = 0;

        % If ANY projections are used, the end points of the range will have
        % intersections. But if GIVEN projections are used, only the values from
        % 2:end-1 will have intersections. The trimEnds function gives values from
        % 2:end-1 for vectors that have length greater than 2.
        flag_search_range_type = thirdFigureNumber;
        switch flag_search_range_type
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

    case {4} % Multi-hit cases
        % All are tested with three walls. The walls are spaced so that,
        % typically, only 1 or 2 are hit. The third wall is spaced far away
        % so that sensor extensions cause it to be hit also.
        wall_starts          = [0 0; 0 0.5; -1 3];
        wall_ends            = [1 0; 1 0.5; 2 3];
        sensorLengthX = 0;
        sensorLengthY = 1;
        totalRangeX = [-1 0 0.5 1 2];
        totalRangeY = [-2 -1 -0.5 -0.25 0 1];

        % If ANY projections are used, the end points of the range will have
        % intersections. But if GIVEN projections are used, only the values from
        % 2:end-1 will have intersections. The trimEnds function gives values from
        % 2:end-1 for vectors that have length greater than 2.

        % intersections    = [nan nan];
        % distances     = nan;
        % wall_segments = nan;
        % tValues       = nan;
        % uValues       = nan;

        flag_search_range_type = thirdFigureNumber;
        switch flag_search_range_type
            case {0}
                % Limited sensor, limited wall
                clear intersection_counts
                intersection_counts = cell(2,3);

                IntersectionNumber = 1;
                N_intersections = 0;
                if flag_testTolerance ~= 2
                    x_values = [0 0.5 1];
                    y_values = -1;
                else
                    x_values =  0.5;
                    y_values = [-0.5 0];
                end
                for jth_y = 1:length(y_values)
                    thisY = y_values(jth_y);
                    for ith_x = 1:length(x_values)
                        thisX = x_values(ith_x);
                        N_intersections = N_intersections+1;
                        intersection_counts{IntersectionNumber,N_intersections} = [thisX   thisY];                        
                    end
                end

                IntersectionNumber = 2;
                N_intersections = 0;
                if flag_testTolerance ~= 2
                    x_values = [0 0.5 1];
                    y_values = [-0.5 -0.25 0];
                else
                    x_values = 0.5;
                    y_values = -0.25;
                end
                for jth_y = 1:length(y_values)
                    thisY = y_values(jth_y);
                    for ith_x = 1:length(x_values)
                        thisX = x_values(ith_x);
                        N_intersections = N_intersections+1;
                        intersection_counts{IntersectionNumber,N_intersections} = [thisX   thisY];
                    end
                end

            case {1}
                % Infinite sensor, limited wall
                clear intersection_counts
                intersection_counts = cell(3,4);

                % 1 intersection case
                IntersectionNumber = 1;
                N_intersections = 0;
                if flag_testTolerance ~= 2
                    x_values = [-1  2];
                    y_values = [-2 -1 -0.5 -0.25 0 1];                    
                else
                    x_values = [0 1];
                    y_values = [-2 -1 -0.5 -0.25 0 1];                    
                end
                for jth_y = 1:length(y_values)
                    thisY = y_values(jth_y);
                    for ith_x = 1:length(x_values)
                        thisX = x_values(ith_x);
                        N_intersections = N_intersections+1;
                        intersection_counts{IntersectionNumber,N_intersections} = [thisX   thisY];                        
                    end
                end

                % 2 intersection cases
                % NONE!

                % 3 intersection cases
                IntersectionNumber = 3;
                N_intersections = 0;
                if flag_testTolerance ~= 2
                    x_values = [0 0.5 1];
                    y_values = [-2 -1 -0.5 -0.25 0 1];
                else
                    x_values = 0.5;
                    y_values = [-2 -1 -0.5 -0.25 0 1];
                end
                for jth_y = 1:length(y_values)
                    thisY = y_values(jth_y);
                    for ith_x = 1:length(x_values)
                        thisX = x_values(ith_x);
                        N_intersections = N_intersections+1;
                        intersection_counts{IntersectionNumber,N_intersections} = [thisX   thisY];
                    end
                end
            case {2}
                % Limited sensor, infinite wall

                % 



                clear intersection_counts
                intersection_counts = cell(3,4);

                % 1 intersection cases
                % wall_starts          = [0 0; 0 0.5; -1 3];
                % wall_ends            = [1 0; 1 0.5; 2 3];
                % sensorLengthX = 0;
                % sensorLengthY = 1;
                % totalRangeX = [-1 0 0.5 1 2];
                % totalRangeY = [-2 -1 -0.5 -0.25 0 1];

                IntersectionNumber = 1;
                N_intersections = 0;
                
                x_values = [-1 0 0.5 1 2];
                if flag_testTolerance ~= 2
                    y_values = -1;
                else
                    y_values = [-0.25 -0.5 0];
                end
                for jth_y = 1:length(y_values)
                    thisY = y_values(jth_y);
                    for ith_x = 1:length(x_values)
                        thisX = x_values(ith_x);
                        N_intersections = N_intersections+1;
                        intersection_counts{IntersectionNumber,N_intersections} = [thisX   thisY];                        
                    end
                end

                % 2 intersection cases
                IntersectionNumber = 2;
                N_intersections = 0;
                x_values = [-1 0 0.5 1 2];
                if flag_testTolerance ~= 2
                    y_values = [-0.5 -0.25 0];
                else
                    y_values = -0.25;
                end

                for jth_y = 1:length(y_values)
                    thisY = y_values(jth_y);
                    for ith_x = 1:length(x_values)
                        thisX = x_values(ith_x);
                        N_intersections = N_intersections+1;
                        intersection_counts{IntersectionNumber,N_intersections} = [thisX   thisY];
                    end
                end

                % 3 intersection cases
                % NONE!
            case {3}
                % Infinite sensor, infinite wall

                clear intersection_counts
                intersection_counts = cell(3,4);

                % 1 intersection cases
                % NONE!

                % 2 intersection cases
                % NONE!

                % 3 intersection cases
                % wall_starts          = [0 0; 0 0.5; -1 3];
                % wall_ends            = [1 0; 1 0.5; 2 3];
                % sensorLengthX = 0;
                % sensorLengthY = 1;
                % totalRangeX = [-1 0 0.5 1 2];
                % totalRangeY = [-2 -1 -0.5 -0.25 0 1];
                IntersectionNumber = 3;
                N_intersections = 0;
                x_values = [-1 0 0.5 1 2];
                y_values = [-2 -1 -0.5 -0.25 0 1];

                for jth_y = 1:length(y_values)
                    thisY = y_values(jth_y);
                    for ith_x = 1:length(x_values)
                        thisX = x_values(ith_x);
                        N_intersections = N_intersections+1;
                        intersection_counts{IntersectionNumber,N_intersections} = [thisX   thisY];
                    end
                end
            otherwise
                warning('on','backtrace');
                warning('Expecting a third figure number with integer values of 0 to 3, but found: %.3f',thirdFigureNumber);
                error('Third figure number not recognized');
        end

    case {5} % Wall as a single point case

        % All are tested one "single point" wall. 
        wall_starts          = [0 0];
        wall_ends            = [0 0];
        sensorLengthX = 0;
        sensorLengthY = 1;
        totalRangeX = [-1 0 1];
        totalRangeY = [-2 -1 -0.5 0 1];

        flag_search_range_type = thirdFigureNumber;
        switch flag_search_range_type
            case {0}
                % Limited sensor, limited wall
                clear intersection_counts
                intersection_counts = cell(2,3);

                IntersectionNumber = 1;
                N_intersections = 0;
                if flag_testTolerance ~= 2
                    x_values = 0;
                    y_values =  [-1 -0.5 0 ];
                else
                    x_values = [];
                    y_values =  [];
                end
                for jth_y = 1:length(y_values)
                    thisY = y_values(jth_y);
                    for ith_x = 1:length(x_values)
                        thisX = x_values(ith_x);
                        N_intersections = N_intersections+1;
                        intersection_counts{IntersectionNumber,N_intersections} = [thisX   thisY];                        
                    end
                end


            case {1}
                % Infinite sensor, limited wall
                clear intersection_counts
                intersection_counts = cell(3,4);

                % 1 intersection case
                IntersectionNumber = 1;
                N_intersections = 0;
                if flag_testTolerance ~= 2
                    x_values = 0;
                    y_values = [-2 -1 -0.5 -0.25 0 1];                    
                else
                    x_values = [];
                    y_values = [];
                end
                for jth_y = 1:length(y_values)
                    thisY = y_values(jth_y);
                    for ith_x = 1:length(x_values)
                        thisX = x_values(ith_x);
                        N_intersections = N_intersections+1;
                        intersection_counts{IntersectionNumber,N_intersections} = [thisX   thisY];                        
                    end
                end

               
            case {2}
                % Limited sensor, infinite wall

                clear intersection_counts
                intersection_counts = cell(3,4);

                % 1 intersection cases
                IntersectionNumber = 1;
                N_intersections = 0;
                
                x_values = [-1 0 1];
                if flag_testTolerance ~= 2
                    y_values =  [-1 -0.5 0 ];
                else
                    y_values = -0.5;
                end
                for jth_y = 1:length(y_values)
                    thisY = y_values(jth_y);
                    for ith_x = 1:length(x_values)
                        thisX = x_values(ith_x);
                        N_intersections = N_intersections+1;
                        intersection_counts{IntersectionNumber,N_intersections} = [thisX   thisY];                        
                    end
                end

            case {3}
                % Infinite sensor, infinite wall

                clear intersection_counts
                intersection_counts = cell(3,4);

                % 1 intersection cases
                IntersectionNumber = 1;
                N_intersections = 0;
                x_values = [-1 0 1];
                y_values = [-2 -1 -0.5 -0.25 0 1];

                for jth_y = 1:length(y_values)
                    thisY = y_values(jth_y);
                    for ith_x = 1:length(x_values)
                        thisX = x_values(ith_x);
                        N_intersections = N_intersections+1;
                        intersection_counts{IntersectionNumber,N_intersections} = [thisX   thisY];
                    end
                end
            otherwise
                warning('on','backtrace');
                warning('Expecting a third figure number with integer values of 0 to 3, but found: %.3f',thirdFigureNumber);
                error('Third figure number not recognized');
        end

    case {6} % Sensor as a single point case

        % All are tested one "single point" sensor.
        wall_starts          = [0 0;  0 -1; -1 0;  0 2; 0 -3; 1 1];
        wall_ends            = [0 1;  0  1; -1 1;  1 2; 0 -2; 1 2];
        sensorLengthX = 0;
        sensorLengthY = 0;
        totalRangeX = 0;
        totalRangeY = 0;

        flag_search_range_type = thirdFigureNumber;
        switch flag_search_range_type
            case {0}
                % Limited sensor, limited wall 
                % 0: GIVEN point of sensor                     is inside GIVEN projection of walls (default)
                clear intersection_counts
                intersection_counts = cell(2,3);

                % wall_starts          = [0 0;  0 -1; -1 0;  0 2; 0 -3; 1 1];
                % wall_ends            = [0 1;  0  1; -1 1;  1 2; 0 -2; 1 2];

                if flag_testTolerance ~= 2
                    IntersectionNumber = 2;
                    intersection_counts{IntersectionNumber,1} = [0 0];
                else
                    % No intersections
                end
                        

            case {1}
                % Infinite sensor, limited wall
                % 1: ORTHO PROJECTION of sensor onto each wall is inside GIVEN projection of walls 
                clear intersection_counts
                intersection_counts = cell(2,3);

                % wall_starts          = [0 0;  0 -1; -1 0;  0 2; 0 -3; 1 1];
                % wall_ends            = [0 1;  0  1; -1 1;  1 2; 0 -2; 1 2];

                if flag_testTolerance ~= 2
                    IntersectionNumber = 4;
                    intersection_counts{IntersectionNumber,1} = [0 0];
                else
                    IntersectionNumber = 1;
                    intersection_counts{IntersectionNumber,1} = [0 0];
                end

            case {2}
                % Limited sensor, infinite wall
                % 2: GIVEN point of sensor                     is inside ANY   projection of walls
                % wall_starts          = [0 0;  0 -1; -1 0;  0 2; 0 -3; 1 1];
                % wall_ends            = [0 1;  0  1; -1 1;  1 2; 0 -2; 1 2];
                clear intersection_counts
                intersection_counts = cell(2,3);

                if flag_testTolerance ~= 2
                    IntersectionNumber = 3;
                    intersection_counts{IntersectionNumber,1} = [0 0];
                else
                    % No intersections
                end

            case {3}
                % Infinite sensor, infinite wall
                % 3: ORTHO PROJECTION of sensor onto each wall is inside ANY   projection of walls
                % wall_starts          = [0 0;  0 -1; -1 0;  0 2; 0 -3; 1 1];
                % wall_ends            = [0 1;  0  1; -1 1;  1 2; 0 -2; 1 2];
                clear intersection_counts
                intersection_counts = cell(2,3);

                IntersectionNumber = 6;
                intersection_counts{IntersectionNumber,1} = [0 0];
                
            otherwise
                warning('on','backtrace');
                warning('Expecting a third figure number with integer values of 0 to 3, but found: %.3f',thirdFigureNumber);
                error('Third figure number not recognized');
        end

    case {7} % All multi-hit cases
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
switch flag_testTolerance
    case {0}
        % do nothing
        tolerance = [];
    case {1}
        % Positive tolerance
        tolerance = 0.001;
    case {2}
        % Negative tolerance
        tolerance = -0.001;
        if firstFigureNumber<4
            % xStarts are trimmed?
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
        % if ith_case == 21
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
                    intersections    = [thisX 0];
                    distances    = 0 - thisY;
                    wall_segments    = 1;
                    testCases(ith_case).expected.wall_segment = 1;
                    tValues     = thisX;
                    uValues     = -thisY;
                else
                    intersections    = [nan nan];
                    distances    = nan;
                    wall_segments = nan;
                    testCases(ith_case).expected.wall_segment = nan;
                    tValues     = nan;
                    uValues     = nan;
                end

            case {3} % Infinite intersect cases

                wall_segments    = [1; 1]; % Default
                if thisX<0
                    distanceSensorStartToWallStart = 0 - thisX;
                elseif thisX>=0 && thisX<=1
                    distanceSensorStartToWallStart = 0;
                else
                    distanceSensorStartToWallStart = 1 - thisX;                    
                end
                maxXstart = max(thisX,0);
                minXend   = min(thisX+1,1);


                % If tolerance is negative, the lines NEVER intersect
                % (fourthFigureNumber == 2)
                if ismember(thisX, xStartsWithIntersections) && fourthFigureNumber~=2
                    switch flag_search_range_type
                        case {0}

                            % Finite sensor, finite walls
                            % if thirdFigureNumber==0 || secondFigureNumber==0
                            %     intersections    = [maxXstart 0; minXend 0];
                            % else
                            %     intersections    = [maxXstart 0; minXend 0];
                            % end
                            intersections    = [maxXstart 0; minXend 0];
                            distances   = sum(([thisX 0; thisX 0] - intersections).^2,2).^0.5;
                            tValues     = [maxXstart; minXend];
                            uValues     = [distanceSensorStartToWallStart; minXend-thisX];
                        case {1}
                            % Infinite sensor, finite walls
                            intersections = [0 0; 1 0];
                            distances     = [0-thisX; 1-thisX];
                            tValues       = [0; 1];
                            uValues       = [0-thisX; 0-thisX+1];
                            if thisX>1 && flag_search_return_type==0
                                % Vector is beyond end of segment, so
                                % closest intersection is the endpoint
                                intersections = flipud(intersections);
                                distances     = flipud(distances);
                                tValues       = flipud(tValues);
                                uValues       = flipud(uValues);
                            end
                        case {2}
                            % Finite sensor, infinite walls
                            intersections = [thisX 0; thisX+1 0];
                            distances     = [0; 1];
                            tValues       = [thisX; thisX+1];
                            uValues       = [0; 1];
                        case {3}
                            % Infinite sensor, infinite walls
                            intersections = [-inf 0; inf 0];
                            distances     = [-inf; inf];
                            tValues       = [-inf; inf];
                            uValues       = [-inf; inf];
                        otherwise

                    end

                else
                    intersections    = [nan nan];
                    distances     = nan;
                    wall_segments = nan;
                    tValues       = nan;
                    uValues       = nan;
                end

            case {4} % Multi-hit cases

                % Find how many intersections there are, and where                
                % intersections are at
                Num_intersections = 0; % Default value
                intersectionCountSize = size(intersection_counts);
                for nth_interection = 1:intersectionCountSize(1)
                    for jth_case = 1:intersectionCountSize(2)
                        this_intersection = intersection_counts{nth_interection,jth_case};
                        if isequal([thisX thisY], this_intersection)
                            Num_intersections = nth_interection;
                        end
                    end
                end
                switch flag_search_range_type
                    case {0, 2, 3}
                        if Num_intersections==1
                            intersections = [thisX 0];
                            distances     = -thisY;
                            wall_segments = 1;
                            tValues       = thisX;
                            uValues       = distances;

                            % This is a very special case where tolerance
                            % changes hits that count as either 1st or 2nd
                            if 2==flag_testTolerance
                                intersections = [thisX thisY+0.5];
                                distances     = 0.5;
                                if thisY==-0.5
                                    wall_segments = 1;
                                else
                                    wall_segments = 2;
                                end
                                tValues       = thisX;
                                uValues       = distances;

                            end
                        elseif Num_intersections==2
                            intersections = [thisX 0; thisX 0.5];
                            distances     = [0-thisY; 0.5-thisY];
                            wall_segments = [1; 2];
                            tValues       = [thisX; thisX];
                            uValues       = distances;
                        elseif Num_intersections==3
                            intersections = [thisX 0; thisX 0.5; thisX 3];
                            distances     = [0-thisY; 0.5-thisY; 3-thisY];
                            wall_segments = [1; 2; 3];
                            tValues       = [thisX; thisX; (thisX+1)/3];
                            uValues       = distances;
                        else
                            intersections    = [nan nan];
                            distances     = nan;
                            wall_segments = nan;
                            tValues       = nan;
                            uValues       = nan;
                        end
                    case{1}
                        if Num_intersections==1
                            intersections = [thisX 3];
                            distances     = 3-thisY;
                            wall_segments = 3;
                            tValues       = (thisX+1)/3;
                            uValues       = distances;
                        elseif Num_intersections==2
                            error('This does not happen?');
                        elseif Num_intersections==3                            
                            intersections = [thisX 0; thisX 0.5; thisX 3];
                            distances     = [0-thisY; 0.5-thisY; 3-thisY];
                            wall_segments = [1; 2; 3];
                            tValues       = [thisX; thisX; (thisX+1)/3];
                            uValues       = distances;
                        else
                            intersections    = [nan nan];
                            distances     = nan;
                            wall_segments = nan;
                            tValues       = nan;
                            uValues       = nan;
                        end
                    otherwise
                end
            case {5} % Single point wall case

                % Find how many intersections there are, and where                
                % intersections are at
                Num_intersections = 0; % Default value
                intersectionCountSize = size(intersection_counts);
                for nth_interection = 1:intersectionCountSize(1)
                    for jth_case = 1:intersectionCountSize(2)
                        this_intersection = intersection_counts{nth_interection,jth_case};
                        if isequal([thisX thisY], this_intersection)
                            Num_intersections = nth_interection;
                        end
                    end
                end

                switch flag_search_range_type
                    case {0, 1, 2, 3}
                        if Num_intersections==1
                            intersections = [thisX 0];
                            distances     = 0-thisY;
                            wall_segments = 1;
                            tValues       = thisX;
                            uValues       = distances;
                        else
                            intersections    = [nan nan];
                            distances     = nan;
                            wall_segments = nan;
                            tValues       = nan;
                            uValues       = nan;
                        end
                    case{7}
                        % Infinite sensor, finite wall
                        if Num_intersections==1
                            intersections = [thisX 0-x];
                            distances     = 3-thisY;
                            wall_segments = 3;
                            tValues       = (thisX+1)/3;
                            uValues       = distances;
                        else
                            intersections    = [nan nan];
                            distances     = nan;
                            wall_segments = nan;
                            tValues       = nan;
                            uValues       = nan;
                        end
                    otherwise
                end

            case {6} % Single point sensor case

                % Find how many intersections there are, and where
                % intersections are at
                Num_intersections = 0; % Default value
                intersectionCountSize = size(intersection_counts);
                for nth_interection = 1:intersectionCountSize(1)
                    for jth_case = 1:intersectionCountSize(2)
                        this_intersection = intersection_counts{nth_interection,jth_case};
                        if isequal([thisX thisY], this_intersection)
                            Num_intersections = nth_interection;
                        end
                    end
                end

                switch Num_intersections
                    case {0}
                        intersections    = [nan nan];
                        distances     = nan;
                        wall_segments = nan;
                        tValues       = nan;
                        uValues       = nan;                        
                    case {1}
                        % Only occurs in the flag_search_range_type=1, flag_testTolerance=2 case
                        assert(1==flag_search_range_type);
                        assert(2==flag_testTolerance);
                        intersections = [0 0];
                        distances     = 0;
                        wall_segments = 2;
                        tValues       = 0.5;
                        uValues       = 0;
                    case{2}
                        % Only occurs in the flag_search_range_type=0 case
                        assert(0==flag_search_range_type);
                        intersections = [0 0; 0 0];
                        distances     = [0; 0];
                        wall_segments = [1; 2];
                        tValues       = [0; 0.5];
                        uValues       = [0; 0];
                    case{3}
                        % Only occurs in the flag_search_range_type=2 case
                        assert(2==flag_search_range_type);
                        intersections = [0 0; 0 0; 0 0];
                        distances     = [0; 0; 0];
                        wall_segments = [1; 2; 5];
                        tValues       = [0; 0.5; 3];
                        uValues       = [0; 0; 0];                        
                    case{4}
                        % Only occurs in the flag_search_range_type=1 case
                        assert(1==flag_search_range_type);
                        intersections = [0 0; 0 0; -1 0; 0 2];
                        distances     = [0; 0; -1; -2];
                        wall_segments = [1; 2; 3; 4];
                        tValues       = [0; 0.5; 0; 0];
                        uValues       = [0; 0; -1; -2];                        
                    case{6}
                        % Only occurs in the flag_search_range_type=3 case
                        assert(3==flag_search_range_type);
                        intersections = [0 0; 0 0; -1 0; 0 2; 0 0; 1 0];
                        distances     = [0; 0; -1; -2; 0; 1];
                        wall_segments = [1; 2; 3; 4; 5; 6];
                        tValues       = [0; 0.5; 0; 0; 3; -1];
                        uValues       = [0; 0; -1; -2; 0; 1];
                    otherwise
                end
            otherwise
                warning('on','backtrace');
                warning('Expecting a first figure number with integer values of 2 to 7, but found: %.3f',firstFigureNumber);
                error('Third figure number not recognized');
        end
        % Check if 2 solutions or 1 solution requested
        if flag_search_return_type==0
            if all(isnan(intersections))
                closestIndex = 1;
            else
                closestIndex = fcn_INTERNAL_selectClosestPoint([thisX thisY], intersections, flag_search_return_type);
            end
            testCases(ith_case).expected.Intersection    = intersections(closestIndex,:);
            testCases(ith_case).expected.Distance        = distances(closestIndex);
            testCases(ith_case).expected.wall_segment    = wall_segments(closestIndex);
            testCases(ith_case).expected.t               = tValues(closestIndex);
            testCases(ith_case).expected.u               = uValues(closestIndex);
        else
            testCases(ith_case).expected.Intersection    = intersections;
            testCases(ith_case).expected.Distance        = distances;
            testCases(ith_case).expected.wall_segment    = wall_segments;
            testCases(ith_case).expected.t               = tValues;
            testCases(ith_case).expected.u               = uValues;
        end
    end
end

fcn_INTERNAL_plotTestCases(fig_num, wall_starts, wall_ends, testCases);

end % Ends fcn_INTERNAL_fillTestCases

%% fcn_INTERNAL_selectClosestPoint 
function [within_indices, distances_squared] = fcn_INTERNAL_selectClosestPoint(sensor_vector_start, intersections, flag_search_return_type)

% Find the distances via Euclidian distance to the sensor's origin
% note: a faster way to do this might be to just
% calculate t*r as a length
distances_squared = sum((intersections - sensor_vector_start).^2,2);
within_indices = find(~isnan(distances_squared));
if ~isempty(within_indices)
    if 0==flag_search_return_type
        % Keep only the minimum distance result
        [~,within_indices] = min(distances_squared);
 
    elseif 1==flag_search_return_type
        % Return all the results by default
    else
        warning('on','backtrace');
        warning('Expecting a flag_search_return_type as integer with values of 0 or 1, but found: %.3f',flag_search_return_type);
        error('Bad flag_search_range_type encountered');
    end
end

end % Ends fcn_INTERNAL_selectClosestPoint