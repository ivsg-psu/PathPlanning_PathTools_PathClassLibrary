% script_test_fcn_Path_fcn_Path_convertPerA2PerB
% This is a script to exercise the function: fcn_Path_fcn_Path_convertPerA2PerB.m
% This function was written on 2025_06_19 by S. Brennan, from original
% Questions or comments? sbrennan@psu.edu

% Modification history:
% 2025_06_19 - S. Brennan
% -- first code write

close all
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
inputs.fig_num   = fig_num;
inputs.start_A   = [0 0];
inputs.start_B   = [0 10];
inputs.vector_A   = [0 5];
inputs.vector_B   = [0 10];
inputs.A_values   = [0; 1];

% SHORT format checking
clear expected
expected.percentage_B = [-1; -0.5];

actual = struct;
actual.percentage_B = ...
    fcn_Path_convertPerA2PerB(...
    inputs.start_A, inputs.start_B, ...
    inputs.vector_A, inputs.vector_B, inputs.A_values,...
    (inputs.fig_num));

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)
%fcn_INTERNAL_printResults(actual.distance,actual.location);

%% Demo of long format checking
fig_num = 10002;
fprintf(1,'Figure: %.0f :Demo of long format checking\n',fig_num);
figure(fig_num); clf;
plotting.FigureExpected = 1;

clear inputs
inputs.fig_num   = fig_num;
inputs.start_A   = [0 0];
inputs.start_B   = [0 10];
inputs.vector_A   = [0 5];
inputs.vector_B   = [0 10];
inputs.A_values   = [0; 1];

% LONG format checking
clear expected
expected.percentage_B.TypeString      = 'numeric'; % See 'help isa' for a listing
expected.percentage_B.Size            = [length(inputs.A_values) 1];
expected.percentage_B.Value           = [-1; -0.5];
expected.percentage_B.ValueDigits     = 4; % How many digits must match

actual = struct;
actual.percentage_B = ...
    fcn_Path_convertPerA2PerB(...
    inputs.start_A, inputs.start_B, ...
    inputs.vector_A, inputs.vector_B, inputs.A_values,...
    (inputs.fig_num));

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)
%fcn_INTERNAL_printResults(actual.distance,actual.location);

%% Demo of negative values
fig_num = 10003;
fprintf(1,'Figure: %.0f : A and B as single vectors, vector A_values \n',fig_num);
figure(fig_num); clf;
plotting.FigureExpected = 1;

clear inputs
inputs.fig_num   = fig_num;
inputs.start_A   = [-2 0];
inputs.start_B   = [4 0];
inputs.vector_A   = [4 0];
inputs.vector_B   = [10 0];
inputs.A_values   = [0; 1];

% SHORT format checking
clear expected
expected.percentage_B = [-0.6; -0.2];

actual = struct;
actual.percentage_B = ...
    fcn_Path_convertPerA2PerB(...
    inputs.start_A, inputs.start_B, ...
    inputs.vector_A, inputs.vector_B, inputs.A_values,...
    (inputs.fig_num));

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)


%% Demo of parallel lines
fig_num = 10004;
fprintf(1,'Figure: %.0f :Demo of parallel lines \n',fig_num);
figure(fig_num); clf;
plotting.FigureExpected = 1;

clear inputs
inputs.fig_num   = fig_num;
inputs.start_A   = [0 0];
inputs.start_B   = [0 5];
inputs.vector_A   = [5 5];
inputs.vector_B   = [10 10];
inputs.A_values   = [0; 1];

% SHORT format checking
clear expected
expected.percentage_B = [-0.25; 0.25];

actual = struct;
actual.percentage_B = ...
    fcn_Path_convertPerA2PerB(...
    inputs.start_A, inputs.start_B, ...
    inputs.vector_A, inputs.vector_B, inputs.A_values,...
    (inputs.fig_num));

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)


%% Demo of vectorized B inputs, A as a single vector
fig_num = 10005;
fprintf(1,'Figure: %.0f :Demo of vectorized B inputs, A as single vector\n',fig_num);
figure(fig_num); clf;
plotting.FigureExpected = 1;

clear inputs
inputs.fig_num   = fig_num;
inputs.start_A   = [0 0];
inputs.start_B   = [0 5; -5 10; 0 15];
inputs.vector_A   = [5 0];
inputs.vector_B   = [10 0; 5 0; 20 0];
inputs.A_values   = [0; 1; 0.5];

% SHORT format checking
clear expected
expected.percentage_B = [0; 2; 0.125];

actual = struct;
actual.percentage_B = ...
    fcn_Path_convertPerA2PerB(...
    inputs.start_A, inputs.start_B, ...
    inputs.vector_A, inputs.vector_B, inputs.A_values,...
    (inputs.fig_num));

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)

%% Demo of vectorized A and B input pairs
fig_num = 10006;
fprintf(1,'Figure: %.0f :Demo of A and B vectors paired together\n',fig_num);
figure(fig_num); clf;
plotting.FigureExpected = 1;

clear inputs
inputs.fig_num   = fig_num;
inputs.start_A   = [0 0; 0 5;  -5 10; 5 15];
inputs.start_B   = [0 1; -5 6; 0 11; -5 16];
inputs.vector_A   = [5 0; 10 0; 5 0; 20 0];
inputs.vector_B   = [10 0; 5 0; 20 0; 5 0];
inputs.A_values   = [0; 1; 0.5; 0.5];

% SHORT format checking
clear expected
expected.percentage_B = [0; 3; -0.125; 4];

actual = struct;
actual.percentage_B = ...
    fcn_Path_convertPerA2PerB(...
    inputs.start_A, inputs.start_B, ...
    inputs.vector_A, inputs.vector_B, inputs.A_values,...
    (inputs.fig_num));

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)


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

%% Basic example - NO FIGURE

fig_num = 80001;
figure(fig_num);
close(fig_num);
plotting.FigureExpected = 0;

clear inputs
inputs.fig_num   = [];
inputs.start_A   = [0 0; 0 5;  -5 10; 5 15];
inputs.start_B   = [0 1; -5 6; 0 11; -5 16];
inputs.vector_A   = [5 0; 10 0; 5 0; 20 0];
inputs.vector_B   = [10 0; 5 0; 20 0; 5 0];
inputs.A_values   = [0; 1; 0.5; 0.5];

% SHORT format checking
clear expected
expected.percentage_B = [0; 3; -0.125; 4];

actual = struct;
actual.percentage_B = ...
    fcn_Path_convertPerA2PerB(...
    inputs.start_A, inputs.start_B, ...
    inputs.vector_A, inputs.vector_B, inputs.A_values,...
    (inputs.fig_num));

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)

%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
figure(fig_num);
close(fig_num);
plotting.FigureExpected = 0;

clear inputs
inputs.fig_num   = -1;
inputs.start_A   = [0 0; 0 5;  -5 10; 5 15];
inputs.start_B   = [0 1; -5 6; 0 11; -5 16];
inputs.vector_A   = [5 0; 10 0; 5 0; 20 0];
inputs.vector_B   = [10 0; 5 0; 20 0; 5 0];
inputs.A_values   = [0; 1; 0.5; 0.5];

% SHORT format checking
clear expected
expected.percentage_B = [0; 3; -0.125; 4];

actual = struct;
actual.percentage_B = ...
    fcn_Path_convertPerA2PerB(...
    inputs.start_A, inputs.start_B, ...
    inputs.vector_A, inputs.vector_B, inputs.A_values,...
    (inputs.fig_num));

fcn_INTERNAL_checkTestCases(inputs, expected, actual, plotting)

%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
figure(fig_num);
close(fig_num);

inputs.start_A   = [0 0; 0 5;  -5 10; 5 15];
inputs.start_B   = [0 1; -5 6; 0 11; -5 16];
inputs.vector_A   = [5 0; 10 0; 5 0; 20 0];
inputs.vector_B   = [10 0; 5 0; 20 0; 5 0];
inputs.A_values   = [0; 1; 0.5; 0.5];

Niterations = 1000;

% Do calculation without pre-calculation
inputs.fig_num   = [];
tic;
for ith_test = 1:Niterations
    % Call the function
    actual.percentage_B = ...
        fcn_Path_convertPerA2PerB(...
        inputs.start_A, inputs.start_B, ...
        inputs.vector_A, inputs.vector_B, inputs.A_values,...
        (inputs.fig_num));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
inputs.fig_num   = -1;
tic;
for ith_test = 1:Niterations
    % Call the function
    actual.percentage_B = ...
        fcn_Path_convertPerA2PerB(...
        inputs.start_A, inputs.start_B, ...
        inputs.vector_A, inputs.vector_B, inputs.A_values,...
        (inputs.fig_num));

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
        assert(isequal(round(thisVariable,expectedValueDigits),expectedValue));
    end

end

% Make sure plot opened up
if isfield(plotting,'FigureExpected') && plotting.FigureExpected==1
    assert(isequal(get(gcf,'Number'),inputs.fig_num));
end

end