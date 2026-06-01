% script_test_fcn_Path_calcDiffAnglesBetweenPathSegments.m
% Tests fcn_Path_calcDiffAnglesBetweenPathSegments
       
% Revision history:
% 2021_01_03
% - first write of the code
% 2021_01_07
% - fixed typos in the comments, minor header clean-ups
% 2025_08_02 by S. Brennan
% - In script_test_fcn_Path_calcDiffAnglesBetweenPathSegments
%   % * added edgeLengths to variable testing
%   % * renamed other variables to match inputs
%   % * added a trival but problem test case from real-world data

close all

%% Calculation of the angle between path segments
figNum = 10001;
titleString = sprintf('Calculation of the angle between path segments');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;


% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;
pathVerticesXY = paths_array{1}; % Pick first path as reference_traversal structure

[changeInAngles, edgeLengths] = fcn_Path_calcDiffAnglesBetweenPathSegments(pathVerticesXY,figNum);

title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(changeInAngles));
assert(isnumeric(edgeLengths));

% Check variable sizes
assert(isequal(size(changeInAngles),[length(pathVerticesXY(:,1))-2 1]));
assert(isequal(size(edgeLengths),[length(pathVerticesXY(:,1))-1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

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
figNum = 80001;
fprintf(1,'Figure: %.0f: Demo of fast mode, empty figNum\n',figNum);
figure(figNum); close(figNum);

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;
pathVerticesXY = paths_array{1}; % Pick first path as reference_traversal structure

[changeInAngles, edgeLengths] = fcn_Path_calcDiffAnglesBetweenPathSegments(pathVerticesXY,[]);

% Check variable types
assert(isnumeric(changeInAngles));
assert(isnumeric(edgeLengths));

% Check variable sizes
assert(isequal(size(changeInAngles),[length(pathVerticesXY(:,1))-2 1]));
assert(isequal(size(edgeLengths),[length(pathVerticesXY(:,1))-1 1]));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Basic fast mode - NO FIGURE, FAST MODE
figNum = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, figNum=-1\n',figNum);
figure(figNum); close(figNum);

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;
pathVerticesXY = paths_array{1}; % Pick first path as reference_traversal structure

[changeInAngles, edgeLengths] = fcn_Path_calcDiffAnglesBetweenPathSegments(pathVerticesXY,-1);

% Check variable types
assert(isnumeric(changeInAngles));
assert(isnumeric(edgeLengths));

% Check variable sizes
assert(isequal(size(changeInAngles),[length(pathVerticesXY(:,1))-2 1]));
assert(isequal(size(edgeLengths),[length(pathVerticesXY(:,1))-1 1]));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
figNum = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',figNum);
figure(figNum);
close(figNum);

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;
pathVerticesXY = paths_array{1}; % Pick first path as reference_traversal structure

Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [changeInAngles, edgeLengths] = fcn_Path_calcDiffAnglesBetweenPathSegments(pathVerticesXY,[]);
end
slow_method = toc;

% Check variable types
assert(isnumeric(changeInAngles));
assert(isnumeric(edgeLengths));

% Check variable sizes
assert(isequal(size(changeInAngles),[length(pathVerticesXY(:,1))-2 1]));
assert(isequal(size(edgeLengths),[length(pathVerticesXY(:,1))-1 1]));

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [changeInAngles, edgeLengths] = fcn_Path_calcDiffAnglesBetweenPathSegments(pathVerticesXY,-1);
end
fast_method = toc;

% Check variable types
assert(isnumeric(changeInAngles));
assert(isnumeric(edgeLengths));

% Check variable sizes
assert(isequal(size(changeInAngles),[length(pathVerticesXY(:,1))-2 1]));
assert(isequal(size(edgeLengths),[length(pathVerticesXY(:,1))-1 1]));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));

% Plot results as bar chart
figure(373737);
clf;
hold on;

X = categorical({'Normal mode','Fast mode'});
X = reordercats(X,{'Normal mode','Fast mode'}); % Forces bars to appear in this exact order, not alphabetized
Y = [slow_method fast_method ]*1000/Niterations;
bar(X,Y)
ylabel('Execution time (Milliseconds)')

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


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

% close all;

%% BUG Point that goes back on itself
figNum = 90001;
titleString = sprintf('Calculation of the angle between path segments');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;


pathVerticesXY = [0 0; 1 0; 0.5 0];

[changeInAngles, edgeLengths] = fcn_Path_calcDiffAnglesBetweenPathSegments(pathVerticesXY,figNum);

title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(changeInAngles));
assert(isnumeric(edgeLengths));

% Check variable sizes
assert(isequal(size(changeInAngles),[length(pathVerticesXY(:,1))-2 1]));
assert(isequal(size(edgeLengths),[length(pathVerticesXY(:,1))-1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));


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
