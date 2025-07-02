% script_test_fcn_Path_findPathWithMostData.m
% tests function: fcn_Path_findPathWithMostData.m

% Revision history
% 2020_11_10
% -- first write of the code
% 2021_01_07
% -- renamed function to change paths to traversals
% 2020_01_09
% -- added more comments during clean-up
% 2025_07_01 - S. Brennan
% -- Removed traversal input type and replaced with cell array of paths
% -- Renamed function from fcn_Path_findTraversalWithMostData

%% BASIC Example: finding longest traversal
fig_num = 10001;
titleString = sprintf('BASIC Example: finding longest path');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Call the function that fills in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;

% Convert paths to path structures
clear cellArrayOfPaths; cellArrayOfPaths = cell(3,1);

for i_Path = 1:length(paths)
    
    cellArrayOfPaths{i_Path,1} = paths{i_Path};
end

% Call the function
index_of_longest = fcn_Path_findPathWithMostData(cellArrayOfPaths, (fig_num));

% Update title
title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(index_of_longest));

% Check variable sizes
assert(isequal(size(index_of_longest),[1 1]));

% Check variable values
assert(isequal(index_of_longest,3));

fprintf(1,'The longest path of the %.0d paths was path %.0d with %.0d elements\n',length(cellArrayOfPaths),index_of_longest,length(cellArrayOfPaths{index_of_longest}(:,1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

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
fprintf(1,'Figure: %.0f: Demo of fast mode, empty fig_num\n',fig_num);
figure(fig_num); close(fig_num);

% Call the function that fills in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
clear cellArrayOfPaths; cellArrayOfPaths = cell(3,1);
for i_Path = 1:length(paths)
    
    cellArrayOfPaths{i_Path,1} = paths{i_Path};
end

% Call the function
index_of_longest = fcn_Path_findPathWithMostData(cellArrayOfPaths, ([]));

% Check variable types
assert(isnumeric(index_of_longest));

% Check variable sizes
assert(isequal(size(index_of_longest),[1 1]));

% Check variable values
assert(isequal(index_of_longest,3));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);


% Call the function that fills in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
clear cellArrayOfPaths; cellArrayOfPaths = cell(3,1);
for i_Path = 1:length(paths)
    
    cellArrayOfPaths{i_Path,1} = paths{i_Path};
end

% Call the function
index_of_longest = fcn_Path_findPathWithMostData(cellArrayOfPaths, (-1));

% Check variable types
assert(isnumeric(index_of_longest));

% Check variable sizes
assert(isequal(size(index_of_longest),[1 1]));

% Check variable values
assert(isequal(index_of_longest,3));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

% Call the function that fills in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
clear cellArrayOfPaths; cellArrayOfPaths = cell(3,1);
for i_Path = 1:length(paths)
    
    cellArrayOfPaths{i_Path,1} = paths{i_Path};
end

Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    index_of_longest = fcn_Path_findPathWithMostData(cellArrayOfPaths, ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    index_of_longest = fcn_Path_findPathWithMostData(cellArrayOfPaths, (-1));
end
fast_method = toc;

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

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

% close all;

%% BUG 



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


