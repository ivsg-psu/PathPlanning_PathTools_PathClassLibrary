% script_test_fcn_Path_findCenterPathBetweenTwoPaths
% Tests the following:
%    center_path = ...
%     fcn_Path_findCenterPathBetweenTwoPaths(...
%     first_path,second_path,(flag_rounding_type),(search_radius),(fig_num))

% Revision history:
% 2023_09_04 by S. Brennan
% -- first write of the code

close all;



%% Basic demo 
fig_num = 10001;
titleString = sprintf('Basic demonstration of fcn_Path_findCenterPathBetweenTwoPaths');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

first_path = [0 0; 1 1; 2 1; 3 2];
second_path   = [0.5 1.5; 1.5 2.1; 4 6];

flag_rounding_type = 1;
search_radius = 10;

center_path = ...
    fcn_Path_findCenterPathBetweenTwoPaths(...
    first_path,second_path,(flag_rounding_type),(search_radius),(fig_num));
title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(center_path));

% Check variable sizes
assert(length(center_path(:,1)) == (length(first_path(:,1)) + length(second_path(:,1))));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Demonstration of fcn_Path_findCenterPathBetweenTwoPaths
fig_num = 10002;
titleString = sprintf('Slightly different example');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

first_path = [0 0; 1 1; 2 1; 3 4];
second_path   = first_path + ones(length(first_path(:,1)),1)*[0 1];

flag_rounding_type = 1;
search_radius = 10;

center_path = ...
    fcn_Path_findCenterPathBetweenTwoPaths(...
    first_path,second_path,(flag_rounding_type),(search_radius),(fig_num));
title(titleString, 'Interpreter','none');


% Check variable types
assert(isnumeric(center_path));

% Check variable sizes
assert(length(center_path(:,1)) == (length(first_path(:,1)) + length(second_path(:,1))));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Demonstration of effect of flag_rounding_type
titleString = sprintf('effect of flag_rounding_type');
fig_num = 10003;
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

first_path = [0 0; 1 1; 2 1; 3 2];
second_path   = [0.5 1.5; 1.5 2.1; 4 6];

search_radius = 10;


for flag_rounding_type = 1:3
    subplot(2,2,flag_rounding_type);
    hold on;
    grid on;
    axis equal;
    center_path = ...
        fcn_Path_findCenterPathBetweenTwoPaths(...
        first_path,second_path,(flag_rounding_type),(search_radius),(fig_num));
    title(sprintf('flag_rounding_type: %.0d',flag_rounding_type),'Interpreter','none');
end
sgtitle(titleString,'Interpreter','none');

% Check variable types
assert(isnumeric(center_path));

% Check variable sizes
assert(length(center_path(:,1)) == (length(first_path(:,1)) + length(second_path(:,1))));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Test case where first path is not intersecting second path with any projections
fig_num = 10004;
titleString = sprintf('first path is not intersecting second path with any projections');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

first_path = [0 0; 10 0];
second_path   = [1 1; 9 1];

flag_rounding_type = 1;
search_radius = 10;

center_path = ...
    fcn_Path_findCenterPathBetweenTwoPaths(...
    first_path,second_path,(flag_rounding_type),(search_radius),(fig_num));
title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(center_path));

% Check variable sizes
assert(length(center_path(:,1)) == (length(first_path(:,1)) + length(second_path(:,1))));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Test case where second path is not intersecting first path with any projections
fig_num = 10005;
titleString = sprintf('second path is not intersecting first path with any projections');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

second_path = [0 0; 10 0];
first_path   = [1 1; 9 1];

flag_rounding_type = 1;
search_radius = 10;

center_path = ...
    fcn_Path_findCenterPathBetweenTwoPaths(...
    first_path,second_path,(flag_rounding_type),(search_radius),(fig_num));
title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(center_path));

% Check variable sizes
assert(length(center_path(:,1)) == (length(first_path(:,1)) + length(second_path(:,1))));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Test case where first path and second path splay
fig_num = 10006;
titleString = sprintf('first path and second path splay');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

first_path = [-5 0; 0 0; 10 0; 30 0];
second_path   = [-5 5; 1 5; 11 11; 30 11];

flag_rounding_type = 1;
search_radius = 20;

center_path = ...
    fcn_Path_findCenterPathBetweenTwoPaths(...
    first_path,second_path,(flag_rounding_type),(search_radius),(fig_num));
title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(center_path));

% Check variable sizes
assert(length(center_path(:,1)) == (length(first_path(:,1)) + length(second_path(:,1))));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Call the center calculation using a real-world path
fig_num = 10007;
titleString = sprintf('center calculation using a real-world path');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

paths = fcn_Path_fillSamplePaths;
first_path = paths{1};
second_path = paths{2};

flag_rounding_type = 1;
search_radius = 10;

center_path = ...
    fcn_Path_findCenterPathBetweenTwoPaths(...
    first_path,second_path,(flag_rounding_type),(search_radius),(fig_num));
title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(center_path));

% Check variable sizes
assert(length(center_path(:,1)) == (length(first_path(:,1)) + length(second_path(:,1))));

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

first_path = [0 0; 1 1; 2 1; 3 2];
second_path   = [0.5 1.5; 1.5 2.1; 4 6];

flag_rounding_type = 1;
search_radius = 10;

center_path = ...
    fcn_Path_findCenterPathBetweenTwoPaths(...
    first_path,second_path,(flag_rounding_type),(search_radius),([]));

% Check variable types
assert(isnumeric(center_path));

% Check variable sizes
assert(length(center_path(:,1)) == (length(first_path(:,1)) + length(second_path(:,1))));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

first_path = [0 0; 1 1; 2 1; 3 2];
second_path   = [0.5 1.5; 1.5 2.1; 4 6];

flag_rounding_type = 1;
search_radius = 10;

center_path = ...
    fcn_Path_findCenterPathBetweenTwoPaths(...
    first_path,second_path,(flag_rounding_type),(search_radius),(-1));

% Check variable types
assert(isnumeric(center_path));

% Check variable sizes
assert(length(center_path(:,1)) == (length(first_path(:,1)) + length(second_path(:,1))));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

first_path = [0 0; 1 1; 2 1; 3 2];
second_path   = [0.5 1.5; 1.5 2.1; 4 6];

flag_rounding_type = 1;
search_radius = 10;

Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    center_path = ...
        fcn_Path_findCenterPathBetweenTwoPaths(...
        first_path,second_path,(flag_rounding_type),(search_radius),([]));

end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    center_path = ...
        fcn_Path_findCenterPathBetweenTwoPaths(...
        first_path,second_path,(flag_rounding_type),(search_radius),(-1));
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


%% FAIL CASES
if 1==0
    %% FAIL CASE: neither path has projections that intersect each other
    fig_num = 90001;
    titleString = sprintf('FAIL CASE: neither path has projections that intersect each other');
    fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
    figure(fig_num); clf;

    first_path = [0 0; 10 0];
    second_path   = [-1 0; -10 10];

    flag_rounding_type = 1;
    search_radius = 10;

    center_path = ...
        fcn_Path_findCenterPathBetweenTwoPaths(...
        first_path,second_path,(flag_rounding_type),(search_radius),(fig_num));
    title(titleString, 'Interpreter','none');

end

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
