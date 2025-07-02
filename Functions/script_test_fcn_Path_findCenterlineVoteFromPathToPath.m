% script_test_fcn_Path_findCenterlineVoteFromPathToPath
% Tests the following:
%    [centerline_points_projected,unit_vectors_orthogonal] = ...
%     fcn_Path_findCenterlineVoteFromPathToPath(...
%     from_path,to_path,(flag_rounding_type),(search_radius),(fig_num))

% Revision history:
% 2023_09_04 by S. Brennan
% -- first write of the code

close all;


%% Basic demonstration of fcn_Path_findCenterlineVoteFromPathToPath
fig_num = 10001;
titleString = sprintf('Finds the center projected from one path toward another');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

from_path = [0 0; 1 1; 2 1; 3 4];
to_path   = from_path + ones(length(from_path(:,1)),1)*[0 1];


flag_rounding_type = 1;
search_radius = 10;
flag_project_full_distance = 0;

[centerline_points_projected,unit_vectors_orthogonal] = ...
    fcn_Path_findCenterlineVoteFromPathToPath(...
    from_path,to_path,(flag_rounding_type),(search_radius),(flag_project_full_distance), (fig_num));
title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(centerline_points_projected));
assert(isnumeric(unit_vectors_orthogonal));

% Check variable sizes
Nfrom = length(from_path(:,1));
assert(isequal(size(centerline_points_projected),[Nfrom 2]));
assert(isequal(size(unit_vectors_orthogonal),[Nfrom 2]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Basic demonstration of fcn_Path_findCenterlineVoteFromPathToPath
fig_num = 10002;
titleString = sprintf('Shows that if the serach distance is too small, nothing is returned');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

from_path = [0 0; 1 1; 2 1; 3 4];
to_path   = from_path + ones(length(from_path(:,1)),1)*[0 1];


flag_rounding_type = 1;
search_radius = 0.1;
flag_project_full_distance = 0;

[centerline_points_projected,unit_vectors_orthogonal] = ...
    fcn_Path_findCenterlineVoteFromPathToPath(...
    from_path,to_path,(flag_rounding_type),(search_radius),(flag_project_full_distance), (fig_num));
title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(centerline_points_projected));
assert(isnumeric(unit_vectors_orthogonal));

% Check variable sizes
Nfrom = length(from_path(:,1));
assert(isequal(size(centerline_points_projected),[Nfrom 2]));
assert(isequal(size(unit_vectors_orthogonal),[Nfrom 2]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Basic demonstration 3 of fcn_Path_findCenterlineVoteFromPathToPath
fig_num = 10003;
titleString = sprintf('Swap from and to paths to show how that full projection returns mapping of 1 onto 2');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;


from_path = [0 0; 1 1; 2 1; 3 4];
to_path   = from_path + ones(length(from_path(:,1)),1)*[0 1];


flag_rounding_type = 1;
search_radius = 10;
flag_project_full_distance = 1;

[centerline_points_projected,unit_vectors_orthogonal] = ...
    fcn_Path_findCenterlineVoteFromPathToPath(...
    from_path,to_path,(flag_rounding_type),(search_radius),(flag_project_full_distance), (fig_num));
title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(centerline_points_projected));
assert(isnumeric(unit_vectors_orthogonal));

% Check variable sizes
Nfrom = length(from_path(:,1));
assert(isequal(size(centerline_points_projected),[Nfrom 2]));
assert(isequal(size(unit_vectors_orthogonal),[Nfrom 2]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Call the center calculation function
fig_num = 10004;
titleString = sprintf('Shows the center calculation for realistic path');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Fill in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;

from_path =  paths{1};
to_path =  paths{2};
flag_rounding_type = 1;
search_radius = 10;
flag_project_full_distance = 0;

[centerline_points_projected,unit_vectors_orthogonal] = ...
    fcn_Path_findCenterlineVoteFromPathToPath(...
    from_path,to_path,(flag_rounding_type),(search_radius),(flag_project_full_distance), (fig_num));
title(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(centerline_points_projected));
assert(isnumeric(unit_vectors_orthogonal));

% Check variable sizes
Nfrom = length(from_path(:,1));
assert(isequal(size(centerline_points_projected),[Nfrom 2]));
assert(isequal(size(unit_vectors_orthogonal),[Nfrom 2]));

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

from_path = [0 0; 1 1; 2 1; 3 4];
to_path   = from_path + ones(length(from_path(:,1)),1)*[0 1];


flag_rounding_type = 1;
search_radius = 0.1;
flag_project_full_distance = 0;

[centerline_points_projected,unit_vectors_orthogonal] = ...
    fcn_Path_findCenterlineVoteFromPathToPath(...
    from_path,to_path,(flag_rounding_type),(search_radius),(flag_project_full_distance), ([]));

% Check variable types
assert(isnumeric(centerline_points_projected));
assert(isnumeric(unit_vectors_orthogonal));

% Check variable sizes
Nfrom = length(from_path(:,1));
assert(isequal(size(centerline_points_projected),[Nfrom 2]));
assert(isequal(size(unit_vectors_orthogonal),[Nfrom 2]));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

from_path = [0 0; 1 1; 2 1; 3 4];
to_path   = from_path + ones(length(from_path(:,1)),1)*[0 1];


flag_rounding_type = 1;
search_radius = 0.1;
flag_project_full_distance = 0;

[centerline_points_projected,unit_vectors_orthogonal] = ...
    fcn_Path_findCenterlineVoteFromPathToPath(...
    from_path,to_path,(flag_rounding_type),(search_radius),(flag_project_full_distance), (-1));

% Check variable types
assert(isnumeric(centerline_points_projected));
assert(isnumeric(unit_vectors_orthogonal));

% Check variable sizes
Nfrom = length(from_path(:,1));
assert(isequal(size(centerline_points_projected),[Nfrom 2]));
assert(isequal(size(unit_vectors_orthogonal),[Nfrom 2]));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

from_path = [0 0; 1 1; 2 1; 3 4];
to_path   = from_path + ones(length(from_path(:,1)),1)*[0 1];


flag_rounding_type = 1;
search_radius = 0.1;
flag_project_full_distance = 0;

Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [centerline_points_projected,unit_vectors_orthogonal] = ...
        fcn_Path_findCenterlineVoteFromPathToPath(...
        from_path,to_path,(flag_rounding_type),(search_radius),(flag_project_full_distance), ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [centerline_points_projected,unit_vectors_orthogonal] = ...
        fcn_Path_findCenterlineVoteFromPathToPath(...
        from_path,to_path,(flag_rounding_type),(search_radius),(flag_project_full_distance), (-1));
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