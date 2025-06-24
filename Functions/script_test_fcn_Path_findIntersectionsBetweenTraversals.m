% script_test_fcn_Path_findIntersectionsBetweenTraversals.m
% This is a script to exercise the function: fcn_Path_findIntersectionsBetweenTraversals
% This function was written on 2020_11_14 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2021_01_23:
% -- first write of the code

close all;


% FORMAT:
% [intersection_points,...
%     s_coordinates_in_traversal_1,...
%     s_coordinates_in_traversal_2] = ...
%     fcn_Path_findIntersectionsBetweenTraversals(...
%     traversal_1,...
%     traversal_2, ...
%     fig_num)

%% BASIC example 1 - perpendicular lines, intersection in middle
fig_num = 10001;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;


% Create a dummy path and convert it to traversal_1
traversal_1_path = [0 0; 10 0];
traversal_1 = fcn_Path_convertPathToTraversalStructure(traversal_1_path);

% Create a dummy path and convert it to traversal_2
traversal_2_path = [5 -5; 5 5];
traversal_2 = fcn_Path_convertPathToTraversalStructure(traversal_2_path);

% Calculate the function output
[intersection_points,...
    s_coordinates_in_traversal_1,...
    s_coordinates_in_traversal_2] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    traversal_1,...
    traversal_2, ...
    fig_num);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);


%% BASIC example 2 - perpendicular lines, intersection in middle
fig_num = 10002;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;


% Create a dummy path and convert it to traversal_1
traversal_1_path = [0 0; 10 0];
traversal_1 = fcn_Path_convertPathToTraversalStructure(traversal_1_path);

% Create a dummy path and convert it to traversal_2
traversal_2_path = [5 -1; 5 1];
traversal_2 = fcn_Path_convertPathToTraversalStructure(traversal_2_path);

% Calculate the function output
[intersection_points,...
    s_coordinates_in_traversal_1,...
    s_coordinates_in_traversal_2] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    traversal_1,...
    traversal_2, ...
    fig_num);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);


%% BASIC example 3 - parallel lines, no intersection
fig_num = 3;

% Create a dummy path and convert it to traversal_1
traversal_1_path = [0 0; 10 0];
traversal_1 = fcn_Path_convertPathToTraversalStructure(traversal_1_path);

% Create a dummy path and convert it to traversal_2
traversal_2_path = [0 5; 10 5];
traversal_2 = fcn_Path_convertPathToTraversalStructure(traversal_2_path);


% Calculate the function output
[intersection_points,...
    s_coordinates_in_traversal_1,...
    s_coordinates_in_traversal_2] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    traversal_1,...
    traversal_2, ...
    fig_num);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);


%% BASIC example 4 - two crossings of 2 onto 1
fig_num = 4;

% Create a dummy path and convert it to traversal_1
traversal_1_path = [0 0; 7 2];
traversal_1 = fcn_Path_convertPathToTraversalStructure(traversal_1_path);

% Create a dummy path and convert it to traversal_2
traversal_2_path = [0 -1; 2 2;  5 -3];
traversal_2 = fcn_Path_convertPathToTraversalStructure(traversal_2_path);


% Calculate the function output
[intersection_points,...
    s_coordinates_in_traversal_1,...
    s_coordinates_in_traversal_2] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    traversal_1,...
    traversal_2, ...
    fig_num);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);


%% BASIC example 3 - two crossings of 1 onto 2
fig_num = 10003;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;


% Create a dummy path and convert it to traversal_1
traversal_1_path = [0 -1; 2 2;  5 -3];
traversal_1 = fcn_Path_convertPathToTraversalStructure(traversal_1_path);

% Create a dummy path and convert it to traversal_2
traversal_2_path = [0 0; 7 2];
traversal_2 = fcn_Path_convertPathToTraversalStructure(traversal_2_path);


% Calculate the function output
[intersection_points,...
    s_coordinates_in_traversal_1,...
    s_coordinates_in_traversal_2] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    traversal_1,...
    traversal_2, ...
    fig_num);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);


%% BASIC example - intersection at vertex of one traversal
fig_num = 10004;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Create a dummy path and convert it to traversal_1
traversal_1_path = [0 2; 7 2];
traversal_1 = fcn_Path_convertPathToTraversalStructure(traversal_1_path);

% Create a dummy path and convert it to traversal_2
traversal_2_path = [0 -1; 2 2;  5 -3];
traversal_2 = fcn_Path_convertPathToTraversalStructure(traversal_2_path);


% Calculate the function output
[intersection_points,...
    s_coordinates_in_traversal_1,...
    s_coordinates_in_traversal_2] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    traversal_1,...
    traversal_2, ...
    fig_num);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);


%% BASIC example - intersection at vertex of both traversals
fig_num = 10005;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Create a dummy path and convert it to traversal_1
traversal_1_path = [2 0; 2 2; 7 2];
traversal_1 = fcn_Path_convertPathToTraversalStructure(traversal_1_path);

% Create a dummy path and convert it to traversal_2
traversal_2_path = [0 -1; 2 2;  5 -3];
traversal_2 = fcn_Path_convertPathToTraversalStructure(traversal_2_path);


% Calculate the function output
[intersection_points,...
    s_coordinates_in_traversal_1,...
    s_coordinates_in_traversal_2] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    traversal_1,...
    traversal_2, ...
    fig_num);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);


%% BASIC example - intersection at both ends of one traversal
fig_num = 10006;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;


% Create a dummy path and convert it to traversal_1
traversal_1_path = [1 0.5; 3 0.5];
traversal_1 = fcn_Path_convertPathToTraversalStructure(traversal_1_path);

% Create a dummy path and convert it to traversal_2
traversal_2_path = [0 -1; 2 2;  4 -1];
traversal_2 = fcn_Path_convertPathToTraversalStructure(traversal_2_path);


% Calculate the function output
[intersection_points,...
    s_coordinates_in_traversal_1,...
    s_coordinates_in_traversal_2] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    traversal_1,...
    traversal_2, ...
    fig_num);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);


%% BASIC example - intersection at both ends of one traversal with two vertices
fig_num = 10007;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Create a dummy path and convert it to traversal_1
traversal_1_path = [0 1; 1 0.5; 3 0.5; 3.5 1];
traversal_1 = fcn_Path_convertPathToTraversalStructure(traversal_1_path);

% Create a dummy path and convert it to traversal_2
traversal_2_path = [0 -1; 2 2;  4 -1];
traversal_2 = fcn_Path_convertPathToTraversalStructure(traversal_2_path);


% Calculate the function output
[intersection_points,...
    s_coordinates_in_traversal_1,...
    s_coordinates_in_traversal_2] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    traversal_1,...
    traversal_2, ...
    fig_num);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);

%% BASIC example - Loop of traversal 1 over same point twice
fig_num = 10008;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Create a dummy path and convert it to traversal_1
traversal_1_path = [0 1; 2 0; 2 1; 0 0];
traversal_1 = fcn_Path_convertPathToTraversalStructure(traversal_1_path);

% Create a dummy path and convert it to traversal_2
traversal_2_path = [0 -1; 2 2;  4 -1];
traversal_2 = fcn_Path_convertPathToTraversalStructure(traversal_2_path);


% Calculate the function output
[intersection_points,...
    s_coordinates_in_traversal_1,...
    s_coordinates_in_traversal_2] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    traversal_1,...
    traversal_2, ...
    fig_num);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);


%% BASIC example - Loop of traversal 2 over same point twice
fig_num = 10009;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;


% Create a dummy path and convert it to traversal_1
traversal_1_path = [0 -1; 2 2;  4 -1];
traversal_1 = fcn_Path_convertPathToTraversalStructure(traversal_1_path);

% Create a dummy path and convert it to traversal_2
traversal_2_path = [0 1; 2 0; 2 1; 0 0];
traversal_2 = fcn_Path_convertPathToTraversalStructure(traversal_2_path);


% Calculate the function output
[intersection_points,...
    s_coordinates_in_traversal_1,...
    s_coordinates_in_traversal_2] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    traversal_1,...
    traversal_2, ...
    fig_num);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);


%% BASIC example - Traversal 2 is on top of traversal 1 for an area
fig_num = 10010;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;


% Create a dummy path and convert it to traversal_1
traversal_1_path = [0 0; 8 2];
traversal_1 = fcn_Path_convertPathToTraversalStructure(traversal_1_path);

% Create a dummy path and convert it to traversal_2
traversal_2_path = [1 1; 2 0.5; 4 1; 5 3];
traversal_2 = fcn_Path_convertPathToTraversalStructure(traversal_2_path);


% Calculate the function output
[intersection_points,...
    s_coordinates_in_traversal_1,...
    s_coordinates_in_traversal_2] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    traversal_1,...
    traversal_2, ...
    fig_num);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);

%% BASIC example - Traversal 1 is on top of traversal 2 for an area
fig_num = 10011;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Create a dummy path and convert it to traversal_1
traversal_1_path = [1 1; 2 0.5; 4 1; 5 3];
traversal_1 = fcn_Path_convertPathToTraversalStructure(traversal_1_path);

% Create a dummy path and convert it to traversal_2
traversal_2_path = [0 0; 8 2];
traversal_2 = fcn_Path_convertPathToTraversalStructure(traversal_2_path);


% Calculate the function output
[intersection_points,...
    s_coordinates_in_traversal_1,...
    s_coordinates_in_traversal_2] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    traversal_1,...
    traversal_2, ...
    fig_num);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);


%% BASIC example - Traversal 2 is on top of traversal 1 for two areas
fig_num = 10012;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;


% Create a dummy path and convert it to traversal_1
traversal_1_path = [0 0; 8 2];
traversal_1 = fcn_Path_convertPathToTraversalStructure(traversal_1_path);

% Create a dummy path and convert it to traversal_2
traversal_2_path = [1 1; 2 0.5; 4 1; 5 3; 6 1.5; 7 1.75; 9 1];
traversal_2 = fcn_Path_convertPathToTraversalStructure(traversal_2_path);


% Calculate the function output
[intersection_points,...
    s_coordinates_in_traversal_1,...
    s_coordinates_in_traversal_2] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    traversal_1,...
    traversal_2, ...
    fig_num);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);

%% BASIC example - Traversal 1 is on top of traversal 2 for two areas
fig_num = 10013;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;


% Create a dummy path and convert it to traversal_1
traversal_1_path = [1 1; 2 0.5; 4 1; 5 3; 6 1.5; 7 1.75; 9 1];
traversal_1 = fcn_Path_convertPathToTraversalStructure(traversal_1_path);

% Create a dummy path and convert it to traversal_2
traversal_2_path = [0 0; 8 2];
traversal_2 = fcn_Path_convertPathToTraversalStructure(traversal_2_path);


% Calculate the function output
[intersection_points,...
    s_coordinates_in_traversal_1,...
    s_coordinates_in_traversal_2] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    traversal_1,...
    traversal_2, ...
    fig_num);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);

%% BASIC example - Traversal 1 is same as Traversal 2 for many segments
fig_num = 10014;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;


% Create a dummy path and convert it to traversal_1
traversal_1_path = [1 1; 2 0.5; 4 1; 5 3; 6 1.5; 7 1.75; 9 1];
traversal_1 = fcn_Path_convertPathToTraversalStructure(traversal_1_path);

% Create a dummy path and convert it to traversal_2
traversal_2_path = [1 1; 2 0.5; 4 1; 5 3; 6 1.5; 7 1.75; 9 1; 10 0];
traversal_2 = fcn_Path_convertPathToTraversalStructure(traversal_2_path);


% Calculate the function output
[intersection_points,...
    s_coordinates_in_traversal_1,...
    s_coordinates_in_traversal_2] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    traversal_1,...
    traversal_2, ...
    fig_num);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);



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

% Create a dummy path and convert it to traversal_1
traversal_1_path = [0 0; 7 2];
traversal_1 = fcn_Path_convertPathToTraversalStructure(traversal_1_path);

% Create a dummy path and convert it to traversal_2
traversal_2_path = [0 -1; 2 2;  5 -3];
traversal_2 = fcn_Path_convertPathToTraversalStructure(traversal_2_path);


% Calculate the function output
[intersection_points,...
    s_coordinates_in_traversal_1,...
    s_coordinates_in_traversal_2] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    traversal_1,...
    traversal_2, ...
    []);

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

% Create a dummy path and convert it to traversal_1
traversal_1_path = [0 0; 7 2];
traversal_1 = fcn_Path_convertPathToTraversalStructure(traversal_1_path);

% Create a dummy path and convert it to traversal_2
traversal_2_path = [0 -1; 2 2;  5 -3];
traversal_2 = fcn_Path_convertPathToTraversalStructure(traversal_2_path);


% Calculate the function output
[intersection_points,...
    s_coordinates_in_traversal_1,...
    s_coordinates_in_traversal_2] = ...
    fcn_Path_findIntersectionsBetweenTraversals(...
    traversal_1,...
    traversal_2, ...
    -1);

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

% Create a dummy path and convert it to traversal_1
traversal_1_path = [0 0; 7 2];
traversal_1 = fcn_Path_convertPathToTraversalStructure(traversal_1_path);

% Create a dummy path and convert it to traversal_2
traversal_2_path = [0 -1; 2 2;  5 -3];
traversal_2 = fcn_Path_convertPathToTraversalStructure(traversal_2_path);


Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [intersection_points,...
        s_coordinates_in_traversal_1,...
        s_coordinates_in_traversal_2] = ...
        fcn_Path_findIntersectionsBetweenTraversals(...
        traversal_1,...
        traversal_2, ...
        []);
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [intersection_points,...
        s_coordinates_in_traversal_1,...
        s_coordinates_in_traversal_2] = ...
        fcn_Path_findIntersectionsBetweenTraversals(...
        traversal_1,...
        traversal_2, ...
        -1);
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
function print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2)
N_chars = 25;


% Print the header
header_1_str = sprintf('%s','Intersection ID');
fixed_header_1_str = fcn_DebugTools_debugPrintStringToNCharacters(header_1_str,N_chars);
header_2_str = sprintf('%s','Location X');
fixed_header_2_str = fcn_DebugTools_debugPrintStringToNCharacters(header_2_str,N_chars);
header_3_str = sprintf('%s','Location Y');
fixed_header_3_str = fcn_DebugTools_debugPrintStringToNCharacters(header_3_str,N_chars);
header_4_str = sprintf('%s','s-coord in traversal 1');
fixed_header_4_str = fcn_DebugTools_debugPrintStringToNCharacters(header_4_str,N_chars);
header_5_str = sprintf('%s','s-coord in traversal 2');
fixed_header_5_str = fcn_DebugTools_debugPrintStringToNCharacters(header_5_str,N_chars);

fprintf(1,'\n\n%s %s %s %s %s\n',...
    fixed_header_1_str,...
    fixed_header_2_str,...
    fixed_header_3_str,...
    fixed_header_4_str,...
    fixed_header_5_str);

% Print the results
if ~isempty(intersection_points)
    for ith_intersection =1:length(intersection_points(:,1))
        results_1_str = sprintf('%.0d',ith_intersection);
        fixed_results_1_str = fcn_DebugTools_debugPrintStringToNCharacters(results_1_str,N_chars);
        results_2_str = sprintf('%.2f',intersection_points(ith_intersection,1));
        fixed_results_2_str = fcn_DebugTools_debugPrintStringToNCharacters(results_2_str,N_chars);
        results_3_str = sprintf('%.2f',intersection_points(ith_intersection,2));
        fixed_results_3_str = fcn_DebugTools_debugPrintStringToNCharacters(results_3_str,N_chars);
        results_4_str = sprintf('%.2f',s_coordinates_in_traversal_1(ith_intersection));
        fixed_results_4_str = fcn_DebugTools_debugPrintStringToNCharacters(results_4_str,N_chars);
        results_5_str = sprintf('%.2f',s_coordinates_in_traversal_2(ith_intersection));
        fixed_results_5_str = fcn_DebugTools_debugPrintStringToNCharacters(results_5_str,N_chars);

        fprintf(1,'%s %s %s %s %s\n',...
            fixed_results_1_str,...
            fixed_results_2_str,...
            fixed_results_3_str,...
            fixed_results_4_str,...
            fixed_results_5_str);

        %         fprintf(1,'\t\t%.0d \t\t\t\t %.2f \t\t\t %.2f \t\t\t\t %.2f \t\t\t\t\t %.2f\n',...
        %             ith_intersection,...
        %             intersection_points(ith_intersection,1),...
        %             intersection_points(ith_intersection,2),...
        %             s_coordinates_in_traversal_1(ith_intersection),...
        %             s_coordinates_in_traversal_2(ith_intersection));

    end % Ends for loop
end % Ends check to see if isempty
end % Ends function

