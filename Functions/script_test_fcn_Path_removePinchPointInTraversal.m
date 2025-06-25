% script_test_fcn_Path_removePinchPointInTraversal.m
% Tests fcn_Path_removePinchPointInTraversal
       
% Revision history:
% 2021_01_24 - S. Brennan
% -- first write of the code

close all

%% Simple test case showing removal of a pinch point
fig_num = 10001;
titleString = sprintf('Simple test case showing removal of a pinch point');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

path = [0 0; 0 4; 2 4; -1 1];
traversal_with_pinch_point = fcn_Path_convertPathToTraversalStructure(path);
[traversal_no_pinch_point] = ...
    fcn_Path_removePinchPointInTraversal(...
    traversal_with_pinch_point,...
    fig_num); 

title(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(traversal_no_pinch_point));
assert(isfield(traversal_no_pinch_point,'X'))
assert(isfield(traversal_no_pinch_point,'Y'))
assert(isfield(traversal_no_pinch_point,'Z'))
assert(isfield(traversal_no_pinch_point,'Diff'))
assert(isfield(traversal_no_pinch_point,'Station'))
assert(isfield(traversal_no_pinch_point,'Yaw'))

% Check variable sizes
NreferencePoints = 3;
assert(isequal(size(traversal_no_pinch_point.X),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Y),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Z),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Diff),[NreferencePoints 2]));
assert(isequal(size(traversal_no_pinch_point.Station),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Yaw),[NreferencePoints-1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Simple test case producing warning
fig_num = 10002;
titleString = sprintf('Simple test case producing warning');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

path = [0 0; 0 4; 4 4; -1 1];
traversal_with_pinch_point = fcn_Path_convertPathToTraversalStructure(path);
[traversal_no_pinch_point] = ...
    fcn_Path_removePinchPointInTraversal(...
    traversal_with_pinch_point,...
    fig_num);

title(titleString, 'Interpreter','none');


% Check variable types
assert(isstruct(traversal_no_pinch_point));
assert(isfield(traversal_no_pinch_point,'X'))
assert(isfield(traversal_no_pinch_point,'Y'))
assert(isfield(traversal_no_pinch_point,'Z'))
assert(isfield(traversal_no_pinch_point,'Diff'))
assert(isfield(traversal_no_pinch_point,'Station'))
assert(isfield(traversal_no_pinch_point,'Yaw'))

% Check variable sizes
NreferencePoints = 3;
assert(isequal(size(traversal_no_pinch_point.X),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Y),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Z),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Diff),[NreferencePoints 2]));
assert(isequal(size(traversal_no_pinch_point.Station),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Yaw),[NreferencePoints-1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Simple test case where there are no intersections, thus no poinch points
fig_num = 10003;
titleString = sprintf('Simple test case where there are no intersections, thus no poinch points');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

path = [0 0; 0 4; 4 5; 5 5];
traversal_with_pinch_point = fcn_Path_convertPathToTraversalStructure(path);
[traversal_no_pinch_point] = ...
    fcn_Path_removePinchPointInTraversal(...
    traversal_with_pinch_point,...
    fig_num);

title(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(traversal_no_pinch_point));
assert(isfield(traversal_no_pinch_point,'X'))
assert(isfield(traversal_no_pinch_point,'Y'))
assert(isfield(traversal_no_pinch_point,'Z'))
assert(isfield(traversal_no_pinch_point,'Diff'))
assert(isfield(traversal_no_pinch_point,'Station'))
assert(isfield(traversal_no_pinch_point,'Yaw'))

% Check variable sizes
NreferencePoints = 4;
assert(isequal(size(traversal_no_pinch_point.X),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Y),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Z),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Diff),[NreferencePoints 2]));
assert(isequal(size(traversal_no_pinch_point.Station),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Yaw),[NreferencePoints-1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Simple test case with multiple crossings
fig_num = 10004;
titleString = sprintf('Simple test case with multiple crossings');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

path = [0 0; 0 4; 2 4; -1 1; -1 3; -0.5 1];
traversal_with_pinch_point = fcn_Path_convertPathToTraversalStructure(path);
[traversal_no_pinch_point] = ...
    fcn_Path_removePinchPointInTraversal(...
    traversal_with_pinch_point,...
    fig_num);
axis([-2 4 -2 6]);

title(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(traversal_no_pinch_point));
assert(isfield(traversal_no_pinch_point,'X'))
assert(isfield(traversal_no_pinch_point,'Y'))
assert(isfield(traversal_no_pinch_point,'Z'))
assert(isfield(traversal_no_pinch_point,'Diff'))
assert(isfield(traversal_no_pinch_point,'Station'))
assert(isfield(traversal_no_pinch_point,'Yaw'))

% Check variable sizes
NreferencePoints = 4;
assert(isequal(size(traversal_no_pinch_point.X),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Y),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Z),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Diff),[NreferencePoints 2]));
assert(isequal(size(traversal_no_pinch_point.Station),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Yaw),[NreferencePoints-1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Simple test case with multiple crossings
fig_num = 10005;
titleString = sprintf('Simple test case with multiple crossings');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

path = [0 0; 0 4; 2 4; -1 1; -1 3; -0.5 1; -0.5 2.5];
traversal_with_pinch_point = fcn_Path_convertPathToTraversalStructure(path);
[traversal_no_pinch_point] = ...
    fcn_Path_removePinchPointInTraversal(...
    traversal_with_pinch_point,...
    fig_num);

title(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(traversal_no_pinch_point));
assert(isfield(traversal_no_pinch_point,'X'))
assert(isfield(traversal_no_pinch_point,'Y'))
assert(isfield(traversal_no_pinch_point,'Z'))
assert(isfield(traversal_no_pinch_point,'Diff'))
assert(isfield(traversal_no_pinch_point,'Station'))
assert(isfield(traversal_no_pinch_point,'Yaw'))

% Check variable sizes
NreferencePoints = 4;
assert(isequal(size(traversal_no_pinch_point.X),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Y),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Z),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Diff),[NreferencePoints 2]));
assert(isequal(size(traversal_no_pinch_point.Station),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Yaw),[NreferencePoints-1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Crossing back toward itself
fig_num = 10006;
titleString = sprintf('Crossing back toward itself');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

path = [0 0; 0 4; 0 2];
traversal_with_pinch_point = fcn_Path_convertPathToTraversalStructure(path);
[traversal_no_pinch_point] = ...
    fcn_Path_removePinchPointInTraversal(...
    traversal_with_pinch_point,...
    fig_num);

title(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(traversal_no_pinch_point));
assert(isfield(traversal_no_pinch_point,'X'))
assert(isfield(traversal_no_pinch_point,'Y'))
assert(isfield(traversal_no_pinch_point,'Z'))
assert(isfield(traversal_no_pinch_point,'Diff'))
assert(isfield(traversal_no_pinch_point,'Station'))
assert(isfield(traversal_no_pinch_point,'Yaw'))

% Check variable sizes
NreferencePoints = 2;
assert(isequal(size(traversal_no_pinch_point.X),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Y),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Z),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Diff),[NreferencePoints 2]));
assert(isequal(size(traversal_no_pinch_point.Station),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Yaw),[NreferencePoints-1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Crossing back toward itself and continuing
fig_num = 10007;
titleString = sprintf('Crossing back toward itself and continuing');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

path = [0 0; 0 4; 0 2; 0 5];
traversal_with_pinch_point = fcn_Path_convertPathToTraversalStructure(path);
[traversal_no_pinch_point] = ...
    fcn_Path_removePinchPointInTraversal(...
    traversal_with_pinch_point,...
    fig_num);

title(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(traversal_no_pinch_point));
assert(isfield(traversal_no_pinch_point,'X'))
assert(isfield(traversal_no_pinch_point,'Y'))
assert(isfield(traversal_no_pinch_point,'Z'))
assert(isfield(traversal_no_pinch_point,'Diff'))
assert(isfield(traversal_no_pinch_point,'Station'))
assert(isfield(traversal_no_pinch_point,'Yaw'))

% Check variable sizes
NreferencePoints = 3;
assert(isequal(size(traversal_no_pinch_point.X),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Y),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Z),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Diff),[NreferencePoints 2]));
assert(isequal(size(traversal_no_pinch_point.Station),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Yaw),[NreferencePoints-1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Crossing back toward itself and continuing elsewhere
fig_num = 10008;
titleString = sprintf('Crossing back toward itself and continuing elsewhere');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

path = [0 0; 0 4; 0 2; 1 1];
traversal_with_pinch_point = fcn_Path_convertPathToTraversalStructure(path);
[traversal_no_pinch_point] = ...
    fcn_Path_removePinchPointInTraversal(...
    traversal_with_pinch_point,...
    fig_num);

title(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(traversal_no_pinch_point));
assert(isfield(traversal_no_pinch_point,'X'))
assert(isfield(traversal_no_pinch_point,'Y'))
assert(isfield(traversal_no_pinch_point,'Z'))
assert(isfield(traversal_no_pinch_point,'Diff'))
assert(isfield(traversal_no_pinch_point,'Station'))
assert(isfield(traversal_no_pinch_point,'Yaw'))

% Check variable sizes
NreferencePoints = 3;
assert(isequal(size(traversal_no_pinch_point.X),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Y),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Z),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Diff),[NreferencePoints 2]));
assert(isequal(size(traversal_no_pinch_point.Station),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Yaw),[NreferencePoints-1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Crossing back toward itself without being in sequence
fig_num = 10009;
titleString = sprintf('Crossing back toward itself without being in sequence');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

path = [0 0; 0 4; 4 4; 0 0; -1 0];
traversal_with_pinch_point = fcn_Path_convertPathToTraversalStructure(path);
[traversal_no_pinch_point] = ...
    fcn_Path_removePinchPointInTraversal(...
    traversal_with_pinch_point,...
    fig_num);

title(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(traversal_no_pinch_point));
assert(isfield(traversal_no_pinch_point,'X'))
assert(isfield(traversal_no_pinch_point,'Y'))
assert(isfield(traversal_no_pinch_point,'Z'))
assert(isfield(traversal_no_pinch_point,'Diff'))
assert(isfield(traversal_no_pinch_point,'Station'))
assert(isfield(traversal_no_pinch_point,'Yaw'))

% Check variable sizes
NreferencePoints = 2;
assert(isequal(size(traversal_no_pinch_point.X),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Y),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Z),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Diff),[NreferencePoints 2]));
assert(isequal(size(traversal_no_pinch_point.Station),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Yaw),[NreferencePoints-1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Crossing back toward itself without being in sequence
fig_num = 10010;
titleString = sprintf('Crossing back toward itself without being in sequence');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

path = [0 0; 0 4; 4 4; 0 0; 0 2];
traversal_with_pinch_point = fcn_Path_convertPathToTraversalStructure(path);
[traversal_no_pinch_point] = ...
    fcn_Path_removePinchPointInTraversal(...
    traversal_with_pinch_point,...
    fig_num);

title(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(traversal_no_pinch_point));
assert(isfield(traversal_no_pinch_point,'X'))
assert(isfield(traversal_no_pinch_point,'Y'))
assert(isfield(traversal_no_pinch_point,'Z'))
assert(isfield(traversal_no_pinch_point,'Diff'))
assert(isfield(traversal_no_pinch_point,'Station'))
assert(isfield(traversal_no_pinch_point,'Yaw'))

% Check variable sizes
NreferencePoints = 2;
assert(isequal(size(traversal_no_pinch_point.X),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Y),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Z),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Diff),[NreferencePoints 2]));
assert(isequal(size(traversal_no_pinch_point.Station),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Yaw),[NreferencePoints-1 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Advanced test case: a pinch point in practice
fig_num = 20001;
titleString = sprintf('Advanced test case: a pinch point in practice');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

figure(fig_num);
clf;
axis equal

% Clear any old variables
clear all_traversals

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1});
all_traversals.traversal{1} = reference_traversal;


% Plot the results? (Note: they are plotted below as well)
if 1==0
    fig_num = 12;
    fcn_Path_plotTraversalsYaw(all_traversals,fig_num);
    fig_num = 13;
    fcn_Path_plotTraversalsXY(all_traversals,fig_num);
end

% Grab the "curve" of the path
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths_array{1}(13:20,:));
offsets = (0:1:10)'; 
flag_rounding_type = []; % Use default
offset_traversals = fcn_Path_fillOffsetTraversalsAboutTraversal(reference_traversal, offsets, flag_rounding_type, fig_num);

% Fill in an array of "fixed" traversals
clear fixed_traversals
for ith_traversal = 1:length(offset_traversals.traversal)
    traversal_with_pinch_point = offset_traversals.traversal{ith_traversal};
    [traversal_no_pinch_point] = ...
        fcn_Path_removePinchPointInTraversal(...
        traversal_with_pinch_point);
    fixed_traversals.traversal{ith_traversal} = traversal_no_pinch_point; 
end

% Plot the results
figure(fig_num);
clf;
axis equal
hold on;
plot(reference_traversal.X,reference_traversal.Y,'b','Linewidth',3);
fcn_Path_plotTraversalsXY(fixed_traversals,fig_num);

title(titleString, 'Interpreter','none');

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

path = [0 0; 0 4; 2 4; -1 1];
traversal_with_pinch_point = fcn_Path_convertPathToTraversalStructure(path);
[traversal_no_pinch_point] = ...
    fcn_Path_removePinchPointInTraversal(...
    traversal_with_pinch_point,...
    []); 

% Check variable types
assert(isstruct(traversal_no_pinch_point));
assert(isfield(traversal_no_pinch_point,'X'))
assert(isfield(traversal_no_pinch_point,'Y'))
assert(isfield(traversal_no_pinch_point,'Z'))
assert(isfield(traversal_no_pinch_point,'Diff'))
assert(isfield(traversal_no_pinch_point,'Station'))
assert(isfield(traversal_no_pinch_point,'Yaw'))

% Check variable sizes
NreferencePoints = 3;
assert(isequal(size(traversal_no_pinch_point.X),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Y),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Z),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Diff),[NreferencePoints 2]));
assert(isequal(size(traversal_no_pinch_point.Station),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Yaw),[NreferencePoints-1 1]));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

path = [0 0; 0 4; 2 4; -1 1];
traversal_with_pinch_point = fcn_Path_convertPathToTraversalStructure(path);
[traversal_no_pinch_point] = ...
    fcn_Path_removePinchPointInTraversal(...
    traversal_with_pinch_point,...
    (-1)); 

% Check variable types
assert(isstruct(traversal_no_pinch_point));
assert(isfield(traversal_no_pinch_point,'X'))
assert(isfield(traversal_no_pinch_point,'Y'))
assert(isfield(traversal_no_pinch_point,'Z'))
assert(isfield(traversal_no_pinch_point,'Diff'))
assert(isfield(traversal_no_pinch_point,'Station'))
assert(isfield(traversal_no_pinch_point,'Yaw'))

% Check variable sizes
NreferencePoints = 3;
assert(isequal(size(traversal_no_pinch_point.X),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Y),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Z),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Diff),[NreferencePoints 2]));
assert(isequal(size(traversal_no_pinch_point.Station),[NreferencePoints 1]));
assert(isequal(size(traversal_no_pinch_point.Yaw),[NreferencePoints-1 1]));


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

path = [0 0; 0 4; 2 4; -1 1];
traversal_with_pinch_point = fcn_Path_convertPathToTraversalStructure(path);


Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [traversal_no_pinch_point] = ...
        fcn_Path_removePinchPointInTraversal(...
        traversal_with_pinch_point,...
        ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [traversal_no_pinch_point] = ...
        fcn_Path_removePinchPointInTraversal(...
        traversal_with_pinch_point,...
        (-1));
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

