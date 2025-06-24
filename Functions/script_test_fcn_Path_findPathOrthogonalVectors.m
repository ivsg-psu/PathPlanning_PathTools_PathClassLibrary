% script_test_fcn_Path_findPathOrthogonalVectors
% This is a script to exercise the function: fcn_Path_findPathOrthogonalVectors
% This function was written on 2023_08_27 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history
%     2023_08_27
%     -- first write of the code

close all;

%% BASIC example 1 - simple horizontal line, default flag (1), no figure
fig_num = 10001;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

path = [0 0; 4 0]; 

% Call the function
[normal_unit_vectors_at_midpoints, normal_unit_vectors_at_joints] = ...
    fcn_Path_findPathOrthogonalVectors(path, [], fig_num);

% Make sure function worked
assert(isequal(round(normal_unit_vectors_at_midpoints,4),[0 1]));
assert(isequal(round(normal_unit_vectors_at_joints,4),[0 1; 0 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example - simple horizontal line, default flag (1), with figure
fig_num = 10002;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

path = [0 0; 4 0];  
flag_rounding_type = 1; % Define the rounding type

% Call the function
[normal_unit_vectors_at_midpoints, normal_unit_vectors_at_joints] = ...
    fcn_Path_findPathOrthogonalVectors(path,flag_rounding_type, fig_num);

% Make sure function worked
assert(isequal(round(normal_unit_vectors_at_midpoints,4),[0 1]));
assert(isequal(round(normal_unit_vectors_at_joints,4),[0 1; 0 1]));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC example - right-angled line segment - flag tests
fig_num = 10003;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;


path = [0 0; 2 0; 2 -2];

% Flag = 1
subplot(1,4,1);
flag_rounding_type = 1; % Define the rounding type

% Call the function
[normal_unit_vectors_at_midpoints, normal_unit_vectors_at_joints] = ...
    fcn_Path_findPathOrthogonalVectors(path,flag_rounding_type, fig_num);
title('flag_rounding_type = 1','Interpreter','none');

% Make sure function worked
assert(isequal(round(normal_unit_vectors_at_midpoints,4),[0 1; 1 0]));
assert(isequal(round(normal_unit_vectors_at_joints,4),[0 1; 0 1; 1 0]));

% Flag = 2
subplot(1,4,2);
flag_rounding_type = 2; % Define the rounding type

% Call the function
[normal_unit_vectors_at_midpoints, normal_unit_vectors_at_joints] = ...
    fcn_Path_findPathOrthogonalVectors(path,flag_rounding_type, fig_num);
title('flag_rounding_type = 2','Interpreter','none');

% Make sure function worked
assert(isequal(round(normal_unit_vectors_at_midpoints,4),[0 1; 1 0]));
assert(isequal(round(normal_unit_vectors_at_joints,4),[0 1; 1 0; 1 0]));

% Flag = 3
subplot(1,4,3);
flag_rounding_type = 3; % Define the rounding type

% Call the function
[normal_unit_vectors_at_midpoints, normal_unit_vectors_at_joints] = ...
    fcn_Path_findPathOrthogonalVectors(path,flag_rounding_type, fig_num);
title('flag_rounding_type = 3','Interpreter','none');

% Make sure function worked
assert(isequal(round(normal_unit_vectors_at_midpoints,4),[0 1; 1 0]));
assert(isequal(round(normal_unit_vectors_at_joints,4),[0 1; round([1 1]./(2^0.5),4); 1 0]));

% Flag = 4
subplot(1,4,4);
flag_rounding_type = 4; % Define the rounding type

% Call the function
[normal_unit_vectors_at_midpoints, normal_unit_vectors_at_joints] = ...
    fcn_Path_findPathOrthogonalVectors(path,flag_rounding_type, fig_num);
title('flag_rounding_type = 4','Interpreter','none');

% Make sure function worked
assert(isequal(round(normal_unit_vectors_at_midpoints,4),[0 1; 1 0]));
assert(isequal(round(normal_unit_vectors_at_joints,4),[-1 0; round([1 1]./(2^0.5),4); 0 -1]));

sgtitle('Function: fcn_Path_findPathOrthogonalVectors, showing effect of flag_rounding_type','Interpreter','none')


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Real path examples
fig_num = 10004;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

flag_rounding_type = 1; % Define the rounding type

% Fill in some dummy data
paths = fcn_Path_fillSamplePaths;


for i_Path = 1:3
    temp_fig_num = fig_num + i_Path - 1;  % Define the figure
    figure(temp_fig_num);
    clf;
    
    % Call the function
    [~, ~] = ...
        fcn_Path_findPathOrthogonalVectors(paths{i_Path},flag_rounding_type, temp_fig_num);
end



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

path = [0 0; 4 0];  
flag_rounding_type = 1; % Define the rounding type

% Call the function
[normal_unit_vectors_at_midpoints, normal_unit_vectors_at_joints] = ...
    fcn_Path_findPathOrthogonalVectors(path,flag_rounding_type, []);

% Make sure function worked
assert(isequal(round(normal_unit_vectors_at_midpoints,4),[0 1]));
assert(isequal(round(normal_unit_vectors_at_joints,4),[0 1; 0 1]));


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

path = [0 0; 4 0];  
flag_rounding_type = 1; % Define the rounding type

% Call the function
[normal_unit_vectors_at_midpoints, normal_unit_vectors_at_joints] = ...
    fcn_Path_findPathOrthogonalVectors(path,flag_rounding_type, -1);

% Make sure function worked
assert(isequal(round(normal_unit_vectors_at_midpoints,4),[0 1]));
assert(isequal(round(normal_unit_vectors_at_joints,4),[0 1; 0 1]));


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

path = [0 0; 4 0];  
flag_rounding_type = 1; % Define the rounding type

Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [normal_unit_vectors_at_midpoints, normal_unit_vectors_at_joints] = ...
        fcn_Path_findPathOrthogonalVectors(path,flag_rounding_type, []);

end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [normal_unit_vectors_at_midpoints, normal_unit_vectors_at_joints] = ...
        fcn_Path_findPathOrthogonalVectors(path,flag_rounding_type, -1);
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
