% script_test_fcn_Path_snapPointOntoNearestPath.m
% This is a script to exercise the function: fcn_Path_snapPointOntoNearestPath.m
% This function was written on 2021_03_06 by Satya Prasad, szm888@psu.edu

% Revision history:


close all;

%% BASIC example 1
fig_num = 10001;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

point = [0.5 0.2];
pathXY = [0 0; 1 0; 2 0; 2 1];


[closest_path_point,s_coordinate,path_point_yaw,first_path_point_index,...
    second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,fig_num);

fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

%% BASIC example 
fig_num = 10002;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

point = [1.4 1.3]; % Define the query point
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 1.5 0.6; 3 0]; % Define an XY path


% Snap the point onto the path
[closest_path_point,s_coordinate,path_point_yaw,first_path_point_index,...
    second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,fig_num);

% Print results to the workspace
fprintf(1,'Figure: %d\n',fig_num);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);


%% BASIC example 1.3 - breaks
fig_num = 10003;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

point = [1.5 1];
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 3 0];


% Snap the point onto the path
[closest_path_point,s_coordinate,path_point_yaw,first_path_point_index,...
    second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,fig_num);

% Print results to the workspace
fprintf(1,'Figure: %d\n',fig_num);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);


%% BASIC example - breaks and shows that it is on neither segement
fig_num = 10004;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

point = [0.9 1.4]; 
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 1.5 0.6; 3 0]; % Define an XY path


% Snap the point onto the path
[closest_path_point,s_coordinate,path_point_yaw,first_path_point_index,...
    second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,fig_num);

% Print results to the workspace
fprintf(1,'Figure: %d\n',fig_num);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);

%% BASIC example - works, but is on BOTH segments
fig_num = 100041;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

point = [1 0.5]; 
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 1.5 0.6; 3 0]; % Define an XY path


% Snap the point onto the path
[closest_path_point,s_coordinate,path_point_yaw,first_path_point_index,...
    second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,fig_num);

% Print results to the workspace
fprintf(1,'Figure: %d\n',fig_num);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);

%% BASIC example - negative s-coords (before path starts)
fig_num = 10005;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

point = [-0.5 0.2];
pathXY = [0 0; 1 0; 2 0; 2 1];


[closest_path_point,s_coordinate,path_point_yaw,first_path_point_index,...
    second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,fig_num);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);


%% BASIC example - positive s-coords (after path ends)
fig_num = 10006;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

point = [4 0.2];
pathXY = [0 0; 1 0; 2 0];


[closest_path_point,s_coordinate,path_point_yaw,first_path_point_index,...
    second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,fig_num);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

%% BASIC example - an example of percentage along segment greater than 100% even though "inside" path
fig_num = 10007;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

point = [0.8 1.3];
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 3 0];


[closest_path_point,s_coordinate,path_point_yaw,first_path_point_index,...
    second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,fig_num); %#ok<*ASGLU>
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fig_num, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);



%% ADVANCED example
fig_num = 20001;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;

% Create some paths:
path1 = [
    8.1797    8.6006
    8.4101   15.8892
    10.0230   25.5102
    10.4839   39.7959
    10.7143   50.8746
    10.9447   61.9534
    13.0184   73.6152
    16.7051   82.9446
    24.7696   89.3586
    35.3687   91.1079
    42.5115   91.3994
    52.8802   89.3586
    58.1797   85.2770
    60.9447   78.5714
    57.9493   72.1574
    57.2581   63.7026
    59.1014   58.1633
    63.0184   57.2886
    67.3963   56.9971
    69.9309   56.7055
    74.3088   56.1224
    78.6866   54.0816
    80.9908   51.4577
    82.3733   49.1254
    84.6774   40.6706
    84.6774   34.5481
    82.8341   28.7172
    80.0691   26.9679
    76.3825   25.2187
    69.2396   20.2624
    65.5530   18.2216
    60.9447   18.8047
    57.4885   22.0117
    50.5760   28.1341
    47.1198   30.7580
    43.8940   34.8397
    39.7465   37.7551
    37.6728   40.9621
    32.6037   42.7114
    30.0691   43.0029
    28.2258   43.2945
    26.3825   43.2945
    24.5392   44.4606
    20.8525   47.3761
    19.2396   49.7085
    16.7051   53.2070
    14.1705   58.1633
    13.0184   62.2449
    10.0230   70.1166
    8.6406   74.4898
    7.7189   79.7376
    6.5668   82.9446
    5.1843   86.7347
    4.2627   88.4840
    3.8018   89.0671];

figure(fig_num); 

plot(path1(:,1),path1(:,2),'r-o');
text(path1(1,1),path1(1,2),'Start');
pathXY = path1;

% Create a query
point = [75 45];
[closest_path_point,s_coordinate,path_point_yaw,first_path_point_index,...
    second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,fig_num);
fprintf(1,'Figure: %d, Closest point is: %.2f %.2f, S-coordinate is: %.2f \n',fig_num, closest_path_point(1,1),closest_path_point(1,2), s_coordinate);


%% ADVANCED example - more complicated path
fig_num = 20002;
fprintf(1,'Figure %.0f: basic demo 1\n',fig_num);
figure(fig_num); clf;


npoints = 20;
Ntests = 1;
seeds{Ntests} = 0;
for i_test = 1:Ntests
    seeds{i_test} = rng;    
    rng(seeds{i_test});
    rand_x = rand(npoints,1);
    rand_y = rand(npoints,1);
    
    pathXY = [cumsum(rand_x),cumsum(rand_y)];
    
    point = mean(pathXY,1);
    
    temp_fig_num = fig_num -1 + i_test;
    [closest_path_point,s_coordinate,path_point_yaw,first_path_point_index,...
        second_path_point_index,percent_along_length] = ...
        fcn_Path_snapPointOntoNearestPath(point, pathXY,temp_fig_num);
    fprintf(1,'Figure: %d, Closest point is: %.2f %.2f, S-coordinate is: %.2f \n',fig_num, closest_path_point(1,1),closest_path_point(1,2), s_coordinate);
end

% % Make sure plot opened up
% assert(isequal(get(gcf,'Number'),fig_num));

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

point = [0.8 1.3];
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 3 0];


[closest_path_point,s_coordinate,path_point_yaw,first_path_point_index,...
    second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,[]); %#ok<*ASGLU>

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

point = [0.8 1.3];
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 3 0];


[closest_path_point,s_coordinate,path_point_yaw,first_path_point_index,...
    second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,-1); %#ok<*ASGLU>

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

point = [0.8 1.3];
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 3 0];


Niterations = 100;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [closest_path_point,s_coordinate,path_point_yaw,first_path_point_index,...
        second_path_point_index,percent_along_length] = ...
        fcn_Path_snapPointOntoNearestPath(point, pathXY,[]); %#ok<*ASGLU>
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [closest_path_point,s_coordinate,path_point_yaw,first_path_point_index,...
        second_path_point_index,percent_along_length] = ...
        fcn_Path_snapPointOntoNearestPath(point, pathXY,-1); %#ok<*ASGLU>
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