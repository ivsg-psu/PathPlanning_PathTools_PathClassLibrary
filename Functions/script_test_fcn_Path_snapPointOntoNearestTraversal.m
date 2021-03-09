% script_test_fcn_Path_snapPointOntoNearestTraversal.m
% This is a script to exercise the function: fcn_Path_snapPointOntoNearestTraversal.m
% This function was written on 2021_03_06 by Satya Prasad, szm888@psu.edu

% Revision history:


close all;

%% BASIC example 1
point = [0.5 0.2];
pathXY = [0 0; 1 0; 2 0; 2 1];
traversal = fcn_Path_convertPathToTraversalStructure(pathXY);

fignum = 111;
[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestTraversal(point, traversal,fignum);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fignum, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

%% BASIC example 1.2 - works
point = [1.4 1.3]; % Define the query point
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 1.5 0.6; 3 0]; % Define an XY path
traversal = fcn_Path_convertPathToTraversalStructure(pathXY);
fignum = 112; % Define the figure number

% Snap the point onto the path
[closest_path_point,s_coordinate,...
    first_path_point_index,second_path_point_index,...
    percent_along_length] = ...
    fcn_Path_snapPointOntoNearestTraversal(point, traversal,fignum);

% Print results to the workspace
fprintf(1,'Figure: %d\n',fignum);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);


%% BASIC example 1.3 - breaks
point = [1.5 1];
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 3 0];
traversal = fcn_Path_convertPathToTraversalStructure(pathXY);
fignum = 113;

% Snap the point onto the path
[closest_path_point,s_coordinate,...
    first_path_point_index,second_path_point_index,...
    percent_along_length] = ...
    fcn_Path_snapPointOntoNearestTraversal(point, traversal,fignum);

% Print results to the workspace
fprintf(1,'Figure: %d\n',fignum);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);


%% BASIC example 1.4 - breaks and shows that it is on neither segement
point = [0.9 1.4]; 
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 1.5 0.6; 3 0]; % Define an XY path
traversal = fcn_Path_convertPathToTraversalStructure(pathXY);
fignum = 114;

% Snap the point onto the path
[closest_path_point,s_coordinate,...
    first_path_point_index,second_path_point_index,...
    percent_along_length] = ...
    fcn_Path_snapPointOntoNearestTraversal(point, traversal,fignum);

% Print results to the workspace
fprintf(1,'Figure: %d\n',fignum);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);

%% BASIC example 1.5 - works, but is on BOTH segments
point = [1 0.5]; 
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 1.5 0.6; 3 0]; % Define an XY path
traversal = fcn_Path_convertPathToTraversalStructure(pathXY);
fignum = 115;

% Snap the point onto the path
[closest_path_point,s_coordinate,...
    first_path_point_index,second_path_point_index,...
    percent_along_length] = ...
    fcn_Path_snapPointOntoNearestTraversal(point, traversal,fignum);

% Print results to the workspace
fprintf(1,'Figure: %d\n',fignum);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);

%% BASIC example 2 - negative s-coords (before path starts)
point = [-0.5 0.2];
pathXY = [0 0; 1 0; 2 0; 2 1];
traversal = fcn_Path_convertPathToTraversalStructure(pathXY);
fignum = 222;
% [closest_path_point,s_coordinate] = ...
%     fcn_Path_snapPointOntoNearestTraversal(point, path,fignum);
% fprintf(1,'Figure: %d, Closest point is: %.2f %.2f, S-coordinate is: %.2f \n',fignum, closest_path_point(1,1),closest_path_point(1,2), s_coordinate);

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestTraversal(point, traversal,fignum);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fignum, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);


%% BASIC example 3 - positive s-coords (after path ends)
point = [4 0.2];
pathXY = [0 0; 1 0; 2 0];
traversal = fcn_Path_convertPathToTraversalStructure(pathXY);
fignum = 333;

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestTraversal(point, traversal,fignum);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fignum, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);

%% BASIC example 4 - an example of percentage along segment greater than 100% even though "inside" path
point = [0.8 1.3];
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 3 0];
traversal = fcn_Path_convertPathToTraversalStructure(pathXY);
fignum = 444;

[closest_path_point,s_coordinate,first_path_point_index,second_path_point_index,percent_along_length] = ...
    fcn_Path_snapPointOntoNearestTraversal(point, traversal,fignum);
fprintf(1,'Figure: %d,\n\t\t Closest point is: %.2f %.2f \n\t\t Matched to the path segment given by indices %d and %d, \n\t\t S-coordinate is: %.2f, \n\t\t percent_along_length is: %.2f\n',...
    fignum, closest_path_point(1,1),closest_path_point(1,2),...
    first_path_point_index,second_path_point_index, ...
    s_coordinate, percent_along_length);



%% ADVANCED example
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

figure(111); plot(path1(:,1),path1(:,2),'r-o');
text(path1(1,1),path1(1,2),'Start');
traversal = fcn_Path_convertPathToTraversalStructure(path1);
% Create a query
fignum = 2222;
point = [75 45];
[closest_path_point,s_coordinate] = ...
    fcn_Path_snapPointOntoNearestTraversal(point, traversal,fignum);
fprintf(1,'Figure: %d, Closest point is: %.2f %.2f, S-coordinate is: %.2f \n',fignum, closest_path_point(1,1),closest_path_point(1,2), s_coordinate);


%% ADVANCED example - more complicated path
npoints = 20;
Ntests = 10;
seeds{Ntests} = 0;
for i_test = 1:Ntests
    seeds{i_test} = rng;    
    rng(seeds{i_test});
    rand_x = rand(npoints,1);
    rand_y = rand(npoints,1);
    
    pathXY = [cumsum(rand_x),cumsum(rand_y)];
    traversal = fcn_Path_convertPathToTraversalStructure(pathXY);
    point = mean(pathXY,1);
    
    fignum = i_test;
    [closest_path_point,s_coordinate] = ...
        fcn_Path_snapPointOntoNearestTraversal(point, traversal,fignum);
    fprintf(1,'Figure: %d, Closest point is: %.2f %.2f, S-coordinate is: %.2f \n',fignum, closest_path_point(1,1),closest_path_point(1,2), s_coordinate);
end

