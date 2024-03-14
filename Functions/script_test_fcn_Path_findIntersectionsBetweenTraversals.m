% script_test_fcn_Path_findIntersectionsBetweenTraversals.m
% This is a script to exercise the function: fcn_Path_findIntersectionsBetweenTraversals
% This function was written on 2020_11_14 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
%      2021_01_23:
%      -- first write of the code

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
fig_num = 1;

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

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);

%% BASIC example 2 - perpendicular lines, intersection in middle
fig_num = 2;

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

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);


%% BASIC example 5 - two crossings of 1 onto 2
fig_num = 5;

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

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);


%% BASIC example 6 - intersection at vertex of one traversal
fig_num = 6;

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

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);


%% BASIC example 7 - intersection at vertex of both traversals
fig_num = 7;

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

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);


%% BASIC example 8 - intersection at both ends of one traversal
fig_num = 8;

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

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);


%% BASIC example 9 - intersection at both ends of one traversal with two vertices
fig_num = 9;

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

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);

%% BASIC example 10 - Loop of traversal 1 over same point twice
fig_num = 10;

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

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);


%% BASIC example 11 - Loop of traversal 2 over same point twice
fig_num = 11;

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

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);


%% BASIC example 12 - Traversal 2 is on top of traversal 1 for an area
fig_num = 12;

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

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);



%% BASIC example 13 - Traversal 1 is on top of traversal 2 for an area
fig_num = 13;

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

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);


%% BASIC example 14 - Traversal 2 is on top of traversal 1 for two areas
fig_num = 14;

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

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);

%% BASIC example 15 - Traversal 1 is on top of traversal 2 for two areas
fig_num = 15;

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

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);

%% BASIC example 16 - Traversal 1 is same as Traversal 2 for many segments
fig_num = 16;

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

print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2);






%% Functions start here


function print_results(intersection_points,s_coordinates_in_traversal_1,s_coordinates_in_traversal_2)
N_chars = 25;


% Print the header
header_1_str = sprintf('%s','Intersection ID');
fixed_header_1_str = fcn_Path_debugPrintStringToNCharacters(header_1_str,N_chars);
header_2_str = sprintf('%s','Location X');
fixed_header_2_str = fcn_Path_debugPrintStringToNCharacters(header_2_str,N_chars);
header_3_str = sprintf('%s','Location Y');
fixed_header_3_str = fcn_Path_debugPrintStringToNCharacters(header_3_str,N_chars);
header_4_str = sprintf('%s','s-coord in traversal 1');
fixed_header_4_str = fcn_Path_debugPrintStringToNCharacters(header_4_str,N_chars);
header_5_str = sprintf('%s','s-coord in traversal 2');
fixed_header_5_str = fcn_Path_debugPrintStringToNCharacters(header_5_str,N_chars);

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
        fixed_results_1_str = fcn_Path_debugPrintStringToNCharacters(results_1_str,N_chars);
        results_2_str = sprintf('%.2f',intersection_points(ith_intersection,1));
        fixed_results_2_str = fcn_Path_debugPrintStringToNCharacters(results_2_str,N_chars);
        results_3_str = sprintf('%.2f',intersection_points(ith_intersection,2));
        fixed_results_3_str = fcn_Path_debugPrintStringToNCharacters(results_3_str,N_chars);
        results_4_str = sprintf('%.2f',s_coordinates_in_traversal_1(ith_intersection));
        fixed_results_4_str = fcn_Path_debugPrintStringToNCharacters(results_4_str,N_chars);
        results_5_str = sprintf('%.2f',s_coordinates_in_traversal_2(ith_intersection));
        fixed_results_5_str = fcn_Path_debugPrintStringToNCharacters(results_5_str,N_chars);
        
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

