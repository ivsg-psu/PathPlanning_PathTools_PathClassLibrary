% script_test_fcn_Path_debugPrintStringToNCharacters.m
% This is a script to exercise the function: fcn_Path_debugPrintStringToNCharacters
% This function was written on 2021_01_24 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
%      2021_01_24:
%      -- first write of the code

close all;
clc;


%% BASIC example 1 - string is too long
test_string = 'This is a really really really long string';
fixed_length_string = fcn_Path_debugPrintStringToNCharacters(test_string,10);
fprintf(1,'The string: %s\nwas converted to: %s\n',test_string,fixed_length_string);

%% BASIC example 2 - string is too short
test_string = 'Short str';
fixed_length_string = fcn_Path_debugPrintStringToNCharacters(test_string,10);
fprintf(1,'The string: %s\nwas converted to: %s\n',test_string,fixed_length_string);

%% Advanced example
% This example shows why the function was written: to show information in a
% delimited format length

N_chars = 15;


% Create dummy data
Npoints = 10;
intersection_points = rand(Npoints,2);
s_coordinates_in_traversal_1 = rand(Npoints,1);
s_coordinates_in_traversal_2 = 1000*rand(Npoints,1);

% Print the header
header_1_str = sprintf('%s','Data ID');
fixed_header_1_str = fcn_Path_debugPrintStringToNCharacters(header_1_str,N_chars);
header_2_str = sprintf('%s','Location X');
fixed_header_2_str = fcn_Path_debugPrintStringToNCharacters(header_2_str,N_chars);
header_3_str = sprintf('%s','Location Y');
fixed_header_3_str = fcn_Path_debugPrintStringToNCharacters(header_3_str,N_chars);
header_4_str = sprintf('%s','s-coord 1');
fixed_header_4_str = fcn_Path_debugPrintStringToNCharacters(header_4_str,N_chars);
header_5_str = sprintf('%s','s-coord 2');
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
        results_2_str = sprintf('%.12f',intersection_points(ith_intersection,1));
        fixed_results_2_str = fcn_Path_debugPrintStringToNCharacters(results_2_str,N_chars);
        results_3_str = sprintf('%.12f',intersection_points(ith_intersection,2));
        fixed_results_3_str = fcn_Path_debugPrintStringToNCharacters(results_3_str,N_chars);
        results_4_str = sprintf('%.12f',s_coordinates_in_traversal_1(ith_intersection));
        fixed_results_4_str = fcn_Path_debugPrintStringToNCharacters(results_4_str,N_chars);
        results_5_str = sprintf('%.12f',s_coordinates_in_traversal_2(ith_intersection));
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

