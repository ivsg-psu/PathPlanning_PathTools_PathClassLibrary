% script_test_fcn_Path_plotTraversalXYWithUpperLowerBands
% Tests fcn_Path_plotTraversalXYWithUpperLowerBands
       
% Revision history:
%      2022_01_03
%      -- first write of the code

close all


%% Test case 1: basic call for one trajectory
fig_num = 23333;
middle_path = [0 0; 1 1; 4 0; 5 0.5];
upper_path = [0 2; 1 2; 3.9 3; 4.9 1.5];
lower_path = [0 -1; 0.9 0; 4.1 -3; 5.1 0];
middle_traversal = fcn_Path_convertPathToTraversalStructure(middle_path);
upper_traversal = fcn_Path_convertPathToTraversalStructure(upper_path);
lower_traversal = fcn_Path_convertPathToTraversalStructure(lower_path);

fcn_Path_plotTraversalXYWithUpperLowerBands( middle_traversal, upper_traversal, lower_traversal, fig_num);
