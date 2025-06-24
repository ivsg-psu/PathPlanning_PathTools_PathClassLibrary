% script_test_fcn_Path_findTraversalStationSegment.m
% This is a script to exercise the function: 
% fcn_Path_findTraversalStationSegment.m
% This function was written on 2020_11_16 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

%      [traversal_trimmed,flag_outside_start, flag_outside_end] = ...
%      fcn_Path_findTraversalStationSegment(...
%      long_traversal, s_coord_start,s_coord_end, 
%      (fig_num))

% Revision history:
%     2021_01_09
%     -- updated name and types to take traversal inputs
%     -- added input checking
%     -- added flag_do_plots    

close all;
clear data;

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
for i_Path = 1:1  % length(paths)
    traversal = fcn_Path_convertPathToTraversalStructure(...
        paths_array{i_Path});
    data.traversal{i_Path} = traversal;
end

% Plot the results?
if 1==1
    fig_num = 12;
    fcn_Path_plotTraversalsYaw(data,fig_num);

    fig_num = 13;
    fcn_Path_plotTraversalsXY(data,fig_num);
end

    
%% BASIC example 1
s_coord_start = 10;
s_coord_end   = 100;
fignum = 111;
[traversal_segment1,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, fignum);
fprintf(1,...
    'Figure: %d, flag_outside_start is: %d, flag_outside_end is: %d \n',...
    fignum, flag_outside_start,flag_outside_end);
title('Normal query - Test case #1');

%% BASIC example 2 - another case
s_coord_start = 40;
s_coord_end   = 100;
fignum = 222;
[traversal_segment2,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, fignum);
fprintf(1,...
    'Figure: %d, flag_outside_start is: %d, flag_outside_end is: %d \n',...
    fignum, flag_outside_start,flag_outside_end);
title('Normal query - Test case #2');

%% BASIC example 3 - another case
s_coord_start = 70;
s_coord_end   = 80;
fignum = 333;
[traversal_segment3,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, fignum);
fprintf(1,...
    'Figure: %d, flag_outside_start is: %d, flag_outside_end is: %d \n',...
    fignum, flag_outside_start,flag_outside_end);
title('Normal query - Test case #3');

%% BASIC example 4 - degenerate within
s_coord_start = 70;
s_coord_end   = 70;
fignum = 444;
[traversal_segment4,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, fignum);
fprintf(1,...
    'Figure: %d, flag_outside_start is: %d, flag_outside_end is: %d \n',...
    fignum, flag_outside_start,flag_outside_end);
title('Query with exact same point');

%% BASIC example 5 - entire traversal within
s_coord_start = -5;
s_coord_end   = 7000;
fignum = 555;
[traversal_segment5,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, fignum);
fprintf(1,...
    'Figure: %d, flag_outside_start is: %d, flag_outside_end is: %d \n',...
    fignum, flag_outside_start,flag_outside_end);
title('Query where entire traversal within start and end of s-limits');

%% BASIC example 6 - end of traversal out
s_coord_start = 200;
s_coord_end   = 7000;
fignum = 666;
[traversal_segment6,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, fignum);
fprintf(1,...
    'Figure: %d, flag_outside_start is: %d, flag_outside_end is: %d \n',...
    fignum, flag_outside_start,flag_outside_end);
title('Query where end point outside of s-limits of traversal');

%% BASIC example 7 - start of traversal out
s_coord_start = -200;
s_coord_end   = 100;
fignum = 777;
[traversal_segment7,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, fignum);
fprintf(1,...
    'Figure: %d, flag_outside_start is: %d, flag_outside_end is: %d \n',...
    fignum, flag_outside_start,flag_outside_end);
title('Query where start point outside of s-limits of traversal');


%% BASIC example 9 - outside of start
s_coord_start = -200;
s_coord_end   = -100;
fignum = 999;
[traversal_segment9,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, fignum);
fprintf(1,...
    'Figure: %d, flag_outside_start is: %d, flag_outside_end is: %d \n',...
    fignum, flag_outside_start,flag_outside_end);
title('Query where start and end points are prior to start of traversal');

%% BASIC example 10 - outside of end
s_coord_start = 2000;
s_coord_end   = 3000;
fignum = 1111;
[traversal_segment10,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, fignum);
fprintf(1,...
    'Figure: %d, flag_outside_start is: %d, flag_outside_end is: %d \n',...
    fignum, flag_outside_start,flag_outside_end);
title('Query where start and end points are after the end of traversal');

%% Intentional warnings
fprintf(1,'\n\nTHE FOLLOWING WILL INTENTIONALLY PRODUCE ERRORS OR WARNINGS\n');

%% BASIC example 8 - incorrect call because start and end are out of order
s_coord_start = 200;
s_coord_end   = 100;
fignum = 888;
[traversal_segment8,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, fignum);
fprintf(1,...
    'Figure: %d, flag_outside_start is: %d, flag_outside_end is: %d \n',...
    fignum, flag_outside_start,flag_outside_end);
title('Query where start and end points are out of order (and auto-fixed)');

