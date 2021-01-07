% script_test_fcn_Path_findPathSXYSegment.m
% This is a script to exercise the function: 
% fcn_Path_findPathSXYSegment.m
% This function was written on 2020_11_16 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

close all;
clear data;

% Fill in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
for i_Path = 1:1  % length(paths)
    traversal = fcn_Path_convertXYtoTraversalStructure(paths{i_Path}(:,1),paths{i_Path}(:,2));
    data.traversal{i_Path} = traversal;
end

% Plot the results?
if 1==1
    % fig_num = 12;
    % fcn_Path_plotPathYaw(data,fig_num);

    fig_num = 13;
    fcn_Path_plotPathXY(data,fig_num);
end
path1 = paths{1};
text(path1(1,1),path1(1,2),'Start');

pathSXY1 = fcn_Path_convertXYtoSXY(path1(:,1),path1(:,2));
    
%% BASIC example 1
s_coord_start = 10;
s_coord_end   = 100;
fignum = 111;
[pathSXY_segment1,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findPathSXYSegment(pathSXY1, s_coord_start,s_coord_end, fignum);
fprintf(1,'Figure: %d, flag_outside_start is: %d, flag_outside_end is: %d \n',fignum, flag_outside_start,flag_outside_end);
title('Normal query - Test case #1');

%% BASIC example 2 - another case
s_coord_start = 40;
s_coord_end   = 100;
fignum = 222;
[pathSXY_segment2,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findPathSXYSegment(pathSXY1, s_coord_start,s_coord_end, fignum);
fprintf(1,'Figure: %d, flag_outside_start is: %d, flag_outside_end is: %d \n',fignum, flag_outside_start,flag_outside_end);
title('Normal query - Test case #2');

%% BASIC example 3 - another case
s_coord_start = 70;
s_coord_end   = 80;
fignum = 333;
[pathSXY_segment3,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findPathSXYSegment(pathSXY1, s_coord_start,s_coord_end, fignum);
fprintf(1,'Figure: %d, flag_outside_start is: %d, flag_outside_end is: %d \n',fignum, flag_outside_start,flag_outside_end);
title('Normal query - Test case #3');

%% BASIC example 4 - degenerate within
s_coord_start = 70;
s_coord_end   = 70;
fignum = 444;
[pathSXY_segment4,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findPathSXYSegment(pathSXY1, s_coord_start,s_coord_end, fignum);
fprintf(1,'Figure: %d, flag_outside_start is: %d, flag_outside_end is: %d \n',fignum, flag_outside_start,flag_outside_end);
title('Query with exact same point');

%% BASIC example 5 - entire path within
s_coord_start = -5;
s_coord_end   = 7000;
fignum = 555;
[pathSXY_segment5,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findPathSXYSegment(pathSXY1, s_coord_start,s_coord_end, fignum);
fprintf(1,'Figure: %d, flag_outside_start is: %d, flag_outside_end is: %d \n',fignum, flag_outside_start,flag_outside_end);
title('Query where entire path within start and end of s-limits');

%% BASIC example 6 - end of path out
s_coord_start = 200;
s_coord_end   = 7000;
fignum = 666;
[pathSXY_segment6,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findPathSXYSegment(pathSXY1, s_coord_start,s_coord_end, fignum);
fprintf(1,'Figure: %d, flag_outside_start is: %d, flag_outside_end is: %d \n',fignum, flag_outside_start,flag_outside_end);
title('Query where end point outside of s-limits of path');

%% BASIC example 7 - start of path out
s_coord_start = -200;
s_coord_end   = 100;
fignum = 777;
[pathSXY_segment7,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findPathSXYSegment(pathSXY1, s_coord_start,s_coord_end, fignum);
fprintf(1,'Figure: %d, flag_outside_start is: %d, flag_outside_end is: %d \n',fignum, flag_outside_start,flag_outside_end);
title('Query where start point outside of s-limits of path');


%% BASIC example 9 - outside of start
s_coord_start = -200;
s_coord_end   = -100;
fignum = 999;
[pathSXY_segment9,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findPathSXYSegment(pathSXY1, s_coord_start,s_coord_end, fignum);
fprintf(1,'Figure: %d, flag_outside_start is: %d, flag_outside_end is: %d \n',fignum, flag_outside_start,flag_outside_end);
title('Query where start and end points are prior to start of path');

%% BASIC example 10 - outside of end
s_coord_start = 2000;
s_coord_end   = 3000;
fignum = 1111;
[pathSXY_segment10,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findPathSXYSegment(pathSXY1, s_coord_start,s_coord_end, fignum);
fprintf(1,'Figure: %d, flag_outside_start is: %d, flag_outside_end is: %d \n',fignum, flag_outside_start,flag_outside_end);
title('Query where start and end points are after the end of path');

%% Intentional warnings
fprintf(1,'\n\nTHE FOLLOWING WILL INTENTIONALLY PRODUCE ERRORS OR WARNINGS\n');

%% BASIC example 8 - incorrect call because start and end are out of order
s_coord_start = 200;
s_coord_end   = 100;
fignum = 888;
[pathSXY_segment8,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findPathSXYSegment(pathSXY1, s_coord_start,s_coord_end, fignum);
fprintf(1,'Figure: %d, flag_outside_start is: %d, flag_outside_end is: %d \n',fignum, flag_outside_start,flag_outside_end);
title('Query where start and end points are out of order (and auto-fixed)');

