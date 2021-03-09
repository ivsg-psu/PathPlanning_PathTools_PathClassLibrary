% script_test_fcn_Path_cleanPathFromForwardBackwardJogs
% Tests the function: fcn_Path_cleanPathFromForwardBackwardJogs

% Revision history
%     2020_01_09
%     -- first write of the code
%     -- need to add more test cases

close all;
clc;

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
for i_Path = 1:length(paths_array)
    traversal = fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
    data.traversal{i_Path} = traversal;
end

% Plot the results?
if 1==0
    fig_num = 12;
    fcn_Path_plotTraversalsYaw(data,fig_num);
    fig_num = 13;
    fcn_Path_plotTraversalsXY(data,fig_num);
end


%% Example 1: Basic call 
fig_num = 1;
path_with_jogs = [0 0; 1 1; 2 2.2; 3.3 3; 2.5 2.7; 3.5 3.6; 5 5];
clean_path = fcn_Path_cleanPathFromForwardBackwardJogs...
    (path_with_jogs,fig_num);

plot(path_with_jogs(:,1), path_with_jogs(:,2),'.-','Linewidth',2,'Markersize',25);
plot(clean_path(:,1), clean_path(:,2),'.-','Linewidth',2,'Markersize',25);
title('Original path with jogs and cleaned path')
xlabel('X [m]')
ylabel('Y [m]')

