% script_test_fcn_Path_findPathOrthogonalVectors
% This is a script to exercise the function: fcn_Path_findPathOrthogonalVectors
% This function was written on 2023_08_27 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history
%     2023_08_27
%     -- first write of the code

close all;

%% BASIC example 1 - simple horizontal line, default flag (1), no figure
path = [0 0; 4 0]; 

% Call the function
[normal_unit_vectors_at_midpoints, normal_unit_vectors_at_joints] = ...
    fcn_Path_findPathOrthogonalVectors(path);

% Make sure function worked
assert(isequal(round(normal_unit_vectors_at_midpoints,4),[0 1]));
assert(isequal(round(normal_unit_vectors_at_joints,4),[0 1; 0 1]));

%% BASIC example 1 - simple horizontal line, default flag (1), with figure
fig_num = 2;
flag_rounding_type = 1; % Define the rounding type
path = [0 0; 4 0];  

% Call the function
[normal_unit_vectors_at_midpoints, normal_unit_vectors_at_joints] = ...
    fcn_Path_findPathOrthogonalVectors(path,flag_rounding_type, fig_num);

% Make sure function worked
assert(isequal(round(normal_unit_vectors_at_midpoints,4),[0 1]));
assert(isequal(round(normal_unit_vectors_at_joints,4),[0 1; 0 1]));


%% BASIC example 3 - right-angled line segment - flag tests
fig_num = 3;
figure(fig_num);
clf;

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



%% Real path examples
flag_rounding_type = 1; % Define the rounding type

% Fill in some dummy data
paths = fcn_Path_fillSamplePaths;


for i_Path = 1:length(paths)
    fig_num = 40+i_Path;  % Define the figure
    figure(fig_num);
    clf;
    
    % Call the function
    [normal_unit_vectors_at_midpoints, normal_unit_vectors_at_joints] = ...
        fcn_Path_findPathOrthogonalVectors(paths{i_Path},flag_rounding_type, fig_num);
end


