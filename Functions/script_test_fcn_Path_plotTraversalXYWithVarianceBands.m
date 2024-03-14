% script_test_fcn_Path_plotTraversalXYWithVarianceBands
% Tests fcn_Path_plotTraversalXYWithVarianceBands
       
% Revision history:
%      2021_01_05
%      -- first write of the code
%      2021_01_07
%      -- fixed naming convention on functions to reflect change from path to
%      traversal notation
%      2021_01_08
%      -- more fixes

close all




%% Test case 1: basic call for one trajectory

% Fill in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths{1});

fcn_Path_plotTraversalXYWithVarianceBands(reference_traversal);


%% Test case 2: advanced call for one trajectory - specify figure
fig_num = 22;


% Fill in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths{1});

std_deviation = [];
fcn_Path_plotTraversalXYWithVarianceBands(reference_traversal,...
    std_deviation,fig_num);

%% Test case 3: advanced call for one trajectory - specify std_deviation
fig_num = 31;


% Fill in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths{1});

std_deviation = 1;
fcn_Path_plotTraversalXYWithVarianceBands(reference_traversal,...
    std_deviation,fig_num);
title(sprintf('Standard deviation: %.0d meters',std_deviation));

fig_num = 32;
std_deviation = 2;
fcn_Path_plotTraversalXYWithVarianceBands(reference_traversal,...
    std_deviation,fig_num);
title(sprintf('Standard deviation: %.0d meters',std_deviation));

fig_num = 35;
std_deviation = 5;
fcn_Path_plotTraversalXYWithVarianceBands(reference_traversal,...
    std_deviation,fig_num);
title(sprintf('Standard deviation: %.0d meters',std_deviation));

%% Test case 4: advanced call for multiple trajectories
fig_num = 4;


% Fill in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths{1});

std_deviation = 2;
for i_Path = 1:length(paths)
    reference_traversal = fcn_Path_convertPathToTraversalStructure(paths{i_Path});
    fcn_Path_plotTraversalXYWithVarianceBands(reference_traversal,...
        std_deviation,fig_num);
end