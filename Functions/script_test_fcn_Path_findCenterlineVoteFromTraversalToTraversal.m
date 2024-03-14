% script_test_fcn_Path_findCenterlineVoteFromTraversalToTraversal
% Tests the following:
%    [centerline_points_projected,unit_vectors_orthogonal] = ...
%     fcn_Path_findCenterlineVoteFromTraversalToTraversal(...
%     from_traversal,to_traversal,(flag_rounding_type),(search_radius),(fig_num))

% Revision history:
% 2023_09_04 by S. Brennan
% -- first write of the code

close all;


% Fill in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
for i_Path = 1:length(paths)
    traversal = fcn_Path_convertPathToTraversalStructure(paths{i_Path});
    data.traversal{i_Path} = traversal;
end

% Plot the results?
if 1==0
    fig_num = 12;
    fcn_Path_plotTraversalsYaw(data,fig_num);
    fig_num = 13;
    fcn_Path_plotTraversalsXY(data,fig_num);
end

%% Initialize the figure numbers
fig_num = 1000;


%% Basic demonstration 1 of fcn_Path_findCenterlineVoteFromTraversalToTraversal
% This function finds the center projected from one traversal toward
% another
from_path = [0 0; 1 1; 2 1; 3 4];
to_path   = from_path + ones(length(from_path(:,1)),1)*[0 1];
from_traversal =  fcn_Path_convertPathToTraversalStructure(from_path);
to_traversal =  fcn_Path_convertPathToTraversalStructure(to_path);
flag_rounding_type = 1;
search_radius = 10;
flag_project_full_distance = 0;
fig_num = 1;

[centerline_points_projected,unit_vectors_orthogonal] = ...
    fcn_Path_findCenterlineVoteFromTraversalToTraversal(...
    from_traversal,to_traversal,(flag_rounding_type),(search_radius),(flag_project_full_distance), (fig_num));

%% Basic demonstration 2 of fcn_Path_findCenterlineVoteFromTraversalToTraversal
% Show how, if the serach distance is too small, nothing is returned
from_path = [0 0; 1 1; 2 1; 3 4];
to_path   = from_path + ones(length(from_path(:,1)),1)*[0 1];
from_traversal =  fcn_Path_convertPathToTraversalStructure(from_path);
to_traversal =  fcn_Path_convertPathToTraversalStructure(to_path);
flag_rounding_type = 1;
search_radius = 0.1;
flag_project_full_distance = 0;
fig_num = 2;

[centerline_points_projected,unit_vectors_orthogonal] = ...
    fcn_Path_findCenterlineVoteFromTraversalToTraversal(...
    from_traversal,to_traversal,(flag_rounding_type),(search_radius),(flag_project_full_distance), (fig_num));

%% Basic demonstration 3 of fcn_Path_findCenterlineVoteFromTraversalToTraversal
% show that full projection returns mapping of 1 onto 2
from_path = [0 0; 1 1; 2 1; 3 4];
to_path   = from_path + ones(length(from_path(:,1)),1)*[0 1];
from_traversal =  fcn_Path_convertPathToTraversalStructure(from_path);
to_traversal =  fcn_Path_convertPathToTraversalStructure(to_path);
flag_rounding_type = 1;
search_radius = 10;
flag_project_full_distance = 1;
fig_num = 3;

[centerline_points_projected,unit_vectors_orthogonal] = ...
    fcn_Path_findCenterlineVoteFromTraversalToTraversal(...
    from_traversal,to_traversal,(flag_rounding_type),(search_radius),(flag_project_full_distance), (fig_num));

%% Call the center calculation function
from_traversal =  data.traversal{1};
to_traversal =  data.traversal{2};
flag_rounding_type = 1;
search_radius = 10;
flag_project_full_distance = 0;
fig_num = 4;

[centerline_points_projected,unit_vectors_orthogonal] = ...
    fcn_Path_findCenterlineVoteFromTraversalToTraversal(...
    from_traversal,to_traversal,(flag_rounding_type),(search_radius),(flag_project_full_distance), (fig_num));

