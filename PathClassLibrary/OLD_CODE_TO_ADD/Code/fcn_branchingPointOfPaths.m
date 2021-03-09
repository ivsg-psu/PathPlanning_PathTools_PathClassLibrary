function [x_branch, y_branch] = fcn_branchingPointOfPaths(x1,y1,x2,y2)
% fcn_branchingPointOfPaths find the point where two paths branch 
L = length(x1);
x_branch = [];
y_branch = [];
d_closest = [];
x2_closest = [];
y2_closest = [];
for i = 1:L
    % distance from the ith point on path 1 (x1_i, y1_i) to all points on path 2
    d_i = ((x1(i)-x2).^2 + (y1(i)-y2).^2).^0.5;
    % find the index of the closest point on path 2 to the ith point on path 1
    idx = find(d_i == min(d_i));
    % there might be several same minimum distance points, pick the first one
    d_closest_i = d_i(idx(1));
    % add this point to the list of all minimum points
    x2_closest =[x2_closest; x2(idx(1))];
    y2_closest = [y2_closest; y2(idx(1))];
    d_closest =  [d_closest; d_closest_i];
end

idx_branch = 0;

for j = 1:10:length(d_closest)-10
    % if (x_j, y_j) is a branching point, then d will gradually increase
    if d_closest(j:j+10) == sort(d_closest(j:j+10)) & d_closest(j)>1.5
        if isempty(x_branch)
            x_branch = [x_branch; (x1(j) + x2_closest(j))/2];
            y_branch = [y_branch; (y1(j) + y2_closest(j))/2]; 
            idx_branch = idx_branch + 1;
        elseif (((x1(j) + x2_closest(j))/2 - x_branch(idx_branch)).^2 + ((y1(j) + y2_closest(j))/2 - y_branch).^2).^0.5 > 100
            x_branch = [x_branch; (x1(j) + x2_closest(j))/2];
            y_branch = [y_branch; (y1(j) + y2_closest(j))/2]; 
            idx_branch = idx_branch + 1;
        end
    end
end

