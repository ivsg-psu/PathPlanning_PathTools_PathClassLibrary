function [PVIP_x, PVIP_y, cluster_idx] = fcn_obtainPVIPs(x1, y1,s1, x2, y2, s2, k)
if nargin < 7
    k = 1;
end
s_coef = (0.5:0.5:5)';
PVIP_x = [];
PVIP_y = [];
[x_branch, y_branch] = fcn_branchingPointOfPaths(x1, y1, x2, y2);
for i = 1:length(s_coef)
    [x1_offset_left, y1_offset_left, ~, ~] = fcn_calculatePathOffsets(x1, y1, s1, s_coef(i));
    [~, ~, x2_offset_right, y2_offset_right] = fcn_calculatePathOffsets(x2, y2, s2, s_coef(i));
    
    [x_offset_int, y_offset_int] = fcn_findIntersectionsOfPaths(x1_offset_left, y1_offset_left, x2_offset_right, y2_offset_right);
    
    for j = 1:length(x_branch)
        d_PVIP_Branch = ((x_offset_int - x_branch(j)).^2 + (y_offset_int - y_branch(j)).^2).^0.5;
        PVIP_idx  =  find(d_PVIP_Branch < 30); % not to far away from the branching point  
        %if (x_offset_int(PVIP_idx) < 1155) | (x_offset_int(PVIP_idx) > 1156)
            PVIP_x = [PVIP_x; x_offset_int(PVIP_idx)];
            PVIP_y = [PVIP_y; y_offset_int(PVIP_idx)];
        %end
        
    end
end
cluster_idx = kmeans([PVIP_x, PVIP_y],k);

%{
figure
fcn_plotPathWithVarianceShading(x1, y1, s1, 1, 'r')
fcn_plotPathWithVarianceShading(x2, y2, s2, 1, 'g')
hold on

for i = 1:k
    plot(PVIP_x(cluster_idx == i), PVIP_y(cluster_idx == i), 'o');
end

xlabel("x [m]")
ylabel("y [m]")
hold off
%} 
end

