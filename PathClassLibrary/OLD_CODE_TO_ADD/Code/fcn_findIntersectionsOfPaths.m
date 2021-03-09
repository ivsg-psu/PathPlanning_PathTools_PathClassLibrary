function [x_intersection, y_intersection] = fcn_findIntersectionsOfPaths(x1, y1, x2, y2)
% This function looks for joint intersections given two paths
% Requires the mapping toolbox
% Input: x1, y1, x2, y2 (each should be 1*n or 1*n array)
% Return: x_intersection, y_intersection
[x_intersection, y_intersection] = polyxpoly(x1,y1,x2,y2);
end
