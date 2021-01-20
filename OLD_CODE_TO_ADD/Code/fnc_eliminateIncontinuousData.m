function [x_new, y_new] = fnc_eliminateIncontinuousData(x, y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

angle = zeros(length(x)-2, 1);
x_new = x;
y_new = y;

for i = 2:length(x)-1
    P0 = [x(i-1), y(i-1)];
    P1 = [x(i), y(i)];
    P2 = [x(i+1), y(i+1)];
    n1 = (P2 - P0) / norm(P2 - P0);  % Normalized vectors
    n2 = (P1 - P0) / norm(P1 - P0);
    angle(i) = atan2d(norm(det([n2; n1])), dot(n1, n2));
end
idx = find(angle > 90);
n = length(idx);
idx = [idx; idx-1; idx+1];
x_new(idx) = [];
y_new(idx) = [];

while n
    angle = zeros(length(x_new)-2, 1);
    for i = 2:length(x_new)-1
        P0 = [x_new(i-1), y_new(i-1)];
        P1 = [x_new(i), y_new(i)];
        P2 = [x_new(i+1), y_new(i+1)];
        n1 = (P2 - P0) / norm(P2 - P0);  % Normalized vectors
        n2 = (P1 - P0) / norm(P1 - P0);
        angle(i) = atan2d(norm(det([n2; n1])), dot(n1, n2));
    end
    idx = find(angle > 90);
    n = length(idx);
    idx = [idx; idx-1; idx+1];
    %idx = [idx];
    x_new(idx) = [];
    y_new(idx) = [];
end
end