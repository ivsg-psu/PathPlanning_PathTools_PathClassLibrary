function [x_offset_left, y_offset_left, x_offset_right, y_offset_right] = fcn_calculatePathOffsets(x, y, s, s_coef)
% This function, given the average path and a specific sigma coefficient 
% value (say, s_coef = 1.7), calculates the left/right offset paths that 
% represent the bounds on 1.7 sigma
% Inputs: x ---- x coordinates of average path
%         y ---- y coordinates of average path
%         s ---- standard deviation of the probability distribution of user
%         input at each x-y point
%         s_coef ---- the number of standard deviations displayed for
%         shading 
% Outputs: x_offset_left --- the x coordinates of the left offset path
%          y_offset_left --- the y coordinates of the left offset path
%          x_offset_right --- the x coordinates of the right offset path
%          y_offset_right --- the y coordinates of the right offset path

x_offset_left = [];
y_offset_left = [];
x_offset_right = [];
y_offset_right = [];

step = 2;
for i = 1: length(x)-step
    % Create outer line buffer for every two points in the data
    x_buff = [x(i); x(i+step)];
    y_buff = [y(i); y(i+step)];
    buffer = polybuffer([x_buff, y_buff],'lines',s(i)*s_coef);
    buffer_pt = buffer.Vertices;
    x_offset_i = buffer_pt(:,1); % x coordinates of the buffer shape
    y_offset_i = buffer_pt(:,2); % y coordinates of the buffer shape
    
    % copy the first value of the array to the end to aviod discontinuity
    x_offset_i = [x_offset_i; x_offset_i(1)];
    y_offset_i = [y_offset_i; y_offset_i(1)];
    
    % MATLAB polybuffer function creates an offset shape from the line that
    % is determined by point i and point i+1. In this shape, there are many points near
    % the two ends of the line, yet the parallel line on each side is
    % only consisted of 2 points. That makes the distance between these 2
    % points significantly larger than that between other points
    
    d_all = []; % distance between point i and point i+1
    for j = 1:length(x_offset_i)-step
        % calculate distance 
        d = norm([x_offset_i(j),y_offset_i(j)]-[x_offset_i(j+step),y_offset_i(j+step)]);
        d_all = [d_all; d]; % append the calculated distance to array
    end
    % Find index of the points that are the starting point of the two parallel lines
    idx = find((d_all < 1.01*max(d_all)) & (d_all > 0.99* max(d_all)));
    
    dy = y(i+step) - y(i);
    dx = x(i+step) - x(i);
    if length(idx) == 4 % should always be true, just in case 
        idx1 = idx(1);
        idx2 = idx(3);
        if dy>0
            x_left = 0.5*(x_offset_i(idx1) + x_offset_i(idx1+step));
            y_left = 0.5*(y_offset_i(idx1) + y_offset_i(idx1+step));
            x_right = 0.5*(x_offset_i(idx2) + x_offset_i(idx2+step));
            y_right = 0.5*(y_offset_i(idx2) + y_offset_i(idx2+step));
        elseif dy<0
            x_right = 0.5*(x_offset_i(idx1) + x_offset_i(idx1+step));
            y_right = 0.5*(y_offset_i(idx1) + y_offset_i(idx1+step));
            x_left = 0.5*(x_offset_i(idx2) + x_offset_i(idx2+step));
            y_left = 0.5*(y_offset_i(idx2) + y_offset_i(idx2+step));
        else % horizontal line
            x_right = 0.5*(x_offset_i(idx1) + x_offset_i(idx1+step));
            x_left = 0.5*(x_offset_i(idx2) + x_offset_i(idx2+step));
            y_right = y_offset_i(idx1);
            y_left = y_offset_i(idx2);

        end
        x_offset_left = [x_offset_left; x_left];
        y_offset_left = [y_offset_left; y_left];
        x_offset_right = [x_offset_right; x_right];
        y_offset_right = [y_offset_right; y_right];
        
    end
end

[x_offset_left, y_offset_left] = fnc_eliminateIncontinuousData(x_offset_left, y_offset_left);
[x_offset_right, y_offset_right] = fnc_eliminateIncontinuousData(x_offset_right, y_offset_right);
end

