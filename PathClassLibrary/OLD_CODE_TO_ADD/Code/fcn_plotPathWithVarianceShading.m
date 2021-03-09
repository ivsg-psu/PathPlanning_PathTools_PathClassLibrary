function fcn_plotPathWithVarianceShading(x, y,s, s_coef,c)
% This function plots the variance shading corresponding to an average path
% Inputs: x ---- x coordinates of path
%         y ---- y coordinates of path
%         s ---- standard deviation of the probability distribution of user
%         input at each x-y point
%         s_coef ---- the number of standard deviations displayed for
%         shading 
% Output: None
if nargin < 4
    s_coef = 1;
    s = zeros(length(x),1)+3;
    for i = 1: length(x)-1
        buffer = polybuffer([x(i:i+1),y(i:i+1)],'lines',s(i)*s_coef);
        buffer_pt = buffer.Vertices;
        x_offset_i = buffer_pt(:,1);
        y_offset_i = buffer_pt(:,2);
        x_offset_i = [x_offset_i; x_offset_i(1)];
        y_offset_i = [y_offset_i; y_offset_i(1)];
        hold on

        p = fill(x_offset_i, y_offset_i, c);
        set(p,'facealpha', 0.1)    
        p.EdgeColor = 'none'; 
    end 
else
    for i = 1: length(x)-1
        buffer = polybuffer([x(i:i+1),y(i:i+1)],'lines',s(i)*s_coef);
        buffer_pt = buffer.Vertices;
        x_offset_i = buffer_pt(:,1);
        y_offset_i = buffer_pt(:,2);
        x_offset_i = [x_offset_i; x_offset_i(1)];
        y_offset_i = [y_offset_i; y_offset_i(1)];
        hold on
        if c == 'r'
            p = fill(x_offset_i, y_offset_i, 'red');
            p.FaceColor = [1 0.8 0.8];
            set(p,'facealpha', 0.4)
            p.EdgeColor = 'none'; 
        elseif c == 'g'
            p = fill(x_offset_i, y_offset_i, 'green');
            set(p,'facealpha', 0.4)
            p.FaceColor = [0.8 1 0.8];      
            p.EdgeColor = 'none'; 
        end
   
    end
    plot(x,y, c, 'LineWidth', 1)
end

end

