function fcn_plotMergedSources(d,fig_number,fig_name,field_name)

%% Set up the figure
h_fig = figure(fig_number);
clf;
set(h_fig,'Name',fig_name);
hold on;

%% Create time vectors
t = d.Clocks.targetTimeVector_GPS{d.(field_name).centiSeconds};
t = t - t(1);

%% Fill in data
C = d.(field_name).Center;
U = d.(field_name).Upper;
L = d.(field_name).Lower;
sigma = d.(field_name).Sigma;

%% Create plots
figure(fig_number);
clf;
hold on;
h_plot = plot(t,C,'Linewidth',1);
plot(t,U,'Linewidth',1,'Color',[0 0 0]);
plot(t,L,'Linewidth',1,'Color',[0 0 0]);

% Plot standard deviation range
plotColor = get(h_plot,'Color');
newColor = plotColor * 0.5 + [1 1 1]*0.5;
data_median = medfilt1(C,7,'truncate');
highsigma = data_median + 2*sigma;
lowsigma  = data_median - 2*sigma;
fcn_plotVarianceBand(t,lowsigma,highsigma,newColor);


xlabel('Time [sec]')   % set  x label
ylabel(field_name); 
grid minor;

title(cat(2,'Plot of ',field_name,' from merged data'));
legend('Central estimate','Upper estimate','Lower estimate','2-sigma range');
end

%% Function to plot the band of variance

function fcn_plotVarianceBand(x,low_y,high_y,color)
% See: https://www.mathworks.com/matlabcentral/fileexchange/58262-shaded-area-error-bar-plot
% options.color_area = [128 193 219]./255;    % Blue theme

if 1==1  % This one looks best, but is memory intensive
    % Plotting the result
    x_vector = [x', fliplr(x')];
    y_vector = [high_y',fliplr(low_y')];
    patch = fill(x_vector, y_vector,color);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', 0.5);
else % Less memory intensive way is here - it just plots lines
    plot(x, low_y, 'r', 'LineWidth', 1);
    hold on;
    plot(x, high_y, 'r', 'LineWidth', 1);
    %legend('Data','median filtered','95% high','95% low');
end
end