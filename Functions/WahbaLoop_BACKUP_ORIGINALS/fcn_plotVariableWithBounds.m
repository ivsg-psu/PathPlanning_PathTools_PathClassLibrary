function fcn_plotVariableWithBounds(d,variableField,fig_number,fig_name,ylabel_string,title_string)
% Plots a generic variable with high/low bounds determined by sigma values
% added and subtracted to the median-filtered result of the data (2 sigma
% is used, and median filter of size 7 points). This is useful to
% illustrate which data would be cut off by a 2-sigma bound on the data.

% Revision history:
% 2019_10_05 - first write of code by S. Brennan, sbrennan@psu.edu
% 2019_10_20 - revised code to use dynamic fieldnames, rather than eval
% command. Also fixed yaw references throughout. And shut off LaTeX
% interpreter on title, etc.
% 2019_10_21 - fixed error if called with a value without sigmas. Added XY
% plotting capability. Added time to plot
% 2019_11_16 - added shaded color bands for sigmas. Fixed bug with data
% types on median calculation.
% 2019_11_26 - added check to prevent plotting errors if fields are
% missing... just throws warnings now rather than errors out of the
% function. This allows the merged data "sensor" to be ignored, if it
% doesn't exist, for example when plotting rawData's yaw_deg field when the
% merged data doesn't yet exist. Also made sure XY plot is square on axes.
% 2019_12_01 - Fixed bug with plotting XY data. This is not a field, so it
% was throwing an error during the input check process. Fixed it.

flag_plotMedian = 0;
flag_do_debug = 1;

%% Let the user know what we are doing
if flag_do_debug
    % Grab function name
    st = dbstack;
    namestr = st.name;
    
    % Show what we are doing
    fprintf(1,'\nWithin function: %s\n',namestr);
    fprintf(1,'\tPlotting data field: %s \n',variableField);
    fprintf(1,'\tPlotting to figure number: %d \n',fig_number);
    fprintf(1,'\tWith figure name: %s \n',fig_name);
    fprintf(1,'\tWith ylabel: %s \n',ylabel_string);
    fprintf(1,'\tWith title: %s \n',title_string);
    
end
%% Check that the variables exist - they might not
flag_do_plotting = 1;
% First, check if data is there, and grab it
if ~isfield(d,variableField) && ~isequal(variableField,'XYplot')
    if flag_do_debug
        fprintf(1,'\t WARNING: field %s not found in structure. Skipping plot.\n',variableField);
    end
    flag_do_plotting = 0;
end

if 1==flag_do_plotting
    
    %% Grab x-data for time plotting
    if ~isequal(variableField,'XYplot')
        % This is a variable versus time plot, time on x-axis
        
        % Decide if GPS or ROS time will be used - depends on sensors listed
        % Check to see if all data is from GPS time, or ROS time
        % Prefer GPS time as it's ahard standard...
        
        % Check to see if the GPS time is good for plotting
        flag_data_from_GPS_Time_is_Good = 1;
        if ~isfield(d,'GPS_Time')
            flag_data_from_GPS_Time_is_Good = 0;
        else
            temp_t = d.GPS_Time;
            if any(isnan(temp_t))
                flag_data_from_GPS_Time_is_Good = 0;
            end
        end

        if 0 == flag_data_from_GPS_Time_is_Good
            % Must use ROS time - recheck here that field exists and there is
            % no NaN. But at this point, throw errors if not.
            if ~isfield(d,'ROS_Time')
                error('Time plotting requested for variable that appears to have no valid time fields...???');
            else
                temp_t = d.ROS_Time;
                if any(isnan(temp_t))
                    error('Time plotting requested for variable that appears to have NaN values...???');
                end
            end
        end

        
        % Find the time with offset
        if flag_data_from_GPS_Time_is_Good
            t = d.GPS_Time - d.GPS_Time(1,1);
        else
            t = d.ROS_Time - d.ROS_Time(1,1);
        end % ends if statement on flag_dat_from_GPS_Time_is_Good
    end
    
    
    %% Set up the y data for time plotting
    if ~isequal(variableField,'XYplot')
        ydata = d.(variableField);
        sigma_name = cat(2,variableField,'_Sigma');
        if isfield(d,sigma_name)
            sigma = d.(sigma_name);
        else
            sigma = std(diff(ydata));
        end
        % x = (1:length(data(:,1)))';
        
        % Calculate median filtered values
        if isa(ydata,'single')||isa(ydata,'double')
            ydata_median = medfilt1(ydata,7,'truncate');
        else
            ydata_median = ydata;
        end
        % Add bounds to the position-based yawRate
        highest_expected_ydata = ydata_median + 2*sigma;
        lowest_expected_ydata  = ydata_median - 2*sigma;
    end
    
    %% Set up the figure
    if 1==flag_do_plotting
        try
            h_fig = figure(fig_number);
        catch
            disp('Debug here');
        end
        clf;
        set(h_fig,'Name',fig_name);
        hold on;
        h_title = title(title_string);
        set(h_title, 'Interpreter', 'none');
        grid minor;
    end
    
    %% Make the plots, depending on whether time plot or XY
    if ~isequal(variableField,'XYplot')
        
        % Insert plots
        plot(t,ydata,'b','LineWidth',1);
        if 1==flag_plotMedian
            plot(t,ydata_median,'c');
        end
        
        % Median filters only work on single or floating point classes
        if isa(ydata,'single')||isa(ydata,'double')
            fcn_plotVarianceBand(t,lowest_expected_ydata,highest_expected_ydata)
            if 1==flag_plotMedian
                legend('Data','median filtered','95% range');
            else
                legend('Data','95% range');
            end
        end
        
        h_ylabel = ylabel(ylabel_string); % set y label
        set(h_ylabel, 'Interpreter', 'none');
        
        % replce t with x if doing indices below...
        if flag_data_from_GPS_Time_is_Good==0
            xlabel('ROS Time [sec]');   % set  x label
        else
            xlabel('GPS Time [sec]');   % set  x label
        end
        
        
    else
        % This is an XY plot
        plot(d.xEast,d.yNorth,'LineWidth',1);
        xlabel('xEast [m]')  %set x label
        ylabel('yNorth [m]') %set y label
        axis square;        
    end
end %Ends the plotting if statement

%% End the function
if flag_do_debug
    % Show what we are doing
    fprintf(1,'Exiting function: %s\n',namestr);
end

return

%% Function to plot the band of variance
function fcn_plotVarianceBand(x,low_y,high_y)
% See: https://www.mathworks.com/matlabcentral/fileexchange/58262-shaded-area-error-bar-plot
% options.color_area = [128 193 219]./255;    % Blue theme

if 1==1  % This one looks best, but is memory intensive
    % Plotting the result with a patch object
    
    % We want to make a patch object, but it doesn't work with NaN values.
    % So we parse out the groups of data that contain NaN values
    [~,indices] = fcn_parseVectorByNaN(low_y);
    
    for i=1:length(indices)
        x_i = x(indices{i});
        high_y_i = high_y(indices{i});
        low_y_i = low_y(indices{i});
        
        % Now make the patch
        x_vector = [x_i', fliplr(x_i')];
        y_vector = [high_y_i',fliplr(low_y_i')];
        patch = fill(x_vector, y_vector,[128 193 219]./255);
        set(patch, 'edgecolor', 'none');
        set(patch, 'FaceAlpha', 0.5);
    end
else % Less memory intensive way is here - it just plots lines
    plot(x, low_y, 'Color',[128 193 219]./255, 'LineWidth', 1);
    hold on;
    plot(x, high_y,  'Color',[128 193 219]./255, 'LineWidth', 1);
    %legend('Data','median filtered','95% high','95% low');
end
return