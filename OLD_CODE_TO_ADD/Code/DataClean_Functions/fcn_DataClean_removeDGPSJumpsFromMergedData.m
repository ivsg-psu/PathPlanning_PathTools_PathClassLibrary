function mergedData = fcn_DataClean_removeDGPSJumpsFromMergedData(mergedData,rawData)

% Revision history:
% 2019_12_01 - first write of the function by sbrennan@psu.edu

flag_do_debug = 1;


data_to_fit = rawData.GPS_Hemisphere;
xEast_pred_increments  = mergedData.MergedGPS.velMagnitude*0.05 .* cos(mergedData.MergedGPS.Yaw_deg*pi/180);
yNorth_pred_increments = mergedData.MergedGPS.velMagnitude*0.05 .* sin(mergedData.MergedGPS.Yaw_deg*pi/180);


%% Let the user know what we are doing
if flag_do_debug
    % Grab function name
    st = dbstack;
    namestr = st.name;

    % Show what we are doing
    fprintf(1,'\nWithin function: %s\n',namestr);
    fprintf(1,'Correct differntial jumps in xEast and yNorth, in merged data.\n');   
    fprintf(1,'Length of data vector going in: %d\n',length(mergedData.MergedGPS.xEast));
end


%% Find valid intervals where DGPS is active
pairings = fcn_DataClean_findStartEndPairsWhereDGPSDrops(data_to_fit.DGPS_is_active);
    
%% Using the valid intervals above, start fixing the data

if ~isempty(pairings)
    for i_pairing = 1:length(pairings(:,1))

        % ONLY NEED THE FOLLOWING IF DOING TIME PLOTS, OR ANALYSIS TO FIND TIME
        % LOCATIONS OF DROP-OUTS
        %     startTime = data_to_fit.GPS_Time(pairings(i_pairing,1)) - data_to_fit.GPS_Time(i_pairing,1);
        %     endTime = data_to_fit.GPS_Time(pairings(i_pairing,2))  - data_to_fit.GPS_Time(i_pairing,1);
        %
        %     plottingFlags.TimeZoomPoint = [startTime endTime]; % Strange jump in xEast data
        %
        %
        %     % First, grab the data
        %     t_min = min(plottingFlags.TimeZoomPoint)+data_to_fit.GPS_Time(1,1);
        %     t_max = max(plottingFlags.TimeZoomPoint)+data_to_fit.GPS_Time(1,1);


        indices_of_interest = pairings(i_pairing,1):pairings(i_pairing,2);
        DGPS_is_active = data_to_fit.DGPS_is_active(indices_of_interest);

        % Fix xEast data
        xEast = data_to_fit.xEast(indices_of_interest);   
        xEast_increment_pred = xEast_pred_increments(indices_of_interest);
        [xEast_clean,~] = fcn_DataClean_medianFilterViaIncrementTemplate(xEast,xEast_increment_pred,DGPS_is_active);
        mergedData.MergedGPS.xEast(indices_of_interest) = xEast_clean;

        % Fix yNorth data
        yNorth = data_to_fit.yNorth(indices_of_interest);   
        yNorth_increment_pred = yNorth_pred_increments(indices_of_interest);
        [yNorth_clean,~] = fcn_DataClean_medianFilterViaIncrementTemplate(yNorth,yNorth_increment_pred,DGPS_is_active);
        mergedData.MergedGPS.yNorth(indices_of_interest) = yNorth_clean;    
    end
end 
if flag_do_debug
    % Show what we are doing
    fprintf(1,'Exiting function: %s\n',namestr);
    fprintf(1,'Length of data vector going out: %d\n',length(mergedData.MergedGPS.xEast));
end

return

function pairings = fcn_DataClean_findStartEndPairsWhereDGPSDrops(DGPS_is_active)
%%
% Find all the locations where DGPS shuts off and then on. Require that
% DGPS needs to be on at least 1 second (20 samples) to be "on". Can set
% this number higher if needed
n_samples_on = 100;

% Set the start and end points of status to zero, to indicate DGPS is
% gained/lost at these points
temp_DGPS_status = [0; DGPS_is_active]; % Artificially shut off at first index
temp_DGPS_status(end,1) = 0; % Artificially shut off at last index

diff_DGPS = diff(temp_DGPS_status);
indices_DGPS_found = find(diff_DGPS==1);
indices_DGPS_lost = find(diff_DGPS==-1);

% Initialize the good_diff vector
good_DGPS_status = temp_DGPS_status;

% To enforce the requirement that the DGPS be on at least n_samples_on
% before a loss valid, we loop through all the locations where DGPS is
% lost, and set all the n_samples prior to this to zero.
for i_loss=1:length(indices_DGPS_lost)
    current_index = indices_DGPS_lost(i_loss); % Grab the current location where DGPS was lost
    start_index = max(1,current_index-n_samples_on);
    good_DGPS_status(start_index:current_index,1) = 0;    
end  

% To enforce the requirement that the DGPS be on at least n_samples_on
% after a loss return, we loop through all the locations where DGPS is
% gained, and set all the n_samples after to this to zero.
for i_found=1:length(indices_DGPS_found)
    current_index = indices_DGPS_found(i_found); % Grab the current location where DGPS was lost
    end_index = min(length(diff_DGPS),current_index+n_samples_on);
    good_DGPS_status(current_index:end_index,1) = 0;    
end  

% Now the vector good_diff_GPS only goes off/on at index intervals that
% satisfy the requirement
diff_DGPS = diff(good_DGPS_status);
indices_DGPS_found = find(diff_DGPS==1);
indices_DGPS_lost = find(diff_DGPS==-1);

% Form pairings:
pairings = [];
for i_lost = 1:length(indices_DGPS_lost)
    
    i_found = find(indices_DGPS_found>indices_DGPS_lost(i_lost),1);
    if ~isempty(i_found) % found a pair
        pairings = [pairings; [indices_DGPS_lost(i_lost) indices_DGPS_found(i_found)]]; %#ok<AGROW>
    end
end

if 1==0
    % Check results
    for i_pair = 1:length(pairings(:,1))
        index_start = pairings(i_pair,1);
        index_end = pairings(i_pair,2);
        ydata = data_to_fit.DGPS_is_active(index_start:index_end,1);
        xdata = (index_start:index_end)'*0.05;
        figure(36363);
        clf;
        plot(xdata,ydata,'k');
        xlabel('Time (sec)');
        ylabel('DGPS status');
        xlim([index_start index_end]*0.05);
        pause;
    end
end