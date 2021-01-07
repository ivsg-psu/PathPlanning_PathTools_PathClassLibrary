function cleanAndTimeAlignedData = fcn_alignToGPSTimeAllData(cleanData)

% This function is used to align in time all the data such that their ROS
% and GPS time stamps are uniform in time, and they all have the same start and
% end times by being trimmed to same points in time. 
%
% Author: Sean Brennan
% First written: 2019_10_06
% Last modify Date: 2019_11_27
%
% Updates:
%  2019_10_06 - first writing by Dr. B
%  2019_10_11 - finished writing
%  2019_10_17 - fixed bug for interpolation of integeters producing a float
%  (navMode, for example)
%  2019_11_18 - fixed bug with extrapolation value being 0, should have
%  been nan. Otherwise, nan fields begin to be filled with 0's.
%  2019_11_25 - fixed bug with sigma values being filled with 0 during
%  extrapolation (above bug fix didn't work - changed to extrapolation now)
%  2019_11_26 - fixed bug with interp1 being passed int64. Used typecast
%  function to convert data to double, if integers are encountered. Also
%  forced the GPS time to be created, in case that the ROS time is the only
%  time vector (for example, with the encoder data).
%
% To do items:
% (as of 2019_10_17) - clean up field checks to do this all at once (see
% fcn_loadSigmaValuesFromRawData for an example)

flag_showMotivation = 0;
flag_do_debug = 1;

fields_to_interpolate_to_nearest_neighbor = [...
    {'navMode'},...
    {'DGPS_is_active'},...
    {'IMUStatus'},...
    {'placeholder'},...
    ];

%% The following code shows why the time alignment is needed, e.g. the motivation
if 1==flag_showMotivation
    % The only source we have (currently) of the GPS to ROS offset is from the
    % IMU_Novatel data. We'll need to fix this in the future.
    
    % For debugging: The following plot shows that the IMU is lock-on at 0.01
    % seconds, but occasionally (surprisingly) misses samples
    figure; histogram(diff(cleanData.IMU_Novatel.GPS_Time),10000);
    
    % However, the following examples show that the raw time data measured in
    % ROS does not actually count up correctly:
    cleanData.IMU_ADIS.ROS_Time(1:10)-cleanData.IMU_ADIS.ROS_Time(1) %#ok<NOPRT>
    cleanData.IMU_Novatel.ROS_Time(1:10)-cleanData.IMU_Novatel.ROS_Time(1) %#ok<NOPRT>
    
    % The histogram of the ROS delta_times show that the vast majority of
    % updates occur at around 50 to 75 ms update rates. Thus, on average, 5 to
    % 7 samples go by per each time step, with each sample not labeled
    % correctly.
    figure(344576); histogram(diff(cleanData.IMU_Novatel.ROS_Time),10000);
    ylim([0 100]);
    % However, the mean is pretty close, at around 0.011 ms. For the ADIS, the
    % data come in at 16 ms. But the mean is lock on at 0.0100 seconds.
    figure(344555); histogram(diff(cleanData.IMU_ADIS.ROS_Time),10000);
    ylim([0 100]);
    %%%%%%%% RESULTS %%%%%%%%%%%%%%
    %mean(diff(cleanData.Encoder_RearWheels.ROS_Time))  % Gives 0.009998581109336 seconds
    mean(diff(cleanData.GPS_Hemisphere.ROS_Time))      % Gives 0.050001513198616 seconds
    mean(diff(cleanData.GPS_Novatel.ROS_Time))         % Gives 0.050000566338432 seconds
    mean(diff(cleanData.GPS_Garmin.ROS_Time))          % Gives 0.040000738800513 seconds
    mean(diff(cleanData.IMU_ADIS.ROS_Time))            % Gives 0.010039343452559 seconds
    mean(diff(cleanData.IMU_Novatel.ROS_Time))         % Gives 0.011675284143736 seconds
    %mean(diff(cleanData.Input_Steering.ROS_Time))      % Gives 0.009998571063095 seconds
    
    % Clearly, the data has to be re-decimated. But to what time base? While it
    % might seem good to use the ROS time, only the GPS time is absolute. Since
    % only a few sensors - GPS - have this, we need to use them to create an
    % absolute decimation for the other sensors.
    
    % Step1: Now we create a uniform sampling in the GPS time sensors. TO do this,
    % compare durations of ROS and GPS time, to make sure they agree.
    
    durationsGPS = cleanData.GPS_Novatel.GPS_Time - cleanData.GPS_Novatel.GPS_Time(1);
    durationsROS = cleanData.GPS_Novatel.ROS_Time - cleanData.GPS_Novatel.ROS_Time(1);
    
    error_in_duration = durationsGPS - durationsROS;
    figure(303); plot(error_in_duration);
    % And the plot shows that they do not agree very well. However, both have
    % very accurate decimation:
    mean(diff(cleanData.GPS_Novatel.GPS_Time)) % Result is 0.050000107343119
    mean(diff(cleanData.GPS_Novatel.ROS_Time)) % Result is 0.050000566338432
    % Further, the decimation on the GPS never falls below 0.04 seconds, so
    % data does not seem to be lost:
    figure; plot(diff(durationsGPS))
    
    % To check for data loss,we can round the data to nearest decimation
    target_decimation = 0.01*round(100*mean(diff(cleanData.GPS_Novatel.GPS_Time)));
    % Now we check to see where the data would end otherwise. We use rounding
    % to see where data would land.
    endDurationGPS = target_decimation*round((1/target_decimation)*durationsGPS(end));
    % Do we get the length of data we expect?
    length(durationsGPS(:,1))
    endDurationGPS/target_decimation + 1 %#ok<NOPRT>
    % We have to add +1 since counting starts at zero.
    % Results here show that they agree.
    
    % Step 2: Force the GPS data to lock into the correct decimation
    durationsGPS_locked = (0:target_decimation:endDurationGPS);
    
    % Step 3: Find the offsets
    offsets(1) = cleanData.GPS_Novatel.GPS_Time(1)...
        - cleanData.GPS_Novatel.ROS_Time(1);
    offsets(2) = cleanData.IMU_Novatel.GPS_Time(1)...
        - cleanData.IMU_Novatel.ROS_Time(1);
    offset_to_add_to_ROS_to_get_GPS_time = mean(offsets);
    
    % Step 4: Decide a start time for all data (should this be a vote?)
    startTimeVotesGPS(1) = cleanData.GPS_Novatel.GPS_Time(1);
    startTimeVotesGPS(2) = cleanData.IMU_Novatel.GPS_Time(1);
    startTimeGPS = target_decimation *round(mean(startTimeVotesGPS)*(1/target_decimation)) %#ok<NOPRT>
    startTimeROS = startTimeGPS - offset_to_add_to_ROS_to_get_GPS_time %#ok<NOPRT>
    
    % Check result
    cleanData.IMU_Novatel.ROS_Time(1)
    cleanData.GPS_Novatel.ROS_Time(1)
    
end

%% Create vectors of trusted GPS data sources
% Step 1: Vote on a start and end time - NOTE: need to automated this
% instead of manually coding it!

% Find the available decimations
decimations(1) = cleanData.GPS_Novatel.GPS_Time_deltaT_target;
decimations(2) = cleanData.IMU_Novatel.GPS_Time_deltaT_target;
decimations(3) = cleanData.GPS_Hemisphere.GPS_Time_deltaT_target;
max_decimation = max(decimations);

% Find the start time in GPS time
startTimeVotesGPS(1) = fcn_findNearestDecimatedTime(cleanData.GPS_Novatel.GPS_Time(1,1),max_decimation);
startTimeVotesGPS(2) = fcn_findNearestDecimatedTime(cleanData.IMU_Novatel.GPS_Time(1,1),max_decimation);
startTimeVotesGPS(3) = fcn_findNearestDecimatedTime(cleanData.GPS_Hemisphere.GPS_Time(1,1),max_decimation);
cleanAndTimeAlignedData.Clocks.startTimeGPS = max(startTimeVotesGPS);  % Align to the GPS referenced sensor that turned on last

% Find the end time in GPS time
endTimeVotesGPS(1) = fcn_findNearestDecimatedTime(cleanData.GPS_Novatel.GPS_Time(end,1),max_decimation);
endTimeVotesGPS(2) = fcn_findNearestDecimatedTime(cleanData.IMU_Novatel.GPS_Time(end,1),max_decimation);
endTimeVotesGPS(3) = fcn_findNearestDecimatedTime(cleanData.GPS_Hemisphere.GPS_Time(end,1),max_decimation);
cleanAndTimeAlignedData.Clocks.endTimeGPS = min(endTimeVotesGPS);  % Align to the GPS referenced sensor that turned off first

% Find the equivalent points in ROS time
startTimeVotesROS(1) = fcn_findNearestDecimatedTime(cleanData.GPS_Novatel.ROS_Time(1,1),max_decimation);
startTimeVotesROS(2) = fcn_findNearestDecimatedTime(cleanData.IMU_Novatel.ROS_Time(1,1),max_decimation);
startTimeVotesROS(3) = fcn_findNearestDecimatedTime(cleanData.GPS_Hemisphere.ROS_Time(1,1),max_decimation);
cleanAndTimeAlignedData.Clocks.startTimeROS = max(startTimeVotesROS);  % Align to the GPS referenced sensor that turned on last

endTimeVotesROS(1) = fcn_findNearestDecimatedTime(cleanData.GPS_Novatel.ROS_Time(end,1),max_decimation);
endTimeVotesROS(2) = fcn_findNearestDecimatedTime(cleanData.IMU_Novatel.ROS_Time(end,1),max_decimation);
endTimeVotesROS(3) = fcn_findNearestDecimatedTime(cleanData.GPS_Hemisphere.ROS_Time(end,1),max_decimation);
cleanAndTimeAlignedData.Clocks.endTimeROS = min(endTimeVotesROS);  % Align to the GPS referenced sensor that turned off first

% Check that they are consistent with each other - and fix them if they
% are not
ROS_Time_target = (cleanAndTimeAlignedData.Clocks.startTimeROS:max_decimation:cleanAndTimeAlignedData.Clocks.endTimeROS)';
GPS_Time_target = (cleanAndTimeAlignedData.Clocks.startTimeGPS:max_decimation:cleanAndTimeAlignedData.Clocks.endTimeGPS)';

target_ROS_Npoints = round((cleanAndTimeAlignedData.Clocks.endTimeROS - cleanAndTimeAlignedData.Clocks.startTimeROS)/max_decimation)+1;
target_GPS_Npoints = round((cleanAndTimeAlignedData.Clocks.endTimeGPS - cleanAndTimeAlignedData.Clocks.startTimeGPS)/max_decimation)+1;

fprintf(1,'\nSTEP 1: we calculate the expected ROS and GPS start and end times, checking consistency:\n');
fprintf(1,'We do this by aggregating all the data that has GPS time, and looking at where\n');
fprintf(1,'each ROS and GPS time vector, for each data set starts and stops.\n');
fprintf(1,' We then take the latest start point for each, and earliest end points.\n');
fprintf(1,'If they are consistent with each other, we should end up with same number of points.\n');
fprintf(1,'\tNumber of expected points for ROS time zone via rounding: %g\n',target_ROS_Npoints);
fprintf(1,'\tNumber of expected points for GPS time zone via rounding: %g\n',target_GPS_Npoints);
fprintf(1,'\tNumber of expected points for ROS time zone via vector length: %g\n',length(ROS_Time_target));
fprintf(1,'\tNumber of expected points for GPS time zone via vector length: %g\n',length(GPS_Time_target));

% Go through each of the four cases

if (target_ROS_Npoints == target_GPS_Npoints)
    % Do nothing - this is good!
elseif (target_ROS_Npoints == (target_GPS_Npoints+2))
    % Case 1: there are 2 more ROS points than GPS points - nudge ROS points
    % inward on both ends
    cleanAndTimeAlignedData.Clocks.startTimeROS = cleanAndTimeAlignedData.Clocks.startTimeROS + max_decimation;
    cleanAndTimeAlignedData.Clocks.endTimeROS   = cleanAndTimeAlignedData.Clocks.endTimeROS - max_decimation;
elseif (target_ROS_Npoints == (target_GPS_Npoints+1))
    % Case 2: there are 1 more ROS points than GPS points - nudge ROS points
    % inward on just end
    cleanAndTimeAlignedData.Clocks.endTimeROS   = cleanAndTimeAlignedData.Clocks.endTimeROS - max_decimation;
elseif (target_ROS_Npoints == (target_GPS_Npoints-1))
    % Case 3: there are 1 more GPS points than ROS points - nudge GPS points
    % inward on just end
    cleanAndTimeAlignedData.Clocks.endTimeGPS   = cleanAndTimeAlignedData.Clocks.endTimeGPS - max_decimation;
elseif (target_ROS_Npoints == (target_GPS_Npoints-2))
    % Case 4: there are 2 more GPS points than ROS points - nudge GPS points
    % inward on both ends
    cleanAndTimeAlignedData.Clocks.startTimeGPS   = cleanAndTimeAlignedData.Clocks.startTimeGPS - max_decimation;
    cleanAndTimeAlignedData.Clocks.endTimeGPS   = cleanAndTimeAlignedData.Clocks.endTimeGPS - max_decimation;
else
    error('Time vectors do not align sufficiently between ROS and GPS time. Unable to continue');
end

% Check result again here - should be fixed now
ROS_Time_target = (cleanAndTimeAlignedData.Clocks.startTimeROS:max_decimation:cleanAndTimeAlignedData.Clocks.endTimeROS)';
GPS_Time_target = (cleanAndTimeAlignedData.Clocks.startTimeGPS:max_decimation:cleanAndTimeAlignedData.Clocks.endTimeGPS)';

target_ROS_Npoints = round((cleanAndTimeAlignedData.Clocks.endTimeROS - cleanAndTimeAlignedData.Clocks.startTimeROS)/max_decimation)+1;
target_GPS_Npoints = round((cleanAndTimeAlignedData.Clocks.endTimeGPS - cleanAndTimeAlignedData.Clocks.startTimeGPS)/max_decimation)+1;

fprintf(1,'\nIf the data is off by just two points, it is likely because the round-off. We then fix this:\n');
fprintf(1,'\tNumber of expected points for ROS time zone via rounding: %g\n',target_ROS_Npoints);
fprintf(1,'\tNumber of expected points for GPS time zone via rounding: %g\n',target_GPS_Npoints);
fprintf(1,'\tNumber of expected points for ROS time zone via vector length: %g\n',length(ROS_Time_target));
fprintf(1,'\tNumber of expected points for GPS time zone via vector length: %g\n',length(GPS_Time_target));


%% Step 2: 
% create decimations and check if vectors are missing data based on
% these start and end points for each data set. Basically, if the ROS time
% vector is not as long as the reference GPS time vector, then data must be
% missing somehwere.

fprintf(1,'\nSTEP 2: We now calculate the ideal time vectors in ROS and GPS time, for each\n');
fprintf(1,'\ncenti-second sampling interval we might see in the data (deltaT = 0.01 to 0.05).\n');
% First, calculate the reference target time vectors in GPS and ROS times
for i_deltaT = 1:5  % only go to 1 ms to 5 ms in current sampling methods
    deltaT = i_deltaT*0.01; % Convert to centi-seconds
    % Based on the delta_t for this data set, calculate what the ROS and
    % GPS time should be, as vectors.
    ROS_Time_target = (cleanAndTimeAlignedData.Clocks.startTimeROS:deltaT:cleanAndTimeAlignedData.Clocks.endTimeROS)';
    GPS_Time_target = (cleanAndTimeAlignedData.Clocks.startTimeGPS:deltaT:cleanAndTimeAlignedData.Clocks.endTimeGPS)';
    cleanAndTimeAlignedData.Clocks.targetTimeVector_ROS{i_deltaT} = ROS_Time_target;
    cleanAndTimeAlignedData.Clocks.targetTimeVector_GPS{i_deltaT} = GPS_Time_target;
    
    fprintf(1,'\nCaclulation results for time vectors assuming sampling time of %g centi-seconds: \n',i_deltaT);
    fprintf(1,'\tNumber of expected points for ROS time: %g\n',length(ROS_Time_target));
    fprintf(1,'\tNumber of expected points for GPS time: %g\n',length(GPS_Time_target));
    
    
end

% Now check how many fields do not have right number of data points
names = fieldnames(cleanData);
for i_data = 1:length(names)
    % Grab the data
    data_name = names{i_data};
    d = eval(cat(2,'cleanData.',data_name));
    
    % Calculate how many points should be there
    if isfield(d,'GPS_Time_deltaT_target')
        deltaT = d.GPS_Time_deltaT_target;
    elseif isfield(d,'ROS_Time_deltaT_target')
        deltaT = d.ROS_Time_deltaT_target;
    else
        error('Time target cannot be determined for: %s ... exiting.',data_name);
    end
    
    % Find the length of the reference vectors
    centi_seconds_samplingRate = round(deltaT*100);
    target_ROS_time = cleanAndTimeAlignedData.Clocks.targetTimeVector_ROS{centi_seconds_samplingRate};
    target_ROS_Npoints = length(target_ROS_time);
    
    % Only count indices within ROS range of times
    indGood = find((d.ROS_Time>=target_ROS_time(1,1))&(d.ROS_Time<=target_ROS_time(end,1)));
    
    % Calculate how many points are there, actually
    fprintf(1,'\nFor data set: %s\n',data_name);
    fprintf(1,'\tNumber of expected points ROS time points: %g\n',target_ROS_Npoints);
    fprintf(1,'\tNumber of actual points: %d\n',length(d.ROS_Time(indGood,1))); %#ok<FNDSB>
    
end

%% Step 3: Trim data to above limits
fprintf(1,'\nSTEP 3: For each data set in cleanData structure, we find if there is missing \n');
fprintf(1,'data relative to the reference GPS or ROS time vectors, given the sampling rate for each \n');
fprintf(1,'data. \n');
names = fieldnames(cleanData);
for i_data = 1:length(names)
    % Grab the data
    data_name = names{i_data};
    d = cleanData.(data_name);
    
    % Update the user
    if flag_do_debug
        fprintf(1,'\n Sensor %d of %d: ',i_data,length(names));
        fprintf(1,'Re-aligning time for sensor: %s\n',data_name);
    end

    
    % Calculate the deltaT we want, by preferring GPS time, but defaulting
    % to ROS time if GPS_time is not available. Use these to calculate the
    % deltaT that we'll use in the steps that follow.
    if isfield(d,'GPS_Time_deltaT_target')
        deltaT = d.GPS_Time_deltaT_target;
    elseif isfield(d,'ROS_Time_deltaT_target')
        deltaT = d.ROS_Time_deltaT_target;
    else
        error('Time target cannot be determined for: %s ... exiting.',data_name);
    end


    % Find the length of the reference vectors
    centi_seconds_samplingRate = round(deltaT*100);
    target_GPS_time = cleanAndTimeAlignedData.Clocks.targetTimeVector_GPS{centi_seconds_samplingRate};
    target_GPS_Npoints = length(target_GPS_time);
    target_ROS_time = cleanAndTimeAlignedData.Clocks.targetTimeVector_ROS{centi_seconds_samplingRate};
    target_ROS_Npoints = length(target_ROS_time);
    
    % Define which vector we are going to use as a reference: GPS time or
    % ROS time. We prefer GPS time, so default to this, but this is not
    % available on all data sources so we need to check.
    flag_using_GPS_time = 0;
    if isfield(d,'GPS_Time') && ~any(isnan(d.GPS_Time))
        % There is no NaN in the GPS data - use this
        time_to_search = fcn_findNearestDecimatedTime(d.GPS_Time,deltaT);
        %time_to_search = d.GPS_Time;
        time_target = target_GPS_time;
        flag_using_GPS_time = 1;
        
    else
        % Use the ROS time data
        time_to_search = fcn_findNearestDecimatedTime(d.ROS_Time,deltaT);
        %time_to_search = d.ROS_Time;
        time_target = target_ROS_time;
    end
    
    % In the search that follows, MATLAB rounds more easily to integers. So
    % we convert the target time (i.e. the time sequence we want) and our
    % search time (i.e. the time sequence we have) to integer values. The
    % results later will simply be indices that match one to the other, so
    % it will not matter that we are using integer times.
    rounded_time_to_target = round(time_target/deltaT);
    rounded_time_to_search = round(time_to_search/deltaT);
    
    %     % For debugging:
    %     rounded_time_to_target = rounded_time_to_target(141570:141580,1);
    %     index_start = find(rounded_time_to_search>rounded_time_to_target(1),1);
    %     index_end = find(rounded_time_to_search>rounded_time_to_target(end,1));
    %     rounded_time_to_search = rounded_time_to_search(index_start:index_end,1);
    %
    
    
    % Report results thus far
    fprintf(1,'\nFor data set: %s\n',data_name);
    fprintf(1,'\tUsing GPS time?: %g\n',flag_using_GPS_time);
    fprintf(1,'\tTargeting deltaT of: %g\n',deltaT);
    fprintf(1,'\tNumber of expected points: %g\n',length(rounded_time_to_target));
    
    %% Target times are determined per above, now shift data
    % METHOD 1: (very fast) use rounding and using array
    % intersection,
    [~,good_target_time_indices_intersect,good_search_time_indices_intersect] = intersect(...
        rounded_time_to_target,...
        rounded_time_to_search);
    zeroed_search = rounded_time_to_search;
    zeroed_search(good_search_time_indices_intersect) = 0;
    bad_search_time_indices_intersect = find(zeroed_search>0);
    
    zeroed_target = rounded_time_to_target;
    zeroed_target(good_target_time_indices_intersect) = 0;
    bad_target_time_indices_intersect = find(zeroed_target>0);
    
    if 1==0  % THIS IS VERY SLOW - DO NOT UNCOMMENT UNLESS VERY PATIENT
        % METHOD 2: (VERY slow)
        % Find the indices of the time_to_search that fits within each
        % sampling interval given by the target_time
        indices_good_time = NaN*rounded_time_to_target; % Initialize indices to NaN
        
        start_of_search = 1;
        for i_time = 1:length(rounded_time_to_target)
            t_start = rounded_time_to_target(i_time);
            t_end = t_start+deltaT;
            t_search = rounded_time_to_search;
            
            % ind_there = find((rounded_time_to_search>=t_start)&(rounded_time_to_search<t_end),1);
            % ind_there = find((t_search>=t_start)&(t_search<t_end),1);
            ind_there = find(t_search>=t_start,1,'first');
            if ~isempty(ind_there)
                if t_search(ind_there)<t_end
                    indices_good_time(i_time) = ind_there;
                    start_of_search = ind_there;
                else
                    ind_there = {};
                end
            end
            if isempty(ind_there)
                if ~ismember(i_time,bad_target_time_indices_intersect)
                    % The only way to enter here is if there's diagreement
                    % between the intersection method and this brute force
                    % method. This allows us to compare Method 1 and 2
                    
                    fprintf(1,'Unable to find samples at i_time: %d\n',i_time);
                    fprintf(1,'seacrhing from: %f\n',t_start);
                    fprintf(1,'seacrhing to: %f\n',t_end);
                    fprintf(1,'first one before is: %f\n,',...
                        rounded_time_to_search(find(rounded_time_to_search<t_end,1,'last')));
                    fprintf(1,'first one after is: %f\n,',...
                        rounded_time_to_search(find(rounded_time_to_search>t_start,1,'first')));
                    fprintf(1,'rounded: first one before is: %f\n,',...
                        rounded_time_to_search(find(rounded_rounded_time_to_search<(t_end/deltaT),1,'last')));
                    fprintf(1,'rounded: first one after is: %f\n,',...
                        rounded_time_to_search(find(rounded_time_to_search>(t_start/deltaT),1,'first')));
                    
                    fprintf(1,'\n');
                end
            end
        end
        
        % Tag the good and bad times
        good_search_time_indices_brute = indices_good_time(indices_good_time>0);
        zeroed_search = rounded_time_to_search;
        zeroed_search(good_search_time_indices_brute) = 0;
        bad_search_time_indices_brute = find(zeroed_search>0);
        
        
        % Debugging display
        fprintf(1,'\tNumber of actual points: %d\n',length(good_search_time_indices_brute));
        if isequal(good_search_time_indices_intersect,good_search_time_indices_brute)
            fprintf(1,'\tThe good indices found by brute force matched that of rounding.\n');
        else
            fprintf(1,'\tThe good indices found by brute force does NOT match that of rounding.\n');
        end
        if isequal(bad_search_time_indices_intersect,bad_search_time_indices_brute)
            fprintf(1,'\tThe bad indices found by brute force matched that of rounding.\n');
        else
            fprintf(1,'\tThe bad indices found by brute force does NOT match that of rounding.\n');
        end
    end
    % The following is for debugging, to literally see the time jumps
    if 1==0
        differences = diff([0; bad_search_time_indices_intersect]) + diff([bad_search_time_indices_intersect; 0]);
        location_of_jumps = find(abs(differences)>2);
        for i_location = 1:length(location_of_jumps)
            i_bad = location_of_jumps(i_location);
            current_bad_index = bad_search_time_indices_intersect(i_bad);
            before = current_bad_index - 3;
            after = current_bad_index + 3;
            
            % Limit before and after to allowable ranges
            before = max(before,1);
            after = min(after,length(rounded_time_to_search));
            
            % Data jump definition
            time_jump = rounded_time_to_search(before:after);
            
            % Show results
            fprintf(1,'\nFor time jump detected at position: %g\n',current_bad_index);
            fprintf(1,'\t\tTimes at time jump: %f\n',time_jump);
            fprintf(1,'\n');
        end
    end
    
    
    % Step 4: now that the indices are defined...
    % Fill the clean and time aligned data structure
    
    % Check to see if the GPS time is missing, but ROS time is there. If
    % so, need to create GPS time.
    if ~isfield(d,'GPS_Time') && isfield(d,'ROS_Time')
        cleanAndTimeAlignedData.(data_name).GPS_Time = target_GPS_time;
        
        cleanAndTimeAlignedData.(data_name).GPS_Time_Sigma =...
            fcn_calcSigmaNoOutliers(d.ROS_Time(good_search_time_indices_intersect));

        if flag_do_debug
            fprintf(1,'\tWARNING: Only ROS time was detected, so GPS time set to locked-in time.\n');
        end
    end
    
    % Loop through remaining subfields, fixing each if necessary
    subfieldNames = fieldnames(d); % Grab all the subfields
    for i_subField = 1:length(subfieldNames)
               
        % Grab the name of the ith subfield
        subFieldName = subfieldNames{i_subField};
        
        if flag_do_debug
            fprintf(1,'\tProcessing subfield: %s ',subFieldName);
        end
        
        if strcmp(subFieldName,'GPS_Time')
            cleanAndTimeAlignedData.(data_name).GPS_Time = target_GPS_time;
            if flag_do_debug
                fprintf(1,'<-- Replaced with locked-in time\n');
            end
        elseif strcmp(subFieldName,'GPS_Time_Sigma')
            GPS_Time_Sigma = fcn_calcSigmaNoOutliers(d.GPS_Time(good_search_time_indices_intersect));
            cleanAndTimeAlignedData.(data_name).GPS_Time_Sigma = GPS_Time_Sigma;
            if flag_do_debug
                fprintf(1,'<-- Replaced with locked-in time\n');
            end
        elseif strcmp(subFieldName,'ROS_Time')
            cleanAndTimeAlignedData.(data_name).ROS_Time = target_ROS_time;
            if flag_do_debug
                fprintf(1,'<-- Replaced with locked-in time\n');
            end
        elseif strcmp(subFieldName,'ROS_Time_Sigma')
            ROS_Time_Sigma = fcn_calcSigmaNoOutliers(d.ROS_Time(good_search_time_indices_intersect));
            cleanAndTimeAlignedData.(data_name).ROS_Time_Sigma = ROS_Time_Sigma;
            if flag_do_debug
                fprintf(1,'<-- Replaced with locked-in time\n');
            end
        else
            
            % Fill in the data from this subfield
            dataInSubfield = d.(subFieldName);            
            
            % If the data is just a scalar, do nothing
            if length(dataInSubfield(:,1))==1
                dataInGPSTime = dataInSubfield;  %#ok<*NASGU>
            else
                % This is a vector of data - need to resize
                timeInUnevenGPSTime = target_GPS_time(good_target_time_indices_intersect);
                dataInUnevenGPSTime = dataInSubfield(good_search_time_indices_intersect);
                
                % Resample the data - put zero in if extrapolating
                if any(strcmp(subFieldName,fields_to_interpolate_to_nearest_neighbor))
                    % This is data that is not a floating type - just use
                    % nearest neighbor. Need to typecast to allow interp to
                    % work. Note: typecasting only works if there's no NaN
                    if any(isnan(dataInUnevenGPSTime))
                        dataInUnevenGPSTime = nan*timeInUnevenGPSTime;
                    else
                        dataInUnevenGPSTime = 1.0*dataInUnevenGPSTime;
                        dataInUnevenGPSTime = cast(dataInUnevenGPSTime,'single');
                    end
                    try
                        dataInGPSTime = interp1(...
                            timeInUnevenGPSTime,...
                            dataInUnevenGPSTime,target_GPS_time,'nearest');
                    
                    catch
                        error('debug here');
                    end

                        
                else  % This is normal data
                    % Use linear interpolation instead
                    
                    if ~all(isnan(dataInUnevenGPSTime))
                        if ~isa(dataInUnevenGPSTime,'float')
                            if isa(dataInUnevenGPSTime,'integer')
                                dataInUnevenGPSTime = typecast(dataInUnevenGPSTime,'double');
                            else
                                error('Unknown data type encountered');
                            end
                        end
                        
                        % Do the interpolation
                        try
                            dataInGPSTime = interp1(...
                                timeInUnevenGPSTime,...
                                dataInUnevenGPSTime,target_GPS_time,'linear','extrap');
                        catch
                            disp('Debug here');
                            pause;
                        end
                    else
                        error('Data is nothing but NaN - unable to interpolate.');
                    end
                end

                
                % Push the data into new data structure - 1 means that the
                % original data was good, 0 means that it was interpolated
                GPS_Time_goodIndices = 0*target_GPS_time;
                GPS_Time_goodIndices(good_target_time_indices_intersect) = 1;
                cleanAndTimeAlignedData.(data_name).GPS_Time_goodIndices = GPS_Time_goodIndices;
                
            end % Ends if check to see if data is a scalar
            
            % Push the subfield data into new data structure
            cleanAndTimeAlignedData.(data_name).(subFieldName) = dataInGPSTime;            
            if flag_do_debug
                fprintf(1,'<-- Replaced with interpolated values\n');
            end
        end
    end
end

%  %disp('Check data here');
%
% %
% % % Step 4: Find the offsets from ROS to GPS
% % % Find offsets using only data that is trusted (need to automate this in
% % % the future)
% % cleanAndTimeAlignedData.offset_to_add_to_ROS_to_get_GPS_time = cleanAndTimeAlignedData.Clocks.startTimeGPS - cleanAndTimeAlignedData.Clocks.startTimeROS;
% % offsets(:,1) = cleanData.GPS_Novatel.GPS_Time - cleanData.GPS_Novatel.ROS_Time;
% %
% % % Low-pass filter the data
% % deltaT = cleanData.GPS_Novatel.GPS_Time_deltaT;  % Update this if above is automated
% % f_sample = 1/deltaT;
% % f_Nyquist = f_sample/2;
% % f_cutoff = 0.5;
% % [b,a] = butter(2,f_cutoff/f_Nyquist);
% % ROS_to_GPS_timeOffset_filtered = filtfilt(b,a,offsets);
% % ROS_to_GPS_timeOffset_filtered_movmeaned = movmean(ROS_to_GPS_timeOffset_filtered,251);
% %
% % % % For debugging: the following shows if the filter worked
% % if 1==0
% %     figure(847464);
% %     plot(offsets - offsets(1),'b');
% %     xlabel('Index'); ylabel('Offset from Ros to GPS (t_ROS - t_GPS)');
% %     hold on;
% %     plot(ROS_to_GPS_timeOffset_filtered - offsets(1),'g');
% %     plot(ROS_to_GPS_timeOffset_filtered_movmeaned - offsets(1),'r');
% %     legend('Raw data','Filtered','Filtered and central averaged');
% % end
% %
% %
% %
% % Step 5: move ROS times to GPS time locked-in decimation. Before doing
% % this, ensure that this conversion aligns correctly by doing a
% % pseudo-conversion of ROS time first, making sure it's still ordered.
%
% deltaT = cleanData.GPS_Novatel.GPS_Time_deltaT;  % Update this if above is automated, but this is most trusted value
%
% if 1==0
%     % Test the result on dummy data (known)
%     target_GPS_time_reference = (cleanAndTimeAlignedData.Clocks.startTimeGPS:deltaT:(cleanAndTimeAlignedData.Clocks.endTimeGPS+deltaT))';
%     Npoints = length(target_GPS_time_reference(:,1));
%     %target_ROS_time_reference = target_GPS_time_reference - ROS_to_GPS_timeOffset_filtered_movmeaned(1:Npoints,1);
%     target_ROS_time_reference = target_GPS_time_reference - ROS_to_GPS_timeOffset_filtered(1:Npoints,1);
%     X  = cleanData.GPS_Novatel.ROS_Time;
%     V  = cleanData.GPS_Novatel.GPS_Time;
%     Xq = target_ROS_time_reference;
%     Vq = interp1(X,V,Xq);
%     predicted_GPS_time = fcn_findNearestDecimatedTime(Vq,deltaT);
%
%     error_in_time = predicted_GPS_time - cleanData.GPS_Novatel.GPS_Time;
%     figure(3838); clf; plot(error_in_time);
%
%     % If the predicted GPS data at this point is correctly sorted, then it
%     % means that our conversion from ROS time to GPS time produced data that is
%     % uniquely increasing, within the decimation accuracy given by deltaT.
%     % Thus, we should use the IDEAL GPS time, rather than the predicted GPS
%     % time
%     if(issorted(predicted_GPS_time))
%         disp('Stopped here');
%     end
% end


end

function goal_decimated = fcn_findNearestDecimatedTime(goal_time,decimation)
tolerance = 0.000001; % If within one microsecond, then it is OK.
goal_decimated = decimation *floor((goal_time+tolerance)/decimation);
end

function real_sigma = fcn_calcSigmaNoOutliers(data)
differences = diff(data);
deviations = differences - mean(differences);
outlier_sigma = std(deviations);
% Reject outliers
deviations_with_no_outliers = deviations(abs(deviations)<(3*outlier_sigma));
real_sigma = std(deviations_with_no_outliers);
end
