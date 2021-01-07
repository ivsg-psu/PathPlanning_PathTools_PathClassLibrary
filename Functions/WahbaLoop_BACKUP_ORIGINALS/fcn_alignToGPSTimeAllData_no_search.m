function cleanAndTimeAlignedData = fcn_alignToGPSTimeAllData(cleanData)

% This function is used to align in time all the data such that their ROS
% and GPS time stamps are uniform in time, and they all have the same start and
% end times by being trimmed to same points in time.
%
% Author: Sean Brennan
% Date: 2019_10_06
% modify Date: 2019_10_06
%
% Updates:
%  2019_10_06 - first writing by Dr. B
% 
% To do items:
% (as of 2019_10_10 - Step 3 can be greatly simplified using an intersect
% command for all the cases.

flag_showMotivation = 0;

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

%% Create a trusted GPS data sources
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
cleanAndTimeAlignedData.startTimeGPS = max(startTimeVotesGPS);  % Align to the GPS referenced sensor that turned on last

% Find the end time in GPS time
endTimeVotesGPS(1) = fcn_findNearestDecimatedTime(cleanData.GPS_Novatel.GPS_Time(end,1),max_decimation);
endTimeVotesGPS(2) = fcn_findNearestDecimatedTime(cleanData.IMU_Novatel.GPS_Time(end,1),max_decimation);
endTimeVotesGPS(3) = fcn_findNearestDecimatedTime(cleanData.GPS_Hemisphere.GPS_Time(end,1),max_decimation);
cleanAndTimeAlignedData.endTimeGPS = min(endTimeVotesGPS);  % Align to the GPS referenced sensor that turned off first

% Find the equivalent points in ROS time
startTimeVotesROS(1) = fcn_findNearestDecimatedTime(cleanData.GPS_Novatel.ROS_Time(1,1),max_decimation);
startTimeVotesROS(2) = fcn_findNearestDecimatedTime(cleanData.IMU_Novatel.ROS_Time(1,1),max_decimation);
startTimeVotesROS(3) = fcn_findNearestDecimatedTime(cleanData.GPS_Hemisphere.ROS_Time(1,1),max_decimation);
cleanAndTimeAlignedData.startTimeROS = max(startTimeVotesROS);  % Align to the GPS referenced sensor that turned on last

endTimeVotesROS(1) = fcn_findNearestDecimatedTime(cleanData.GPS_Novatel.ROS_Time(end,1),max_decimation);
endTimeVotesROS(2) = fcn_findNearestDecimatedTime(cleanData.IMU_Novatel.ROS_Time(end,1),max_decimation);
endTimeVotesROS(3) = fcn_findNearestDecimatedTime(cleanData.GPS_Hemisphere.ROS_Time(end,1),max_decimation);
cleanAndTimeAlignedData.endTimeROS = min(endTimeVotesROS);  % Align to the GPS referenced sensor that turned off first

% Check that they are consistent with each other - and fix them if they
% are not
ROS_Time_target = (cleanAndTimeAlignedData.startTimeROS:max_decimation:cleanAndTimeAlignedData.endTimeROS)';
GPS_Time_target = (cleanAndTimeAlignedData.startTimeGPS:max_decimation:cleanAndTimeAlignedData.endTimeGPS)';

target_ROS_Npoints = round((cleanAndTimeAlignedData.endTimeROS - cleanAndTimeAlignedData.startTimeROS)/max_decimation)+1;
target_GPS_Npoints = round((cleanAndTimeAlignedData.endTimeGPS - cleanAndTimeAlignedData.startTimeGPS)/max_decimation)+1;

fprintf(1,'\nFinished calculating ROS and GPS start and end times, checking consistency:\n');
fprintf(1,'\tNumber of expected points for ROS time zone via rounding: %g\n',target_ROS_Npoints);
fprintf(1,'\tNumber of expected points for GPS time zone via rounding: %g\n',target_GPS_Npoints);
fprintf(1,'\tNumber of expected points for ROS time zone via vector length: %g\n',length(ROS_Time_target));
fprintf(1,'\tNumber of expected points for GPS time zone via vector length: %g\n',length(GPS_Time_target));

% Go through each of the cases

if (target_ROS_Npoints == target_GPS_Npoints)
    % Do nothing - this is good!
elseif (target_ROS_Npoints == (target_GPS_Npoints+2))
    % Case 1: there are 2 more ROS points than GPS points - nudge ROS points
    % inward on both ends
    cleanAndTimeAlignedData.startTimeROS = cleanAndTimeAlignedData.startTimeROS + max_decimation;
    cleanAndTimeAlignedData.endTimeROS   = cleanAndTimeAlignedData.endTimeROS - max_decimation;
elseif (target_ROS_Npoints == (target_GPS_Npoints+1))
    % Case 2: there are 1 more ROS points than GPS points - nudge ROS points
    % inward on just end
    cleanAndTimeAlignedData.endTimeROS   = cleanAndTimeAlignedData.endTimeROS - max_decimation;
elseif (target_ROS_Npoints == (target_GPS_Npoints-1))
    % Case 3: there are 1 more GPS points than ROS points - nudge GPS points
    % inward on just end
    cleanAndTimeAlignedData.endTimeGPS   = cleanAndTimeAlignedData.endTimeGPS - max_decimation;
elseif (target_ROS_Npoints == (target_GPS_Npoints-2))
    % Case 4: there are 2 more GPS points than ROS points - nudge GPS points
    % inward on both ends
    cleanAndTimeAlignedData.startTimeGPS   = cleanAndTimeAlignedData.startTimeGPS - max_decimation;
    cleanAndTimeAlignedData.endTimeGPS   = cleanAndTimeAlignedData.endTimeGPS - max_decimation;
else
    error('Time vectors do not align sufficiently between ROS and GPS time. Unable to continue');
end

% Check result again here - should be fixed now
ROS_Time_target = (cleanAndTimeAlignedData.startTimeROS:max_decimation:cleanAndTimeAlignedData.endTimeROS)';
GPS_Time_target = (cleanAndTimeAlignedData.startTimeGPS:max_decimation:cleanAndTimeAlignedData.endTimeGPS)';

target_ROS_Npoints = round((cleanAndTimeAlignedData.endTimeROS - cleanAndTimeAlignedData.startTimeROS)/max_decimation)+1;
target_GPS_Npoints = round((cleanAndTimeAlignedData.endTimeGPS - cleanAndTimeAlignedData.startTimeGPS)/max_decimation)+1;

fprintf(1,'\nFinished corrections for consistency of ROS and GPS start and end times, checking consistency again:\n');
fprintf(1,'\tNumber of expected points for ROS time zone via rounding: %g\n',target_ROS_Npoints);
fprintf(1,'\tNumber of expected points for GPS time zone via rounding: %g\n',target_GPS_Npoints);
fprintf(1,'\tNumber of expected points for ROS time zone via vector length: %g\n',length(ROS_Time_target));
fprintf(1,'\tNumber of expected points for GPS time zone via vector length: %g\n',length(GPS_Time_target));




% Step 2: create decimations and check if vectors are missing data based on
% these start and end points for each data set. Basically, if the ROS time
% vector is not as long as the reference GPS time vector, then data must be
% missing somehwere.

names = fieldnames(cleanData);
for i_data = 1:length(names)
    % Grab the data
    data_name = names{i_data};
    d = eval(cat(2,'cleanData.',data_name));
    
    % Calculate how many points should be there
    deltaT = d.GPS_Time_deltaT_target;
    if isnan(deltaT)  % If the GPS time is empty, it will give an NaN
        deltaT = d.ROS_Time_deltaT_target;
    end
    
    % Based on the delta_t for this data set, calculate what the ROS and
    % GPS time should be, as vectors.
    ROS_Time_target = (cleanAndTimeAlignedData.startTimeROS:deltaT:cleanAndTimeAlignedData.endTimeROS)';
    GPS_Time_target = (cleanAndTimeAlignedData.startTimeGPS:deltaT:cleanAndTimeAlignedData.endTimeGPS)';
    
    % Find the length of these vectors
    target_ROS_Npoints = length(ROS_Time_target);
    target_GPS_Npoints = length(GPS_Time_target);
    
    % Calculate how many points are there
    fprintf(1,'\nFor data set: %s\n',data_name);
    fprintf(1,'\tNumber of expected points GPS time points: %g\n',target_GPS_Npoints);
    fprintf(1,'\tNumber of expected points ROS time points: %g\n',target_ROS_Npoints);
    fprintf(1,'\tNumber of actual points: %d\n',length(d.ROS_Time(:,1)));
    
end

% Step 3: Trim data to above limits
names = fieldnames(cleanData);
for i_data = 1:length(names)
    % Grab the data
    data_name = names{i_data};
    d = eval(cat(2,'cleanData.',data_name));
    
    % Calculate how many points should be there
    deltaT = d.GPS_Time_deltaT_target;
    if isnan(deltaT)
        deltaT = d.ROS_Time_deltaT_target;
    end
    
    % Create vectors representing desired ROS and GPD data
    ROS_Time_target = (cleanAndTimeAlignedData.startTimeROS:deltaT:cleanAndTimeAlignedData.endTimeROS)';
    GPS_Time_target = (cleanAndTimeAlignedData.startTimeGPS:deltaT:cleanAndTimeAlignedData.endTimeGPS)';
    
    target_ROS_Npoints = length(ROS_Time_target);
    target_GPS_Npoints = length(GPS_Time_target);
    
    
    % Calculate how many points are actually there in the ROS time data,
    % using the start and end ROS time trim points
    [closest_ROS_time_start_index] = find(d.ROS_Time>=cleanAndTimeAlignedData.startTimeROS,1);
    [closest_ROS_time_end_index] = find(d.ROS_Time>=cleanAndTimeAlignedData.endTimeROS,1);
    
    if isempty(closest_ROS_time_start_index)
        closest_ROS_time_start_index = 1;
    end
    if isempty(closest_ROS_time_end_index)
        closest_ROS_time_end_index = length(d.ROS_Time(:,1));
    end
    
    actual_Npoints = closest_ROS_time_end_index - closest_ROS_time_start_index+1;
    
    if 1==1  % Print results?
        fprintf(1,'\nFor data set: %s\n',data_name);
        fprintf(1,'\tNumber of expected GPS points: %d\n',target_GPS_Npoints);
        fprintf(1,'\tNumber of expected ROS points: %d\n',target_ROS_Npoints);
        fprintf(1,'\tNumber of actual points via ROS time after trimming: %d\n',actual_Npoints);
    end
    
    flag_dataGood = 0;
    indices_to_use = (closest_ROS_time_start_index:closest_ROS_time_end_index)';
    if (target_ROS_Npoints == actual_Npoints)
        % Perfect match - data is good as is
        flag_dataGood = 1;
        
    elseif (target_ROS_Npoints == (actual_Npoints-2))
        % There are two points too many in the ROS time. Trim from both
        % ends.
        flag_dataGood = 1;
        indices_to_use = indices_to_use(2:end-1);
    elseif (target_ROS_Npoints == (actual_Npoints-1))
        % There are one points too many in the ROS time. Trim from both
        % ends. Need to take one away. Both errors below will be negative.
        % The one that is least negative should be removed.
        
        error_startTime = cleanAndTimeAlignedData.startTimeROS - d.ROS_Time(closest_ROS_time_start_index);
        error_endTime = cleanAndTimeAlignedData.endTimeROS - d.ROS_Time(closest_ROS_time_end_index);
        
        if error_startTime < error_endTime
            % Move the  end time forward
            indices_to_use = indices_to_use(1:end-1);
        else
            indices_to_use = indices_to_use(2:end);
        end
        flag_dataGood = 1;
        
    elseif (target_ROS_Npoints == (actual_Npoints+1))
        % There is one point too few in the ROS time vector. Need to add.
        % Find which end of the ROS time is furthest away from the goal,
        % and trim that end. Both of these errors will be negative, because
        % the find function takes the first actual ROS time point after the
        % desired mark.
        error_startTime = cleanAndTimeAlignedData.startTimeROS - d.ROS_Time(closest_ROS_time_start_index);
        error_endTime = cleanAndTimeAlignedData.endTimeROS - d.ROS_Time(closest_ROS_time_end_index);
        
        if error_startTime < error_endTime
            % Move the start time forward
            indices_to_use = [max(1,closest_ROS_time_start_index-1);
                indices_to_use];
        else
            indices_to_use = [ indices_to_use;
                min(length(d.ROS_Time),closest_ROS_time_end_index+1);];
        end
        
        flag_dataGood = 1;
        
    elseif (target_ROS_Npoints == (actual_Npoints+2))
        % There is two points too fewmin the ROS time vector. Add one from
        % both sides
        flag_dataGood = 1;
        indices_to_use = [...
            max(1,closest_ROS_time_start_index-1);
            indices_to_use;
            min(length(d.ROS_Time),closest_ROS_time_end_index+1);]; %#ok<*AGROW>
    else
        warning('severe mismatch detected in data');
        if any(isnan(d.GPS_Time))
            %             start_index = find(d.ROS_Time>cleanAndTimeAlignedData.startTimeROS,1);
            %             end_index = length(d.ROS_Time);
            
            offset = d.ROS_Time(1,1) - cleanAndTimeAlignedData.startTimeROS;
            offset_decimated = fcn_findNearestDecimatedTime(offset,d.ROS_Time_deltaT_target);
            endTime = d.ROS_Time_deltaT_target*(length(d.ROS_Time(:,1))-1);
            d.GPS_Time = cleanAndTimeAlignedData.startTimeGPS + ...
                offset_decimated + ...
                + (0:d.ROS_Time_deltaT_target:endTime)';
            d.GPS_Time_deltaT_target = d.ROS_Time_deltaT_target;
        end
        
        roundedDataTime = round(d.GPS_Time*(1/d.GPS_Time_deltaT_target));
        roundedTargetTime = round((1/d.GPS_Time_deltaT_target)*GPS_Time_target);
        
        % For debugging
        %[roundedTime(1:30) roundedTargetTime(1:30)]
        
        %         % Find the values of existing GPS time that are within the target
        %         % time
        %         index_start = find(roundedTime==cleanAndTimeAlignedData.startTimeGPS);
        %         index_end   = find(roundedTime==cleanAndTimeAlignedData.endTimeGPS);
        
        % Reconstruct a time vector by storing which indices we'll use
        [mergedTime,indices_GPS_target,indices_to_use] = ...
            intersect(roundedTargetTime,roundedDataTime);
        mergedTime = mergedTime*d.GPS_Time_deltaT_target;

        if 1==0
            length(mergedTime);
            length(indices_GPS_target);
        end
        
        
    end
    fprintf(1,'\tNumber of in indices after trimming: %d\n',length(indices_to_use));
    
    % Step 4: now that the indices are defined...
    % Fill the clean and time aligned data structure
    
    subfieldNames = fieldnames(d);
    for i_subField = 1:length(subfieldNames)
        commandString2 = cat(2,'d.',subfieldNames{i_subField});
        dataInSubfield = eval(commandString2);
        if length(dataInSubfield(:,1))==1
            dataInGPSTIme = dataInSubfield; 
        else
            % This is a vector of data
            dataInGPSTIme = dataInSubfield(indices_to_use);
            
            if 1==flag_dataGood
                indices_GPS_target = (1:length(GPS_Time_target))'; %#ok<*NASGU>
            else
                % If enter here, then data is not uniformly spaced and
                % needs to be interpolated.
                dataInGPSTIme = interp1(mergedTime, dataInGPSTIme,GPS_Time_target,'linear','extrap');

                
            end
            % Push the data into new data structure
            goodIndices = 0*GPS_Time_target;
            commandString5 = cat(2,...
                'cleanAndTimeAlignedData.',...
                data_name,...
                '.GPS_Time_goodIndices = indices_GPS_target;');
            eval(commandString5);
 
        end
        % Push the data into new data structure
        commandString3 = cat(2,...
            'cleanAndTimeAlignedData.',...
            data_name,...
            '.',subfieldNames{i_subField},...
            '= dataInGPSTIme;');
        eval(commandString3);
    end
    commandString4 = cat(2,...
        'cleanAndTimeAlignedData.',...
        data_name,...
        '.GPS_Time = GPS_Time_target;');
    eval(commandString4);
    
end

 %disp('Check data here');

%
% % Step 4: Find the offsets from ROS to GPS
% % Find offsets using only data that is trusted (need to automate this in
% % the future)
% cleanAndTimeAlignedData.offset_to_add_to_ROS_to_get_GPS_time = cleanAndTimeAlignedData.startTimeGPS - cleanAndTimeAlignedData.startTimeROS;
% offsets(:,1) = cleanData.GPS_Novatel.GPS_Time - cleanData.GPS_Novatel.ROS_Time;
%
% % Low-pass filter the data
% deltaT = cleanData.GPS_Novatel.GPS_Time_deltaT;  % Update this if above is automated
% f_sample = 1/deltaT;
% f_Nyquist = f_sample/2;
% f_cutoff = 0.5;
% [b,a] = butter(2,f_cutoff/f_Nyquist);
% ROS_to_GPS_timeOffset_filtered = filtfilt(b,a,offsets);
% ROS_to_GPS_timeOffset_filtered_movmeaned = movmean(ROS_to_GPS_timeOffset_filtered,251);
%
% % % For debugging: the following shows if the filter worked
% if 1==0
%     figure(847464);
%     plot(offsets - offsets(1),'b');
%     xlabel('Index'); ylabel('Offset from Ros to GPS (t_ROS - t_GPS)');
%     hold on;
%     plot(ROS_to_GPS_timeOffset_filtered - offsets(1),'g');
%     plot(ROS_to_GPS_timeOffset_filtered_movmeaned - offsets(1),'r');
%     legend('Raw data','Filtered','Filtered and central averaged');
% end
%
%
%
% Step 5: move ROS times to GPS time locked-in decimation. Before doing
% this, ensure that this conversion aligns correctly by doing a
% pseudo-conversion of ROS time first, making sure it's still ordered.

deltaT = cleanData.GPS_Novatel.GPS_Time_deltaT;  % Update this if above is automated, but this is most trusted value

if 1==0
    % Test the result on dummy data (known)
    target_GPS_time_reference = (cleanAndTimeAlignedData.startTimeGPS:deltaT:(cleanAndTimeAlignedData.endTimeGPS+deltaT))';
    Npoints = length(target_GPS_time_reference(:,1));
    %target_ROS_time_reference = target_GPS_time_reference - ROS_to_GPS_timeOffset_filtered_movmeaned(1:Npoints,1);
    target_ROS_time_reference = target_GPS_time_reference - ROS_to_GPS_timeOffset_filtered(1:Npoints,1);
    X  = cleanData.GPS_Novatel.ROS_Time;
    V  = cleanData.GPS_Novatel.GPS_Time;
    Xq = target_ROS_time_reference;
    Vq = interp1(X,V,Xq);
    predicted_GPS_time = fcn_findNearestDecimatedTime(Vq,deltaT);
    
    error_in_time = predicted_GPS_time - cleanData.GPS_Novatel.GPS_Time;
    figure(3838); clf; plot(error_in_time);
    
    % If the predicted GPS data at this point is correctly sorted, then it
    % means that our conversion from ROS time to GPS time produced data that is
    % uniquely increasing, within the decimation accuracy given by deltaT.
    % Thus, we should use the IDEAL GPS time, rather than the predicted GPS
    % time
    if(issorted(predicted_GPS_time))
        disp('Stopped here');
    end
end


% Step 6: Fix missing time steps?




end

function goal_decimated = fcn_findNearestDecimatedTime(goal_time,decimation)
goal_decimated = decimation *round(goal_time/decimation);
end
