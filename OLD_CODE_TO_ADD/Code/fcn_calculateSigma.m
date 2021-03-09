function sigma = fcn_calculateSigma(closestXs,closestYs,mean_x, mean_y)

%fcn_lateralOffset Summary of this function goes here
% input:  aligned_Data_ByStation: laps data
%             mean_Data: the average data of route    
%               numLaps: the number of laps
% output: the maginitude and direction(left or right) of lateral offfset of
% laps data relative to average data.
%   Detailed explanation goes here
% update: 10/21, calculate the direction using satya's function 
% update: 10/23, calculate the magnitude using projection 

% issues: the maginitude calculated by just corresponding index has errors

%--------------------------------------------------------------------------------------------------------------------------
numLaps = size(closestXs,2);
for i_Laps=1:numLaps
            
        point_difference = [closestXs(:,i_Laps); closestYs(:, i_Laps)]-[mean_x; mean_y];
           
        lateral_offset_magnitude(i_Laps,:)=sqrt(point_difference(1,:).^2+point_difference(2,:).^2);  %calculate the magintude 
       
        %calculate the offset magintude 
        east_gps_interp = closestXs(:,i_Laps);
        north_gps_interp = closestYs(:, i_Laps);
        
        targetENUArray = [mean_x  mean_y zeros(length(mean_x),1)];
        queryArray_mag = [east_gps_interp, north_gps_interp, zeros(length(mean_x), 1)];
        
        [lateral_offset_magnitude(i_Laps,:)] = fcn_latealOffsetMagnitude(targetENUArray,queryArray_mag);
   
        
        % calculate the direction(left or right) lateral offset
        ENU = [mean_x  mean_y  zeros(length(mean_x),1)];
        queryArray = [east_gps_interp north_gps_interp zeros(length(mean_x),1)];
        
        lateral_offset_direction(i_Laps,:) = fcn_Locate_Points_Path(queryArray, ENU);
        
        % combine the magnitude and direction 
        lateral_offset(i_Laps,:)= lateral_offset_magnitude(i_Laps,:).*lateral_offset_direction(i_Laps,:);
         
         
        %h_fig = figure (200);
        %set(h_fig,'Name','lateral_offset_magnitude');
        %hold on 
        %plot(aligned_Data_ByStation.traversal{i_Laps}.station,lateral_offset_magnitude(i_Laps,:))
        %legend_string_offset{i_Laps} = sprintf('Lap %d',i_Laps);
           
end


sigma = (std(lateral_offset))';
%figure
%plot(aligned_Data_ByStation.traversal{i_Laps}.station, sigma,'r','LineWidth',2)
%hold on
%plot(aligned_Data_ByStation.traversal{i_Laps}.station, lateral_offset_mean,'r','LineWidth',2)
end

