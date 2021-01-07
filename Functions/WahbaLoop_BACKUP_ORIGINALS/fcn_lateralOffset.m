function [lateral_offset,lateral_offset_direction,lateral_offset_magnitude] = fcn_lateralOffset(aligned_Data_ByStation,mean_Data,numLaps)
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
    
for i_Laps=1:numLaps
            
        point_difference=[aligned_Data_ByStation.traversal{i_Laps}.X ;aligned_Data_ByStation.traversal{i_Laps}.Y]...
            -[mean_Data.mean_xEast'; mean_Data.mean_yNorth'];
           
        lateral_offset_magnitude(i_Laps,:)=sqrt(point_difference(1,:).^2+point_difference(2,:).^2);  %calculate the magintude 
       
        %calculate the offset magintude 
        east_gps_interp = aligned_Data_ByStation.traversal{i_Laps}.X;
        north_gps_interp = aligned_Data_ByStation.traversal{i_Laps}.Y;
        
        targetENUArray = [mean_Data.mean_xEast  mean_Data.mean_yNorth zeros(length(mean_Data.mean_xEast),1)];
        queryArray_mag = [east_gps_interp' north_gps_interp' zeros(length(mean_Data.mean_xEast),1)];
        
        [lateral_offset_magnitude(i_Laps,:)] = fcn_latealOffsetMagnitude(targetENUArray,queryArray_mag);
   
        
        % calculate the direction(left or right) lateral offset
        ENU = [mean_Data.mean_xEast  mean_Data.mean_yNorth zeros(length(mean_Data.mean_xEast),1)];
        queryArray = [east_gps_interp' north_gps_interp' zeros(length(mean_Data.mean_xEast),1)];
        
        lateral_offset_direction(i_Laps,:) = fcn_Locate_Points_Path(queryArray, ENU);
        
        % combine the magnitude and direction 
        lateral_offset(i_Laps,:)= lateral_offset_magnitude(i_Laps,:).*lateral_offset_direction(i_Laps,:);
         
         
        h_fig = figure (200);
        set(h_fig,'Name','lateral_offset_magnitude');
        hold on 
        plot(aligned_Data_ByStation.traversal{i_Laps}.station,lateral_offset_magnitude(i_Laps,:))
        legend_string_offset{i_Laps} = sprintf('Lap %d',i_Laps);
           
end

hold on 
lateral_offset_mean_magnitude=mean(lateral_offset_magnitude);
plot(aligned_Data_ByStation.traversal{i_Laps}.station,lateral_offset_mean_magnitude,'r','LineWidth',2)
grid on
xlabel('Station Distance [m]')
ylabel('Lateral Offset magnitude [m]')
legend([legend_string_offset,'mean']);

%%
h_fig = figure (2032);
set(h_fig,'Name','lateral_offset_magnitude_direction');

legend_string = ''; % Initialize an empty string
for i_Laps = 1:numLaps
    plot(aligned_Data_ByStation.traversal{i_Laps}.station,lateral_offset(i_Laps,:),'LineWidth', 1);
    hold on
    legend_string{i_Laps} = sprintf('Lap %d',i_Laps);
end

lateral_offset_mean=mean(lateral_offset);
plot(aligned_Data_ByStation.traversal{i_Laps}.station,lateral_offset_mean,'r','LineWidth',2)
grid on
xlabel('Station Distance [m]')
ylabel('Lateral Offset [m]')
legend([legend_string_offset,'mean_route']);

end

