function [start_point] = fcn_defineStartPoint(route_name,timeFilteredData)
%FCN_DEFINESTARTPOINT Summary of this function goes here
%   Detailed explanation goes here
     %plot and define the start point by observing 
     start_index= 2400;
     %[start_longitude,start_latitude, start_xEast,start_yNorth] = fcn_plotraw_lla(timeFilteredData.GPS_Hemisphere, 1254, 'lla_raw_data', start_index);
     %fcn_googleEarth('test_track',timeFilteredData.GPS_Hemisphere)
     
  if route_name== 1 % 1 means 'test_track';  2 means wahba_loop; 
        start_point.start_longitude=-77.833842140800000;  %deg
        start_point.start_latitude =40.862636161300000;   %deg
        start_point.start_xEast=1345.204537286125; % meters
        start_point.start_yNorth=6190.884280063217; % meters
        
        start_point.end_longitude=-77.833842140800000;  %deg
        start_point.end_latitude =40.862636161300000;   %deg
        start_point.end_xEast=1345.204537286125; % meters
        start_point.end_yNorth=6190.884280063217; % meters
        
        start_point.start_yaw_angle = 37.38; %deg
        start_point.expectedRouteLength = 1555.5; % meters
        start_point.direction = CCW; % meters
  elseif route_name == 2
        start_point.start_longitude=-77.87652037222755;
        start_point.start_latitude =40.828390558947870;
        start_point.start_xEast=-2254.319012077573;
        start_point.start_yNorth=2387.887818394200;
        
        start_point.end_longitude=-77.87652037222755;
        start_point.end_latitude =40.828390558947870;
        start_point.end_xEast=-2254.319012077573;
        start_point.end_yNorth=2387.887818394200;
        
        start_point.start_yaw_angle = 190;
        start_point.expectedRouteLength = 11265.5;
        start_point.direction = CCW; % meters
       
  end
  
end

