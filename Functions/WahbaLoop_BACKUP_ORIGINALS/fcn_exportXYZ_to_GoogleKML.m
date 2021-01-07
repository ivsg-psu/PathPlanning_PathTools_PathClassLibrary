function fcn_exportXYZ_to_GoogleKML(data,filename)
% Revision history:
% 2019_11_21 - first write of the code
% data: struct

%% Start by checking to see if the appropriate fields exist

% Initialize flags
flag_use_LLA = 1;
flag_use_ENU = 1;


% Check to see if we should use LLA
if ~isfield(data,'Latitude') || ~isfield(data,'Longitude') || ~isfield(data,'Altitude')
    flag_use_LLA = 0;
else
    flag_use_ENU = 0;
end

% If can't use LLA, check if we can use ENU
 if 0==flag_use_LLA
    if ~isfield(data,'xEast') || ~isfield(data,'yNorth') || ~isfield(data,'zUp')
        flag_use_ENU = 0;
    end    
 end

 % If using ENU, need to convert it

if 1==flag_use_LLA 
     LAT = data.Latitude;
     LON = data.Longitude;
     ALT = data.Altitude;
elseif 1==flag_use_ENU
    lat0 =40+48/60+ 24.81098/3600; %cenvert to degree units 
    lon0 = -77 - 50/60 - 59.26859/3600;
    h0 = 337.6654968261719; %
    spheroid = referenceEllipsoid('wgs84');
    [LAT,LON,ALT] = enu2geodetic(data.xEast,data.yNorth,data.zUp,lat0,lon0,h0,spheroid);
else 
     error('No available ENU or LLA data! \n')
end
 %% Convert results out

% kmlwrite(FILENAME, LAT, LON, ALT)
kmlwriteline(filename,LAT,LON,ALT, ...
       'Description', 'this is a description',...
       'Name', filename,...
       'Color','r',...
       'Width',1,...
       'AltitudeMode','clampToGround');
%
% filename2 = 'DGPS_2019_09_17_WahbaLoop_StartTerminal_clampToGround.kml';
% name2 = 'right_GPS';
% kmlwriteline(filename2,GPS_right(start_index:end_index,1), GPS_right(start_index:end_index,2), GPS_right(start_index:end_index,3),...
%     'Name',name2,'Color','b','Width',4, ...
%     'AltitudeMode','clampToGround');% relativeToSeaLevel,clampToGround

end

