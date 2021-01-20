function [RouteStructure] = fcn_defineRouteStructure(route_name)
    
if strcmpi(route_name,'loop_1')
    RouteStructure.start_longitude= -77.836323768500000;  %deg
    RouteStructure.start_latitude = 40.863483410000000;   %deg
    RouteStructure.start_xEast = 1.135953863897120e+03; % meters
    RouteStructure.start_yNorth = 6.284941486180607e+03; % meters
     
    RouteStructure.end_longitude= -77.836323768500000;  %deg
    RouteStructure.end_latitude = 40.863483410000000;   %deg
    RouteStructure.end_xEast = 1.135953863897120e+03; % meters
    RouteStructure.end_yNorth = 6.284941486180607e+03; % meters
    
    %{
    RouteStructure.end_longitude= -77.836323768500000;  %deg
    RouteStructure.end_latitude = 40.863483410000000;   %deg
    RouteStructure.end_xEast = 1.135953863897120e+03; % meters
    RouteStructure.end_yNorth = 6.284941486180607e+03; % meters
    %}
elseif strcmpi(route_name,'loop_2')
    RouteStructure.start_longitude= -77.836040952100000;  %deg
    RouteStructure.start_latitude = 40.864000106900000;   %deg
    RouteStructure.start_xEast = 1.159790178414793e+03; % meters
    RouteStructure.start_yNorth = 6.342329343260531e+03; % meters
    
    
    RouteStructure.end_longitude= -77.836483314500000;  %deg
    RouteStructure.end_latitude = 40.863850037700000;   %deg
    RouteStructure.end_xEast = 1.122495845967782e+03; % meters
    RouteStructure.end_yNorth = 6.325655854464573e+03; % meters
    
elseif strcmpi(route_name,'loop_3')
    RouteStructure.start_longitude= -77.836413936400000;  %deg
    RouteStructure.start_latitude = 40.863765895400000;   %deg
    RouteStructure.start_xEast = 1.128346757673815e+03; % meters
    RouteStructure.start_yNorth = 6.316312386800244e+03; % meters
        
    RouteStructure.end_longitude= 77.836393631400000;  %deg
    RouteStructure.end_latitude = 40.863750265000000;   %deg
    RouteStructure.end_xEast = 1.130059007539598e+03; % meters
    RouteStructure.end_yNorth = 6.314576875127237e+03; % meters
    
elseif strcmpi(route_name,'loop_4')
    RouteStructure.start_longitude= -77.836465213100000;  %deg
    RouteStructure.start_latitude = 40.864258590000000;   %deg
    RouteStructure.start_xEast = 1.124015034803460e+03; % meters
    RouteStructure.start_yNorth = 6.371028266226611e+03; % meters
    
    RouteStructure.end_longitude= -77.836465213100000;  %deg
    RouteStructure.end_latitude = 40.864258590000000;   %deg
    RouteStructure.end_xEast = 1.124015034803460e+03; % meters
    RouteStructure.end_yNorth = 6.371028266226611e+03; % meters
%{
    RouteStructure.end_longitude= -77.836836310600000;  %deg
    RouteStructure.end_latitude = 40.863853863400000;   %deg
    RouteStructure.end_xEast = 1.092733612507107e+03; % meters
    RouteStructure.end_yNorth = 6.487610698927554e+03; % meters
    %}
elseif strcmpi(route_name,'loop_5')
    RouteStructure.start_longitude= -77.836349935900000;  %deg
    RouteStructure.start_latitude = 40.863459680000000;   %deg
    RouteStructure.start_xEast = 1.134362081897091e+03; % meters
    RouteStructure.start_yNorth = 6.282305881113656e+03; % meters

    RouteStructure.end_longitude= -77.836324165200000;  %deg
    RouteStructure.end_latitude = 40.863564240900000;   %deg
    RouteStructure.end_xEast = 1.135919087827313e+03; % meters
    RouteStructure.end_yNorth = 6.293918612868068e+03; % meters

elseif strcmpi(route_name,'lane_change_full_loop')
    RouteStructure.start_longitude= -77.834870624600000;  %deg
    RouteStructure.start_latitude = 40.864918977200000;   %deg
    RouteStructure.start_xEast = 1.258446196989194e+03; % meters
    RouteStructure.start_yNorth = 6.444394489180419e+03; % meters

    RouteStructure.end_longitude= -77.834870624600000;  %deg
    RouteStructure.end_latitude =  40.864918977200000;   %deg
    RouteStructure.end_xEast = 1.258446196989194e+03; % meters
    RouteStructure.end_yNorth = 6.444394489180419e+03; % meters
    
elseif strcmpi(route_name,'lane_change')  % Added by SNB on 2020_05_21
    RouteStructure.start_longitude= -77.834870624600000;  %deg
    RouteStructure.start_latitude = 40.864918977200000;   %deg
    RouteStructure.start_xEast = 1190; % meters
    RouteStructure.start_yNorth = 6455; % meters

    RouteStructure.end_longitude= -77.834870624600000;  %deg
    RouteStructure.end_latitude =  40.864918977200000;   %deg
    RouteStructure.end_xEast = 1390; % meters
    RouteStructure.end_yNorth = 6512; % meters

end

