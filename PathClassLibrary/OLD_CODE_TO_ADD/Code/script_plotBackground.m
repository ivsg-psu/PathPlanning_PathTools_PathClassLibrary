figure
xlim = [min(x) max(x)];
ylim = [max(y) min(y)];
hold on;
plot(loops(1).MergedGPS.xEast, loops(1).MergedGPS.yNorth, 'r')
length(lapData{1}.MergedGPS.xEast)

plot(loops(2).MergedGPS.xEast, loops(2).MergedGPS.yNorth, 'g')
plot(loops(3).MergedGPS.xEast, loops(3).MergedGPS.yNorth, 'm')
plot(loops(4).MergedGPS.xEast, loops(4).MergedGPS.yNorth, 'y')
plot(loops(5).MergedGPS.xEast, loops(5).MergedGPS.yNorth, 'w')
plot(rawData.GPS_Hemisphere.xEast(91500:100000), rawData.GPS_Hemisphere.yNorth(91500:100000), 'b')

legend('loop 1', 'loop 2', 'loop 3', 'loop 4', 'loop 5', 'lane changes', 'location', 'southeast')
axis tight; 

I = imread('test_track_map.png'); 
h = image(xlim,ylim,I); 
uistack(h,'bottom')

xlabel("x [m]")
ylabel("y [m]")