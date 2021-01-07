% script_test_fcn_Path_FindOrthogonalHitFromPathToPath.m
% This is a script to exercise the function: fcn_Path_FindOrthogonalHitFromPathToPath.m
% This function was written on 2020_12_30 by S. Brennan
%     Modified on 2020_12_30 to write original form of code
% Questions or comments? sbrennan@psu.edu


%% Real path examples
close all;

% Fill in some dummy data
paths = fcn_Path_fillSamplePaths;

% Convert paths to traversal structures
clear data;
for i_Path = 1:1 %length(paths)
    traversal = fcn_Path_convertXYtoTraversalStructure(paths{i_Path}(:,1),paths{i_Path}(:,2));
    data.traversal{i_Path} = traversal;
end

% Call the plot command to show results in XY
fcn_Path_plotPathXY(data,fig_num);

% Calculate the distances
traversal_to_check = data.traversal{1};

% Fill in the traversal information
stations = traversal_to_check.Station;
X = traversal_to_check.X;
Y = traversal_to_check.Y;
points = [X Y];

origin = [X(1) Y(1)];
distances = sum((points-origin).^2,2).^0.5;
differences = [0; diff(distances)];
sign_changes = differences(1:end-1).*differences(2:end);
extrema_indices = find(sign_changes<0);

figure(46464);
clf;
hold on;
grid on;
grid minor;

plot(stations,distances,'k-','Linewidth',3);
plot(stations(extrema_indices,1),distances(extrema_indices,1),'ro','Markersize',20)

figure(fig_num);
hold on;
plot(X(extrema_indices,1),Y(extrema_indices,1),'ro','Markersize',20)


function print_results(stations,closest_path_point,distances)
fprintf(1,'\n\nStation \t Location X \t Location Y \t Distance \n');
for i_station =1:length(stations)
    fprintf(1,'%.2f \t\t %.2f \t\t\t %.2f \t\t\t %.2f\n',stations(i_station),closest_path_point(i_station,1),closest_path_point(i_station,2),distances(i_station));
end
end
