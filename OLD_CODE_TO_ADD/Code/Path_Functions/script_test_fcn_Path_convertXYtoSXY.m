% script_test_fcn_Path_convertXYtoSXY.m
% Tests fcn_Path_convertXYtoSXY

% Revision history
%      2020_11_10
%      -- first write of the code
%      2021_01_07
%      -- more comments (slightly)

% Fill in some dummy data
paths = fcn_Path_fillSamplePaths;

figure(1);
clf;
hold on;
grid on;
grid minor;

Station_step = 30;

% Show how the function works and plot result
for i_traveral = 1:length(paths)
    pathSXY = fcn_Path_convertXYtoSXY(paths{i_traveral}(:,1),paths{i_traveral}(:,2));
    plot(pathSXY(:,2),pathSXY(:,3),'.-','Linewidth',3,'Markersize',20)
end

% Plot station markers
for i_traveral = 1:length(paths)
    pathSXY = fcn_Path_convertXYtoSXY(paths{i_traveral}(:,1),paths{i_traveral}(:,2));
    for i_station = Station_step:Station_step:pathSXY(end,1)
        index = find(pathSXY(:,1)>=i_station,1);
        plot(pathSXY(index,2),pathSXY(index,3),'k.','Markersize',15);
        text(pathSXY(index,2),pathSXY(index,3),sprintf('Station: %.2f',pathSXY(index,1)));
    end
end



