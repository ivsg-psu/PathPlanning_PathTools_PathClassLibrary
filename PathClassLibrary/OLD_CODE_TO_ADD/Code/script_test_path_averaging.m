% script_test_path_averaging

% Steps:
% 1. Create some dummy test paths
% 2. Define a way of searching out projection points from reference path(e.g. first travesal) to other travesals, but avoiding cross-overs
% 3. Average these projection points to generate an initial average path
% 4. Project from the initial average path to nearby trajectories to find projections
% 5. Average projections to find new average path. Repeat until average path is not changing
% 6. Use final average path to define "true" s-coordinates of the original
% trajectories, using projection


% Steps:
% 1. Create some dummy test paths

clear all
close all
figure(1);
clf;
hold on;
grid minor;

if ~exist('path1')
    
    % Build or load Path1
    if 1==0
        % Grab some points manually to create a starter path
        figure(1);
        axis([0 100 0 100]);
        
        % Get a set of starting points
        [X1,Y1] = ginput;
        path1 = [X1 Y1];
    else
        path1 = [
            8.1797    8.6006
            8.4101   15.8892
            10.0230   25.5102
            10.4839   39.7959
            10.7143   50.8746
            10.9447   61.9534
            13.0184   73.6152
            16.7051   82.9446
            24.7696   89.3586
            35.3687   91.1079
            42.5115   91.3994
            52.8802   89.3586
            58.1797   85.2770
            60.9447   78.5714
            57.9493   72.1574
            57.2581   63.7026
            59.1014   58.1633
            63.0184   57.2886
            67.3963   56.9971
            69.9309   56.7055
            74.3088   56.1224
            78.6866   54.0816
            80.9908   51.4577
            82.3733   49.1254
            84.6774   40.6706
            84.6774   34.5481
            82.8341   28.7172
            80.0691   26.9679
            76.3825   25.2187
            69.2396   20.2624
            65.5530   18.2216
            60.9447   18.8047
            57.4885   22.0117
            50.5760   28.1341
            47.1198   30.7580
            43.8940   34.8397
            39.7465   37.7551
            37.6728   40.9621
            32.6037   42.7114
            30.0691   43.0029
            28.2258   43.2945
            26.3825   43.2945
            24.5392   44.4606
            20.8525   47.3761
            19.2396   49.7085
            16.7051   53.2070
            14.1705   58.1633
            13.0184   62.2449
            10.0230   70.1166
            8.6406   74.4898
            7.7189   79.7376
            6.5668   82.9446
            5.1843   86.7347
            4.2627   88.4840
            3.8018   89.0671];
        X1 = path1(:,1);
        Y1 = path1(:,2);
    end
    % Show result
    figure(1); plot(path1(:,1),path1(:,2),'r-o');
    text(path1(1,1),path1(1,2),'Start');
    % Build or load Path 2
    if 1==0
        % Grab other points manually to create a test path
        figure(1);
        axis([0 100 0 100]);
        
        % Get a set of starting points
        [X2,Y2] = ginput;
        plot(X2,Y2)
        path2 = [X2 Y2];
    else
        path2 = [
            9.1014    9.7668
            8.6406   17.0554
            8.6406   26.0933
            11.4055   35.7143
            10.9447   44.1691
            10.7143   49.7085
            10.7143   53.4985
            10.7143   59.6210
            11.4055   65.1603
            12.7880   71.5743
            14.1705   78.5714
            19.0092   85.8601
            24.3088   88.4840
            31.6820   88.7755
            40.8986   90.2332
            51.4977   91.3994
            56.3364   87.9009
            63.4793   82.3615
            61.8664   77.4052
            59.3318   71.5743
            58.1797   65.7434
            57.7189   62.2449
            59.7926   54.9563
            63.4793   54.9563
            69.2396   57.2886
            76.1521   55.2478
            81.4516   50.5831
            82.3733   49.4169
            84.9078   45.0437
            85.8295   39.7959
            85.3687   36.0058
            84.2166   31.3411
            82.1429   25.5102
            75.0000   23.1778
            71.7742   21.7201
            65.0922   19.0962
            62.7880   18.8047
            57.7189   20.8455
            56.1060   22.8863
            54.2627   25.5102
            50.8065   29.5918
            48.0415   30.4665
            43.4332   31.3411
            40.6682   35.1312
            38.3641   38.9213
            36.7512   40.9621
            34.2166   43.2945
            31.2212   43.5860
            26.8433   44.1691
            24.0783   45.3353
            21.7742   47.3761
            17.8571   51.4577
            15.0922   55.5394
            12.5576   59.9125
            11.1751   63.9942
            10.4839   68.3673
            10.9447   76.2391
            11.1751   84.9854
            8.6406   90.8163
            5.4147   94.6064];
    end
    % Show result
    figure(1); plot(path2(:,1),path2(:,2),'g-o');
    
    % Build or load Path 2
    if 1==0
        % Grab other points manually to create a test path
        figure(1);
        axis([0 100 0 100]);
        
        % Get a set of starting points
        [X3,Y3] = ginput;
        plot(X3,Y3)
        path3 = [X3 Y3];
    else
        path3 = [            
        9.7926   10.6414
        9.7926   15.8892
        9.5622   19.0962
        9.3318   19.9708
        8.8710   21.7201
        8.8710   22.8863
        9.3318   24.0525
        10.2535   25.2187
        10.9447   26.6764
        10.9447   26.9679
        11.4055   28.4257
        11.4055   29.8834
        11.4055   31.0496
        11.1751   31.6327
        10.9447   33.6735
        10.2535   35.7143
        10.2535   36.8805
        10.2535   38.3382
        10.0230   40.3790
        10.0230   41.5452
        10.0230   43.2945
        10.2535   45.9184
        10.4839   49.4169
        11.6359   52.9155
        12.0968   55.5394
        12.3272   58.4548
        12.3272   60.4956
        12.5576   64.2857
        13.4793   69.2420
        13.7097   71.8659
        13.7097   74.4898
        14.1705   76.8222
        14.1705   77.9883
        14.4009   79.7376
        14.6313   82.3615
        15.3226   83.8192
        16.2442   85.2770
        18.7788   88.1924
        19.9309   89.3586
        21.7742   89.9417
        26.6129   89.9417
        29.6083   90.2332
        31.9124   90.5248
        36.7512   91.1079
        37.9032   91.1079
        42.5115   91.1079
        45.7373   91.1079
        49.1935   91.6910
        52.6498   91.6910
        56.7972   90.5248
        62.7880   90.5248
        64.1705   89.6501
        66.9355   85.5685
        68.0876   82.6531
        67.6267   80.3207
        66.0138   77.1137
        63.7097   75.0729
        62.5576   73.6152
        60.0230   69.5335
        57.4885   65.4519
        55.8756   62.8280
        55.4147   59.6210
        55.6452   58.4548
        56.1060   57.5802
        57.9493   56.1224
        59.5622   55.5394
        62.3272   56.4140
        65.3226   57.8717
        71.3134   57.8717
        75.6912   57.5802
        79.1475   56.4140
        80.7604   55.8309
        82.3733   54.9563
        85.8295   45.6268
        87.4424   44.1691
        87.6728   40.6706
        86.5207   39.5044
        85.5991   32.2157
        83.9862   28.7172
        83.0645   26.6764
        80.9908   24.0525
        76.8433   21.7201
        73.6175   20.5539
        71.0829   19.0962
        66.7051   18.8047
        63.4793   19.6793
        62.3272   20.8455
        59.7926   22.8863
        58.8710   24.6356
        55.6452   26.9679
        54.9539   28.7172
        52.1889   31.9242
        50.8065   33.3819
        50.1152   33.6735
        47.5806   34.5481
        43.6636   34.5481
        37.2120   35.1312
        34.4470   35.1312
        33.7558   37.4636
        32.6037   41.2536
        31.6820   43.8776
        24.0783   47.3761
        20.8525   47.9592
        13.2488   48.5423
        12.5576   49.1254
        11.8664   52.9155
        12.3272   57.2886
        12.3272   59.9125
        12.5576   62.2449
        12.5576   65.4519
        12.0968   69.2420
        10.7143   75.6560
        9.3318   78.5714
        8.8710   82.6531
        8.6406   84.4023
        7.9493   86.4431
        6.5668   88.4840
        5.8756   89.6501];
    end
    % Show result
    figure(1); plot(path3(:,1),path3(:,2),'b-o');
end

data.traversal{1}.X = path1(:,1);
data.traversal{1}.Y = path1(:,2);
data.traversal{1}.Z = 0*path1(:,1);
data.traversal{1}.Diff = [[0 0]; diff(path1)];
data.traversal{1}.Station = cumsum(sqrt(sum(data.traversal{1}.Diff.^2,2)));
data.traversal{1}.Yaw = fnc_yaw(path1(:,1),path1(:,2));

data.traversal{2}.X = path2(:,1);
data.traversal{2}.Y = path2(:,2);
data.traversal{2}.Z = 0*path2(:,1);
data.traversal{2}.Diff = [[0 0]; diff(path2)];
data.traversal{2}.Station = cumsum(sqrt(sum(data.traversal{2}.Diff.^2,2)));
data.traversal{2}.Yaw = fnc_yaw(path2(:,1),path2(:,2));

data.traversal{3}.X = path3(:,1);
data.traversal{3}.Y = path3(:,2);
data.traversal{3}.Z = 0*path3(:,1);
data.traversal{3}.Diff = [[0 0]; diff(path3)];
data.traversal{3}.Station = cumsum(sqrt(sum(data.traversal{3}.Diff.^2,2)));
data.traversal{3}.Yaw = fnc_yaw(path3(:,1),path3(:,2));

figure(11)
hold on
for i_path= 1:length(data.traversal)
    plot(data.traversal{i_path}.Station,data.traversal{i_path}.Yaw)
end
title('station vs yaw')
xlabel('station [m]')
ylabel('yaw [degree]')

%% 2. Define a way of searching out projection points from reference path(e.g. first travesal) to other travesals, but avoiding cross-overs


figure(2)
hold on
for i_path= 1:length(data.traversal)
    plot(data.traversal{i_path}.X,data.traversal{i_path}.Y,'-o')
end
title('path')
xlabel('x [m]')
ylabel('y [m]')

% choose the initial reference path, choose the path with most data points 
data_length = zeros(1,length(data.traversal));
for i_path = 1:length(data.traversal)
    data_length(i_path) = length(data.traversal{i_path}.X);
end
[~,initital_reference_path_id] = max(data_length);
path = data.traversal{initital_reference_path_id}; %initial reference path


% path_last  = path; % store the last step data in the loop
num_iteration = 10;  % the number of iteration to fidn the average path 
iteration_error  = zeros(num_iteration,3);
for i =1:num_iteration 
    % the loop:
    % 4. Project from the average path to nearby trajectories to find projections
    % 5. Average projections to find new average path. Repeat until average
    % path is not changing

    % step2: searching out nerest points from reference path to other travesals, but avoiding cross-overs
    [closestXs,closestYs,closestZs,closestYaws] = findClosestPointsFromPath(path, data,1);
    
    % last path 
    path_last  = path;
    
    figure(22)
    hold on
    for i_path= 1:length(data.traversal)
        plot(data.traversal{i_path}.X,data.traversal{i_path}.Y,'-o')
        plot(closestXs(:,i_path),closestYs(:,i_path),'b-o')
    end
    title('original path and nearest')
    xlabel('x [m]')
    ylabel('y [m]')
    
    %%3. Average these projection points to generate an initial average path
    
    path_average = [mean(closestXs,2) mean(closestYs,2) mean(closestZs,2)];
    path_average_station  = [0; cumsum(sqrt(sum(diff(path_average).^2,2)))];
    path_average_yaw = mean(closestYaws,2);
    
    figure(2)
    hold on
    plot(path_average(:,1), path_average(:,2),'LineWidth',3)
    for i_path= 1:length(data.traversal)
        plot(data.traversal{i_path}.X,data.traversal{i_path}.Y)
    end
    title('original path and average')
    xlabel('x [m]')
    ylabel('y [m]')
    
    %%3.5 Interpolation of mean data by equal interval 
    interval = 1; % meters
    nb_points = round(path_average_station(end)/interval); % the number of points after interplation 
    % interplate the X,Y,Z
    path_average_interp = fcn_interparc(nb_points,path_average(:,1), path_average(:,2),path_average(:,3),'spline');
    % interplate the yaw
    path_average_interp_station = [0; cumsum(sqrt(sum(diff(path_average_interp).^2,2)))];
    path_average_interp_yaw = interp1(path_average_station, path_average_yaw, path_average_interp_station,'spline','extrap');
         
    % update path 
    path.X = path_average_interp(:,1);
    path.Y = path_average_interp(:,2);
    path.Z = path_average_interp(:,3);
    path.Yaw = path_average_interp_yaw;
    path.Station = path_average_interp_station;
    
    % error, Notes: more eidt to make sure their data length are the same.
    if i>=2
%         iteration_error(i) = mean(path_average - path_average_last);
    end

    %figure(4)
    hold on
    plot(path_average_interp(:,1), path_average_interp(:,2),'.','LineWidth',1)
    
    pause(0.1)
end

%% 6. Use final average path to define "true" s-coordinates of the original trajectories, using projection
path_average_final  = path;

% plot the result
figure(111)
hold on
plot(path.Station, path.Yaw,'r-' ,'LineWidth',2)
plot(path.Station,fnc_yaw(path.X,path.Y),'b-' ,'LineWidth',2)
for i_path= 1:length(data.traversal)
    plot(data.traversal{i_path}.Station,data.traversal{i_path}.Yaw)
end
legend('result of interplation', 'result of final calculation')
title('station vs yaw')
xlabel('station [m]')
ylabel('yaw [degree]')

figure(222)
clf
hold on
for i_path= 1:length(data.traversal)
    %plot(data.traversal{i_path}.X,data.traversal{i_path}.Y,'-')
    plot(closestXs(:,i_path),closestYs(:,i_path),'-o')
end
plot(path_last.X,path_last.Y,'r-o','LineWidth',2)
plot(path_average(:,1),path_average(:,2),'b-o','LineWidth',2)

title('original, closest points and mean path')
xlabel('x [m]')
ylabel('y [m]')


% ===============================================
%  nested function for yaw angle calculation (NOT needed for real data)
% ===============================================
function lane_yaw = fnc_yaw(X,Y)
% Yaw,north is zero, clockwise is positive direction, range 0-360 degrees

% Make sure that x and y are column vectors.
X=X(:);
Y=Y(:);
% calculate the atan2d and convert it to Yaw
yaw = 90 - atan2d(diff(Y), diff(X)); %
lane_yaw = [yaw(1); yaw];
n_yaw  = lane_yaw < 0;
lane_yaw(n_yaw) = lane_yaw(n_yaw)+360;
end