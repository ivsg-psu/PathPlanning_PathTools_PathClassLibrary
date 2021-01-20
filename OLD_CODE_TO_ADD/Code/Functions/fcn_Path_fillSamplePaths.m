function paths = fcn_Path_fillSamplePaths
% fcn_Path_fillSamplePaths
% Produces dummy sample paths. Note: can go into the function and change
% flag to allow user-selected paths.
%
% FORMAT:
%
%       paths = fcn_Path_fillSamplePaths
%
% INPUTS:
%
%      (none)
%
% OUTPUTS:
%
%      X: an N x 1 vector of X positions
%      Y: an N x 1 vector of Y positions
%
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_Path_fillSamplePaths.m for a full
%       test suite.
%
% This function was written on 2020_11_12 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
%      2020_11_12 
%      -- wrote the code
%      2021_01_07 
%      -- minor updates to comments


flag_do_debug = 0; % Flag to plot the results for debugging
flag_grab_user_inputs = 0; % Flag to allow user to click to draw paths
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
end


%% check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag_check_inputs == 1
    % Are there the right number of inputs?
    if  nargin > 0
        error('Incorrect number of input arguments')
    end
    
end

%% Solve for the circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Build or load Path1
if 1==flag_grab_user_inputs
    % Grab some points manually to create a starter path
    figure(1);
    axis([0 100 0 100]);
    
    % Get a set of starting points
    [X1,Y1] = ginput;
    paths{1} = [X1 Y1];
else
    paths{1} = [
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
end


% Build or load Path 2
if 1==flag_grab_user_inputs
    % Grab other points manually to create a test path
    figure(1);

    % Show prior results
    clf; hold on;
    for i_Path = 1:length(paths)
        plot(paths{i_Path}(:,1),paths{i_Path}(:,2),'-');
        text(paths{i_Path}(1,1),paths{i_Path}(1,2),'Start');
    end
    
    axis([0 100 0 100]);

    % Get a set of starting points
    [X2,Y2] = ginput;
    plot(X2,Y2)
    paths{2} = [X2 Y2];
else
    paths{2} = [
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

% Build or load Path 3
if 1==flag_grab_user_inputs
    % Grab other points manually to create a test path
    figure(1);

    % Show prior results
    clf; hold on;
    for i_Path = 1:length(paths)
        plot(paths{i_Path}(:,1),paths{i_Path}(:,2),'-');
        text(paths{i_Path}(1,1),paths{i_Path}(1,2),'Start');
    end
   
    axis([0 100 0 100]);
    

    % Get a set of starting points
    [X3,Y3] = ginput;
    plot(X3,Y3)
    paths{3} = [X3 Y3];
else
    paths{3} = [
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


%% Any debugging?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _
%  |  __ \     | |
%  | |  | | ___| |__  _   _  __ _
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_do_debug
    % Prep a figure location
    close all
    figure(1);
    clf;
    hold on;
    grid minor;
    
    % Show result
    for i_Path = 1:length(paths)
        plot(paths{i_Path}(:,1),paths{i_Path}(:,2),'-');
        text(paths{i_Path}(1,1),paths{i_Path}(1,2),'Start');
    end

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end
end
