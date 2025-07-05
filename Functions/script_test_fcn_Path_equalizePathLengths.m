% script_test_fcn_Path_equalizePathLengths
% Tests the function: fcn_Path_equalizePathLengths

% Revision history
% 2025_06_26 - Sean Brennan
% -- first write of the code

close all;


%% BASIC CALL: Very simple example where first segment is reference
fig_num = 10001;
titleString = sprintf('BASIC CALL: Very simple example where first segment is reference');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

cellArrayOfUnequalPaths{1,1} = [1 4; 10 4];
cellArrayOfUnequalPaths{2,1} = [0 2; 7 2];
cellArrayOfUnequalPaths{3,1} = [3 0; 6 0];

[cellArrayOfEqualizedPaths, leastExtensionIndex, bestStartIndex, bestEndIndex] = fcn_Path_equalizePathLengths(cellArrayOfUnequalPaths,(fig_num));

title(titleString, 'Interpreter','none');

% Check variable types
assert(iscell(cellArrayOfEqualizedPaths));
assert(isnumeric(leastExtensionIndex));
assert(isnumeric(bestStartIndex));
assert(isnumeric(bestEndIndex));

% Check variable sizes
% Show that the cell array has same number of paths after extensions
assert(isequal(size(cellArrayOfEqualizedPaths),[length(cellArrayOfUnequalPaths) 1]));
% Show extensions make the cell array longer
for ith_cell = 1:length(cellArrayOfUnequalPaths)
    assert(length(cellArrayOfEqualizedPaths{ith_cell}(:,1))>=length(cellArrayOfUnequalPaths{ith_cell}(:,1)))
end
assert(isequal(size(leastExtensionIndex),[1 1]));
assert(isequal(size(bestStartIndex),[1 1]));
assert(isequal(size(bestEndIndex),[1 1]));

% Check variable values
assert(isequal(cellArrayOfEqualizedPaths{1},[0 4; cellArrayOfUnequalPaths{1}]));
assert(isequal(cellArrayOfEqualizedPaths{2},[cellArrayOfUnequalPaths{2}; 10 2]));
assert(isequal(cellArrayOfEqualizedPaths{3},[0 0; cellArrayOfUnequalPaths{3}; 10 0]));
assert(isequal(leastExtensionIndex,1));
assert(isequal(bestStartIndex,2));
assert(isequal(bestEndIndex,1));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC CALL: Very simple example where second segment is reference
fig_num = 10002;
titleString = sprintf('BASIC CALL: Very simple example where second segment is reference');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

cellArrayOfUnequalPaths{1,1} = [0 4; 7 4];
cellArrayOfUnequalPaths{2,1} = [1 2; 10 2];
cellArrayOfUnequalPaths{3,1} = [3 0; 6 0];

[cellArrayOfEqualizedPaths, leastExtensionIndex, bestStartIndex, bestEndIndex] = fcn_Path_equalizePathLengths(cellArrayOfUnequalPaths,(fig_num));

title(titleString, 'Interpreter','none');

% Check variable types
assert(iscell(cellArrayOfEqualizedPaths));
assert(isnumeric(leastExtensionIndex));
assert(isnumeric(bestStartIndex));
assert(isnumeric(bestEndIndex));

% Check variable sizes
% Show that the cell array has same number of paths after extensions
assert(isequal(size(cellArrayOfEqualizedPaths),[length(cellArrayOfUnequalPaths) 1]));
% Show extensions make the cell array longer
for ith_cell = 1:length(cellArrayOfUnequalPaths)
    assert(length(cellArrayOfEqualizedPaths{ith_cell}(:,1))>=length(cellArrayOfUnequalPaths{ith_cell}(:,1)))
end
assert(isequal(size(leastExtensionIndex),[1 1]));
assert(isequal(size(bestStartIndex),[1 1]));
assert(isequal(size(bestEndIndex),[1 1]));

% Check variable values
assert(isequal(cellArrayOfEqualizedPaths{1},[cellArrayOfUnequalPaths{1}; 10 4]));
assert(isequal(cellArrayOfEqualizedPaths{2},[0 2; cellArrayOfUnequalPaths{2}]));
assert(isequal(cellArrayOfEqualizedPaths{3},[0 0; cellArrayOfUnequalPaths{3}; 10 0]));
assert(isequal(leastExtensionIndex,2));
assert(isequal(bestStartIndex,1));
assert(isequal(bestEndIndex,2));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC CALL: Very simple example where third segment is reference
fig_num = 10003;
titleString = sprintf('BASIC CALL: Very simple example where third segment is reference');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

cellArrayOfUnequalPaths{1,1} = [0 4; 7 4];
cellArrayOfUnequalPaths{2,1} = [3 2; 6 2];
cellArrayOfUnequalPaths{3,1} = [1 0; 10 0];

[cellArrayOfEqualizedPaths, leastExtensionIndex, bestStartIndex, bestEndIndex] = fcn_Path_equalizePathLengths(cellArrayOfUnequalPaths,(fig_num));

title(titleString, 'Interpreter','none');

% Check variable types
assert(iscell(cellArrayOfEqualizedPaths));
assert(isnumeric(leastExtensionIndex));
assert(isnumeric(bestStartIndex));
assert(isnumeric(bestEndIndex));

% Check variable sizes
% Show that the cell array has same number of paths after extensions
assert(isequal(size(cellArrayOfEqualizedPaths),[length(cellArrayOfUnequalPaths) 1]));
% Show extensions make the cell array longer
for ith_cell = 1:length(cellArrayOfUnequalPaths)
    assert(length(cellArrayOfEqualizedPaths{ith_cell}(:,1))>=length(cellArrayOfUnequalPaths{ith_cell}(:,1)))
end
assert(isequal(size(leastExtensionIndex),[1 1]));
assert(isequal(size(bestStartIndex),[1 1]));
assert(isequal(size(bestEndIndex),[1 1]));

% Check variable values
assert(isequal(cellArrayOfEqualizedPaths{1},[cellArrayOfUnequalPaths{1}; 10 4]));
assert(isequal(cellArrayOfEqualizedPaths{2},[0 2; cellArrayOfUnequalPaths{2}; 10 2]));
assert(isequal(cellArrayOfEqualizedPaths{3},[0 0; cellArrayOfUnequalPaths{3}]));
assert(isequal(leastExtensionIndex,3));
assert(isequal(bestStartIndex,1));
assert(isequal(bestEndIndex,3));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC CALL: Angled paths
fig_num = 10004;
titleString = sprintf('BASIC CALL: Angled paths');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

cellArrayOfUnequalPaths{1,1} = [1 4; 10 6];
cellArrayOfUnequalPaths{2,1} = [0 2; 7 2];
cellArrayOfUnequalPaths{3,1} = [3 0; 6 -1];

[cellArrayOfEqualizedPaths, leastExtensionIndex, bestStartIndex, bestEndIndex] = fcn_Path_equalizePathLengths(cellArrayOfUnequalPaths,(fig_num));

title(titleString, 'Interpreter','none');

% Check variable types
assert(iscell(cellArrayOfEqualizedPaths));
assert(isnumeric(leastExtensionIndex));
assert(isnumeric(bestStartIndex));
assert(isnumeric(bestEndIndex));

% Check variable sizes
% Show that the cell array has same number of paths after extensions
assert(isequal(size(cellArrayOfEqualizedPaths),[length(cellArrayOfUnequalPaths) 1]));
% Show extensions make the cell array longer
for ith_cell = 1:length(cellArrayOfUnequalPaths)
    assert(length(cellArrayOfEqualizedPaths{ith_cell}(:,1))>=length(cellArrayOfUnequalPaths{ith_cell}(:,1)))
end
assert(isequal(size(leastExtensionIndex),[1 1]));
assert(isequal(size(bestStartIndex),[1 1]));
assert(isequal(size(bestEndIndex),[1 1]));

% Check variable values
% assert(isequal(cellArrayOfEqualizedPaths{1},[0 4; cellArrayOfUnequalPaths{1}]));
% assert(isequal(cellArrayOfEqualizedPaths{2},[cellArrayOfUnequalPaths{2}; 10 2]));
% assert(isequal(cellArrayOfEqualizedPaths{3},[0 0; cellArrayOfUnequalPaths{3}; 10 0]));
assert(isequal(leastExtensionIndex,1));
assert(isequal(bestStartIndex,2));
assert(isequal(bestEndIndex,1));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC CALL: Multisegmented, angled paths
fig_num = 10005;
titleString = sprintf('BASIC CALL: Multisegmented, angled paths');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

cellArrayOfUnequalPaths{1,1} = [-5 3; 1 4; 10 6];
cellArrayOfUnequalPaths{2,1} = [-6 5; 0 2; 7 2];
cellArrayOfUnequalPaths{3,1} = [-4 1; 3 0; 6 -1];

[cellArrayOfEqualizedPaths, leastExtensionIndex, bestStartIndex, bestEndIndex] = fcn_Path_equalizePathLengths(cellArrayOfUnequalPaths,(fig_num));

title(titleString, 'Interpreter','none');

% Check variable types
assert(iscell(cellArrayOfEqualizedPaths));
assert(isnumeric(leastExtensionIndex));
assert(isnumeric(bestStartIndex));
assert(isnumeric(bestEndIndex));

% Check variable sizes
% Show that the cell array has same number of paths after extensions
assert(isequal(size(cellArrayOfEqualizedPaths),[length(cellArrayOfUnequalPaths) 1]));
% Show extensions make the cell array longer
for ith_cell = 1:length(cellArrayOfUnequalPaths)
    assert(length(cellArrayOfEqualizedPaths{ith_cell}(:,1))>=length(cellArrayOfUnequalPaths{ith_cell}(:,1)))
end
assert(isequal(size(leastExtensionIndex),[1 1]));
assert(isequal(size(bestStartIndex),[1 1]));
assert(isequal(size(bestEndIndex),[1 1]));

% Check variable values
% assert(isequal(cellArrayOfEqualizedPaths{1},[0 4; cellArrayOfUnequalPaths{1}]));
% assert(isequal(cellArrayOfEqualizedPaths{2},[cellArrayOfUnequalPaths{2}; 10 2]));
% assert(isequal(cellArrayOfEqualizedPaths{3},[0 0; cellArrayOfUnequalPaths{3}; 10 0]));
assert(isequal(leastExtensionIndex,1));
assert(isequal(bestStartIndex,2));
assert(isequal(bestEndIndex,1));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC CALL: Best end is not same as least changed
fig_num = 10006;
titleString = sprintf('BASIC CALL: Best end is not same as least changed');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

cellArrayOfUnequalPaths{1,1} = [0 5; 0 10];
cellArrayOfUnequalPaths{2,1} = [2 0; 2 9];
cellArrayOfUnequalPaths{3,1} = [4 1; 4 7];

[cellArrayOfEqualizedPaths, leastExtensionIndex, bestStartIndex, bestEndIndex] = fcn_Path_equalizePathLengths(cellArrayOfUnequalPaths,(fig_num));

title(titleString, 'Interpreter','none');

% Check variable types
assert(iscell(cellArrayOfEqualizedPaths));
assert(isnumeric(leastExtensionIndex));
assert(isnumeric(bestStartIndex));
assert(isnumeric(bestEndIndex));

% Check variable sizes
% Show that the cell array has same number of paths after extensions
assert(isequal(size(cellArrayOfEqualizedPaths),[length(cellArrayOfUnequalPaths) 1]));
% Show extensions make the cell array longer
for ith_cell = 1:length(cellArrayOfUnequalPaths)
    assert(length(cellArrayOfEqualizedPaths{ith_cell}(:,1))>=length(cellArrayOfUnequalPaths{ith_cell}(:,1)))
end
assert(isequal(size(leastExtensionIndex),[1 1]));
assert(isequal(size(bestStartIndex),[1 1]));
assert(isequal(size(bestEndIndex),[1 1]));

% Check variable values
assert(isequal(cellArrayOfEqualizedPaths{1},[0 0; cellArrayOfUnequalPaths{1}]));
assert(isequal(cellArrayOfEqualizedPaths{2},[cellArrayOfUnequalPaths{2}; 2 10]));
assert(isequal(cellArrayOfEqualizedPaths{3},[4 0; cellArrayOfUnequalPaths{3}; 4 10]));
assert(isequal(leastExtensionIndex,2));
assert(isequal(bestStartIndex,2));
assert(isequal(bestEndIndex,1));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% BASIC CALL: No ends hit each other
fig_num = 10007;
titleString = sprintf('BASIC CALL: No ends hit each other');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

cellArrayOfUnequalPaths{1,1} = [0 5; 0 10; -0.5 10.9];
cellArrayOfUnequalPaths{2,1} = [2 0; 2 10; 2 11];
cellArrayOfUnequalPaths{3,1} = [4 1; 4 10; 4.5 10.9];

[cellArrayOfEqualizedPaths, leastExtensionIndex, bestStartIndex, bestEndIndex] = fcn_Path_equalizePathLengths(cellArrayOfUnequalPaths,(fig_num));

title(titleString, 'Interpreter','none');

% Check variable types
assert(iscell(cellArrayOfEqualizedPaths));
assert(isnumeric(leastExtensionIndex));
assert(isnumeric(bestStartIndex));
assert(isnumeric(bestEndIndex));

% Check variable sizes
% Show that the cell array has same number of paths after extensions
assert(isequal(size(cellArrayOfEqualizedPaths),[length(cellArrayOfUnequalPaths) 1]));
% Show extensions make the cell array longer
for ith_cell = 1:length(cellArrayOfUnequalPaths)
    assert(length(cellArrayOfEqualizedPaths{ith_cell}(:,1))>=length(cellArrayOfUnequalPaths{ith_cell}(:,1)))
end
assert(isequal(size(leastExtensionIndex),[1 1]));
assert(isequal(size(bestStartIndex),[1 1]));
assert(isequal(size(bestEndIndex),[1 1]));

% Check variable values
assert(isequal(round(cellArrayOfEqualizedPaths{1},4),round([0 0; cellArrayOfUnequalPaths{1}(1:end,:); -0.5556   11.000],4)));
assert(isequal(round(cellArrayOfEqualizedPaths{2},4),round([cellArrayOfUnequalPaths{2}],4)));
assert(isequal(round(cellArrayOfEqualizedPaths{3},4),round([4 0; cellArrayOfUnequalPaths{3}; 4.5556   11.0000],4)));
assert(isequal(leastExtensionIndex,2));
assert(isequal(bestStartIndex,2));
assert(isequal(bestEndIndex,2));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% BASIC CALL: All ends hit each other
fig_num = 10008;
titleString = sprintf('BASIC CALL: All ends hit each other');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

cellArrayOfUnequalPaths{1,1} = [0 5; 0 10; 0.5 11.1];
cellArrayOfUnequalPaths{2,1} = [2 0; 2 10; 2 11];
cellArrayOfUnequalPaths{3,1} = [4 1; 4 10; 3.5 11.1];

[cellArrayOfEqualizedPaths, leastExtensionIndex, bestStartIndex, bestEndIndex] = fcn_Path_equalizePathLengths(cellArrayOfUnequalPaths,(fig_num));

title(titleString, 'Interpreter','none');

% Check variable types
assert(iscell(cellArrayOfEqualizedPaths));
assert(isnumeric(leastExtensionIndex));
assert(isnumeric(bestStartIndex));
assert(isnumeric(bestEndIndex));

% Check variable sizes
% Show that the cell array has same number of paths after extensions
assert(isequal(size(cellArrayOfEqualizedPaths),[length(cellArrayOfUnequalPaths) 1]));
% Show extensions make the cell array longer
for ith_cell = 1:length(cellArrayOfUnequalPaths)
    assert(length(cellArrayOfEqualizedPaths{ith_cell}(:,1))>=length(cellArrayOfUnequalPaths{ith_cell}(:,1)))
end
assert(isequal(size(leastExtensionIndex),[1 1]));
assert(isequal(size(bestStartIndex),[1 1]));
assert(isequal(size(bestEndIndex),[1 1]));

% Check variable values
assert(isequal(round(cellArrayOfEqualizedPaths{1},4),round([0 0; cellArrayOfUnequalPaths{1}(1:end-1,:); 0.4545   11.000],4)));
assert(isequal(round(cellArrayOfEqualizedPaths{2},4),round([cellArrayOfUnequalPaths{2}],4)));
assert(isequal(round(cellArrayOfEqualizedPaths{3},4),round([4 0; cellArrayOfUnequalPaths{3}(1:end-1,:); 3.5455   11.0000],4)));
assert(isequal(leastExtensionIndex,2));
assert(isequal(bestStartIndex,2));
assert(isequal(bestEndIndex,2));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% BASIC CALL: Two of three ends hit each other
fig_num = 10009;
titleString = sprintf('BASIC CALL: Two of three ends hit each other');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

cellArrayOfUnequalPaths{1,1} = [0 5; 0 10; -0.5 10.9]; % Does not hit any
cellArrayOfUnequalPaths{2,1} = [2 0; 2 10; 2 11]; % Does not hit any
cellArrayOfUnequalPaths{3,1} = [4 1; 4 10; 4 10.9]; % Hits all

[cellArrayOfEqualizedPaths, leastExtensionIndex, bestStartIndex, bestEndIndex] = fcn_Path_equalizePathLengths(cellArrayOfUnequalPaths,(fig_num));

title(titleString, 'Interpreter','none');

% Check variable types
assert(iscell(cellArrayOfEqualizedPaths));
assert(isnumeric(leastExtensionIndex));
assert(isnumeric(bestStartIndex));
assert(isnumeric(bestEndIndex));

% Check variable sizes
% Show that the cell array has same number of paths after extensions
assert(isequal(size(cellArrayOfEqualizedPaths),[length(cellArrayOfUnequalPaths) 1]));
% Show extensions make the cell array longer
for ith_cell = 1:length(cellArrayOfUnequalPaths)
    assert(length(cellArrayOfEqualizedPaths{ith_cell}(:,1))>=length(cellArrayOfUnequalPaths{ith_cell}(:,1)))
end
assert(isequal(size(leastExtensionIndex),[1 1]));
assert(isequal(size(bestStartIndex),[1 1]));
assert(isequal(size(bestEndIndex),[1 1]));

% Check variable values
assert(isequal(round(cellArrayOfEqualizedPaths{1},4),round([0 0; cellArrayOfUnequalPaths{1}(1:end,:); -0.5556   11.0000],4)));
assert(isequal(round(cellArrayOfEqualizedPaths{2},4),round([cellArrayOfUnequalPaths{2}],4)));
assert(isequal(round(cellArrayOfEqualizedPaths{3},4),round([4 0; cellArrayOfUnequalPaths{3}(1:end-1,:); 4.000   11.0000],4)));
assert(isequal(leastExtensionIndex,2));
assert(isequal(bestStartIndex,2));
assert(isequal(bestEndIndex,2));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% ADVANCED CALL: using a realistic test path
fig_num = 20001;
titleString = sprintf('ADVANCED CALL: using a realistic test path');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

cellArrayOfUnequalPaths = fcn_INTERNAL_loadData;

[cellArrayOfEqualizedPaths, leastExtensionIndex, bestStartIndex, bestEndIndex] = fcn_Path_equalizePathLengths(cellArrayOfUnequalPaths,(fig_num));

title(titleString, 'Interpreter','none');

% Check variable types
assert(iscell(cellArrayOfEqualizedPaths));
assert(isnumeric(leastExtensionIndex));
assert(isnumeric(bestStartIndex));
assert(isnumeric(bestEndIndex));

% Check variable sizes
% Show that the cell array has same number of paths after extensions
assert(isequal(size(cellArrayOfEqualizedPaths),[length(cellArrayOfUnequalPaths) 1]));
% Show extensions make the cell array longer
for ith_cell = 1:length(cellArrayOfUnequalPaths)
    assert(length(cellArrayOfEqualizedPaths{ith_cell}(:,1))>=length(cellArrayOfUnequalPaths{ith_cell}(:,1)))
end
assert(isequal(size(leastExtensionIndex),[1 1]));
assert(isequal(size(bestStartIndex),[1 1]));
assert(isequal(size(bestEndIndex),[1 1]));

% Check variable values
% assert(isequal(cellArrayOfEqualizedPaths{1},[0 4; cellArrayOfUnequalPaths{1}]));
% assert(isequal(cellArrayOfEqualizedPaths{2},[cellArrayOfUnequalPaths{2}; 10 2]));
% assert(isequal(cellArrayOfEqualizedPaths{3},[0 0; cellArrayOfUnequalPaths{3}; 10 0]));
assert(isequal(leastExtensionIndex,2));
assert(isequal(bestStartIndex,1));
assert(isequal(bestEndIndex,2));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Fast Mode Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ______        _     __  __           _        _______        _
% |  ____|      | |   |  \/  |         | |      |__   __|      | |
% | |__ __ _ ___| |_  | \  / | ___   __| | ___     | | ___  ___| |_ ___
% |  __/ _` / __| __| | |\/| |/ _ \ / _` |/ _ \    | |/ _ \/ __| __/ __|
% | | | (_| \__ \ |_  | |  | | (_) | (_| |  __/    | |  __/\__ \ |_\__ \
% |_|  \__,_|___/\__| |_|  |_|\___/ \__,_|\___|    |_|\___||___/\__|___/
%
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Fast%20Mode%20Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures start with 8

close all;
fprintf(1,'Figure: 8XXXXXX: FAST mode cases\n');

%% Basic example - NO FIGURE
fig_num = 80001;
fprintf(1,'Figure: %.0f: Demo of fast mode, empty fig_num\n',fig_num);
figure(fig_num); close(fig_num);

cellArrayOfUnequalPaths{1,1} = [0 4; 7 4];
cellArrayOfUnequalPaths{2,1} = [3 2; 6 2];
cellArrayOfUnequalPaths{3,1} = [1 0; 10 0];

[cellArrayOfEqualizedPaths, leastExtensionIndex, bestStartIndex, bestEndIndex] = fcn_Path_equalizePathLengths(cellArrayOfUnequalPaths,([]));

% Check variable types
assert(iscell(cellArrayOfEqualizedPaths));
assert(isnumeric(leastExtensionIndex));
assert(isnumeric(bestStartIndex));
assert(isnumeric(bestEndIndex));

% Check variable sizes
% Show that the cell array has same number of paths after extensions
assert(isequal(size(cellArrayOfEqualizedPaths),[length(cellArrayOfUnequalPaths) 1]));
% Show extensions make the cell array longer
for ith_cell = 1:length(cellArrayOfUnequalPaths)
    assert(length(cellArrayOfEqualizedPaths{ith_cell}(:,1))>=length(cellArrayOfUnequalPaths{ith_cell}(:,1)))
end
assert(isequal(size(leastExtensionIndex),[1 1]));
assert(isequal(size(bestStartIndex),[1 1]));
assert(isequal(size(bestEndIndex),[1 1]));

% Check variable values
assert(isequal(cellArrayOfEqualizedPaths{1},[cellArrayOfUnequalPaths{1}; 10 4]));
assert(isequal(cellArrayOfEqualizedPaths{2},[0 2; cellArrayOfUnequalPaths{2}; 10 2]));
assert(isequal(cellArrayOfEqualizedPaths{3},[0 0; cellArrayOfUnequalPaths{3}]));
assert(isequal(leastExtensionIndex,3));
assert(isequal(bestStartIndex,1));
assert(isequal(bestEndIndex,3));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

cellArrayOfUnequalPaths{1,1} = [0 4; 7 4];
cellArrayOfUnequalPaths{2,1} = [3 2; 6 2];
cellArrayOfUnequalPaths{3,1} = [1 0; 10 0];

[cellArrayOfEqualizedPaths, leastExtensionIndex, bestStartIndex, bestEndIndex] = fcn_Path_equalizePathLengths(cellArrayOfUnequalPaths,(-1));

% Check variable types
assert(iscell(cellArrayOfEqualizedPaths));
assert(isnumeric(leastExtensionIndex));
assert(isnumeric(bestStartIndex));
assert(isnumeric(bestEndIndex));

% Check variable sizes
% Show that the cell array has same number of paths after extensions
assert(isequal(size(cellArrayOfEqualizedPaths),[length(cellArrayOfUnequalPaths) 1]));
% Show extensions make the cell array longer
for ith_cell = 1:length(cellArrayOfUnequalPaths)
    assert(length(cellArrayOfEqualizedPaths{ith_cell}(:,1))>=length(cellArrayOfUnequalPaths{ith_cell}(:,1)))
end
assert(isequal(size(leastExtensionIndex),[1 1]));
assert(isequal(size(bestStartIndex),[1 1]));
assert(isequal(size(bestEndIndex),[1 1]));

% Check variable values
assert(isequal(cellArrayOfEqualizedPaths{1},[cellArrayOfUnequalPaths{1}; 10 4]));
assert(isequal(cellArrayOfEqualizedPaths{2},[0 2; cellArrayOfUnequalPaths{2}; 10 2]));
assert(isequal(cellArrayOfEqualizedPaths{3},[0 0; cellArrayOfUnequalPaths{3}]));
assert(isequal(leastExtensionIndex,3));
assert(isequal(bestStartIndex,1));
assert(isequal(bestEndIndex,3));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

cellArrayOfUnequalPaths{1,1} = [0 4; 7 4];
cellArrayOfUnequalPaths{2,1} = [3 2; 6 2];
cellArrayOfUnequalPaths{3,1} = [1 0; 10 0];


Niterations = 50;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [cellArrayOfEqualizedPaths, leastExtensionIndex, bestStartIndex, bestEndIndex] = fcn_Path_equalizePathLengths(cellArrayOfUnequalPaths,([])); %#ok<ASGLU>
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [cellArrayOfEqualizedPaths, leastExtensionIndex, bestStartIndex, bestEndIndex] = fcn_Path_equalizePathLengths(cellArrayOfUnequalPaths,(-1)); %#ok<ASGLU>
end
fast_method = toc;

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

% Plot results as bar chart
figure(373737);
clf;
hold on;

X = categorical({'Normal mode','Fast mode'});
X = reordercats(X,{'Normal mode','Fast mode'}); % Forces bars to appear in this exact order, not alphabetized
Y = [slow_method fast_method ]*1000/Niterations;
bar(X,Y)
ylabel('Execution time (Milliseconds)')


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% BUG cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ____  _    _  _____
% |  _ \| |  | |/ ____|
% | |_) | |  | | |  __    ___ __ _ ___  ___  ___
% |  _ <| |  | | | |_ |  / __/ _` / __|/ _ \/ __|
% | |_) | |__| | |__| | | (_| (_| \__ \  __/\__ \
% |____/ \____/ \_____|  \___\__,_|___/\___||___/
%
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=BUG%20cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All bug case figures start with the number 9

close all;
fprintf(1,'Figure: 9XXXXXX: BUG mode cases\n');

%% Bug found when processing HSOV data
fig_num = 90001;
titleString = sprintf('BUG mode: Bug found when processing HSOV data');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

load('testData1_fcn_Path_equalizePathLengths','cellArrayOfPaths');

[cellArrayOfEqualizedPaths, leastExtensionIndex, bestStartIndex, bestEndIndex] = fcn_Path_equalizePathLengths(cellArrayOfPaths,(fig_num));

title(titleString, 'Interpreter','none');

% Check variable types
assert(iscell(cellArrayOfEqualizedPaths));
assert(isnumeric(leastExtensionIndex));
assert(isnumeric(bestStartIndex));
assert(isnumeric(bestEndIndex));

% Check variable sizes
% Show that the cell array has same number of paths after extensions
assert(isequal(size(cellArrayOfEqualizedPaths),[length(cellArrayOfUnequalPaths) 1]));
% Show extensions make the cell array longer
for ith_cell = 1:length(cellArrayOfUnequalPaths)
    assert(length(cellArrayOfEqualizedPaths{ith_cell}(:,1))>=length(cellArrayOfUnequalPaths{ith_cell}(:,1)))
end
assert(isequal(size(leastExtensionIndex),[1 1]));
assert(isequal(size(bestStartIndex),[1 1]));
assert(isequal(size(bestEndIndex),[1 1]));

% % Check variable values
% % assert(isequal(cellArrayOfEqualizedPaths{1},[0 4; cellArrayOfUnequalPaths{1}]));
% % assert(isequal(cellArrayOfEqualizedPaths{2},[cellArrayOfUnequalPaths{2}; 10 2]));
% % assert(isequal(cellArrayOfEqualizedPaths{3},[0 0; cellArrayOfUnequalPaths{3}; 10 0]));
% assert(isequal(leastExtensionIndex,2));
% assert(isequal(bestStartIndex,1));
% assert(isequal(bestEndIndex,2));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§

%% fcn_INTERNAL_loadData
function cellArrayOfUnequalPaths = fcn_INTERNAL_loadData

% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;
temp = paths_array';
cellArrayOfUnequalPaths = temp(1:3,:);

end % Ends fcn_INTERNAL_loadData