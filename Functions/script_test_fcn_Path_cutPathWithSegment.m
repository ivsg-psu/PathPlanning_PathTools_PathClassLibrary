% script_test_fcn_Path_cutPathWithSegment.m
% This is a script to exercise the function: fcn_Path_cutPathWithSegment.m
% This function was written on 2023_08_26 by S. Brennan, sbrennan@psu.edu
% FORMAT:
%
%    [cut_path_before, cut_path_after] = fcn_Path_cutPathWithSegment(pathToCut,cutting_segment,(figNum));

% Revision history:
% 2023_09_26 by S. Brennan
% - first write of the code

close all;

%% Basic Example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   ____            _        ______                           _      
%  |  _ \          (_)      |  ____|                         | |     
%  | |_) | __ _ ___ _  ___  | |__  __  ____ _ _ __ ___  _ __ | | ___ 
%  |  _ < / _` / __| |/ __| |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \
%  | |_) | (_| \__ \ | (__  | |____ >  < (_| | | | | | | |_) | |  __/
%  |____/ \__,_|___/_|\___| |______/_/\_\__,_|_| |_| |_| .__/|_|\___|
%                                                      | |           
%                                                      |_|          
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Basic%20Example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§



%% BASIC example
% A simple line segment, a simple query, zero distance in rear segments
figNum = 10001;
fprintf(1,'Figure: %.0f: Demo of fast mode, empty figNum\n',figNum);
figure(figNum); close(figNum);

pathToCut = [0 -1; 0 5];
cutting_segment = [-1 2; 1 2];

[cut_path_before, cut_path_after] = fcn_Path_cutPathWithSegment(pathToCut,cutting_segment,(figNum));

expected_cut_path_before = [0 -1;0 2];
expected_cut_path_after  = [0 2;0 5];
assert(max(sum((expected_cut_path_before - cut_path_before).^2,2).^0.5)<1E-10);
assert(max(sum((expected_cut_path_after - cut_path_after).^2,2).^0.5)<1E-10);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% BASIC example - cut line on end
% A simple line segment, a simple query, zero distance in rear segments
figNum = 10002;
fprintf(1,'Figure: %.0f: Demo of fast mode, empty figNum\n',figNum);
figure(figNum); close(figNum);

pathToCut = [0 -1; 0 5];
cutting_segment = [-1 5; 1 5];

[cut_path_before, cut_path_after] = fcn_Path_cutPathWithSegment(pathToCut,cutting_segment,(figNum));

expected_cut_path_before = [0 -1;0 5];
expected_cut_path_after  = [0 5];
assert(max(sum((expected_cut_path_before - cut_path_before).^2,2).^0.5)<1E-10);
assert(max(sum((expected_cut_path_after - cut_path_after).^2,2).^0.5)<1E-10);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% BASIC example - cut line on start
% A simple line segment, a simple query, zero distance in rear segments
figNum = 10003;
fprintf(1,'Figure: %.0f: Demo of fast mode, empty figNum\n',figNum);
figure(figNum); close(figNum);

pathToCut = [0 -1; 0 5];
cutting_segment = [-1 -1; 1 -1];

[cut_path_before, cut_path_after] = fcn_Path_cutPathWithSegment(pathToCut,cutting_segment,(figNum));

expected_cut_path_before = [0 -1];
expected_cut_path_after  = [0 -1; 0 5];
assert(max(sum((expected_cut_path_before - cut_path_before).^2,2).^0.5)<1E-10);
assert(max(sum((expected_cut_path_after - cut_path_after).^2,2).^0.5)<1E-10);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% BASIC example - cut line that misses
% A simple line segment, a simple query, zero distance in rear segments
figNum = 10004;
fprintf(1,'Figure: %.0f: Demo of fast mode, empty figNum\n',figNum);
figure(figNum); close(figNum);

pathToCut = [0 -1; 0 5];
cutting_segment = [-1 6; 1 6];

[cut_path_before, cut_path_after] = fcn_Path_cutPathWithSegment(pathToCut,cutting_segment,(figNum));

expected_cut_path_before = [];
expected_cut_path_after  = [];
assert(isempty(cut_path_before));
assert(isempty(cut_path_after));


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

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
fprintf(1,'Figure: 8XXXXXX: Demo of fast mode cases\n');

%% Basic example - NO FIGURE
figNum = 80001;
fprintf(1,'Figure: %.0f: Demo of fast mode, empty figNum\n',figNum);
figure(figNum); close(figNum);

pathToCut = [0 -1; 0 5];
cutting_segment = [-1 2; 1 2];

[cut_path_before, cut_path_after] = fcn_Path_cutPathWithSegment(pathToCut,cutting_segment,([]));

expected_cut_path_before = [0 -1;0 2];
expected_cut_path_after  = [0 2;0 5];
assert(max(sum((expected_cut_path_before - cut_path_before).^2,2).^0.5)<1E-10);
assert(max(sum((expected_cut_path_after - cut_path_after).^2,2).^0.5)<1E-10);

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Basic fast mode - NO FIGURE, FAST MODE
figNum = 80002;
fprintf(1,'Figure: %.0f: Demo of fast mode, figNum=-1\n',figNum);
figure(figNum); close(figNum);

pathToCut = [0 -1; 0 5];
cutting_segment = [-1 2; 1 2];

[cut_path_before, cut_path_after] = fcn_Path_cutPathWithSegment(pathToCut,cutting_segment,(-1));

expected_cut_path_before = [0 -1;0 2];
expected_cut_path_after  = [0 2;0 5];
assert(max(sum((expected_cut_path_before - cut_path_before).^2,2).^0.5)<1E-10);
assert(max(sum((expected_cut_path_after - cut_path_after).^2,2).^0.5)<1E-10);


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
figNum = 80003;
fprintf(1,'Figure: %.0f: Fast mode comparisons\n',figNum);
figure(figNum);
close(figNum);

pathToCut = [0 -1; 0 5];
cutting_segment = [-1 2; 1 2];

Niterations = 20;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [cut_path_before, cut_path_after] = fcn_Path_cutPathWithSegment(pathToCut,cutting_segment,([]));

end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [cut_path_before, cut_path_after] = fcn_Path_cutPathWithSegment(pathToCut,cutting_segment,(-1));

end
fast_method = toc;

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));

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
assert(~any(figHandles==figNum));


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

% close all;

%% BUG 



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
