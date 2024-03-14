% script_test_fcn_Path_cutPathWithSegment.m
% This is a script to exercise the function: fcn_Path_cutPathWithSegment.m
% This function was written on 2023_08_26 by S. Brennan, sbrennan@psu.edu
% FORMAT:
%
%    [cut_path_before, cut_path_after] = fcn_Path_cutPathWithSegment(pathToCut,cutting_segment,(fig_num));

% Revision history:
% 2023_09_26 by S. Brennan
% -- first write of the code



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



%% BASIC example 1
% A simple line segment, a simple query, zero distance in rear segments
pathToCut = [0 -1; 0 5];
cutting_segment = [-1 2; 1 2];

core_fig_num = 2323;
fig_num = core_fig_num + 1;

[cut_path_before, cut_path_after] = fcn_Path_cutPathWithSegment(pathToCut,cutting_segment,(fig_num));

expected_cut_path_before = [0 -1;0 2];
expected_cut_path_after  = [0 2;0 5];
assert(max(sum((expected_cut_path_before - cut_path_before).^2,2).^0.5)<1E-10);
assert(max(sum((expected_cut_path_after - cut_path_after).^2,2).^0.5)<1E-10);

%% BASIC example 2 - cut line on end
% A simple line segment, a simple query, zero distance in rear segments
pathToCut = [0 -1; 0 5];
cutting_segment = [-1 5; 1 5];

core_fig_num = 2323;
fig_num = core_fig_num + 2;

[cut_path_before, cut_path_after] = fcn_Path_cutPathWithSegment(pathToCut,cutting_segment,(fig_num));

expected_cut_path_before = [0 -1;0 5];
expected_cut_path_after  = [0 5];
assert(max(sum((expected_cut_path_before - cut_path_before).^2,2).^0.5)<1E-10);
assert(max(sum((expected_cut_path_after - cut_path_after).^2,2).^0.5)<1E-10);


%% BASIC example 3 - cut line on start
% A simple line segment, a simple query, zero distance in rear segments
pathToCut = [0 -1; 0 5];
cutting_segment = [-1 -1; 1 -1];

core_fig_num = 2323;
fig_num = core_fig_num + 3;

[cut_path_before, cut_path_after] = fcn_Path_cutPathWithSegment(pathToCut,cutting_segment,(fig_num));

expected_cut_path_before = [0 -1];
expected_cut_path_after  = [0 -1; 0 5];
assert(max(sum((expected_cut_path_before - cut_path_before).^2,2).^0.5)<1E-10);
assert(max(sum((expected_cut_path_after - cut_path_after).^2,2).^0.5)<1E-10);

%% BASIC example 4 - cut line that misses
% A simple line segment, a simple query, zero distance in rear segments
pathToCut = [0 -1; 0 5];
cutting_segment = [-1 6; 1 6];

core_fig_num = 2323;
fig_num = core_fig_num + 4;

[cut_path_before, cut_path_after] = fcn_Path_cutPathWithSegment(pathToCut,cutting_segment,(fig_num));

expected_cut_path_before = [];
expected_cut_path_after  = [];
assert(isempty(cut_path_before));
assert(isempty(cut_path_after));
