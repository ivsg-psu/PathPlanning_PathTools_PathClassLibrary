% script_test_fcn_Path_calcSingleTraversalStandardDeviation.m
% Tests fcn_Path_calcSingleTraversalStandardDeviation
       
% Revision history:
% 2021_01_05
% -- first write of the code

close all
clc

% Clear any old variables
clear all_traversals

% Fill in sample paths (as a starter)
paths = fcn_Path_fillSamplePaths;

% Pick first path 1s reference_traversal structure
reference_traversal = fcn_Path_convertPathToTraversalStructure(paths{1});
all_traversals.traversal{1} = reference_traversal;


%% Test case 1: basic call for one trajectory
std_deviation = fcn_Path_calcSingleTraversalStandardDeviation(reference_traversal);


%% Test case 2: advanced call for one trajectory - specify figure
fig_num = 22;
std_deviation = [];
std_deviation = fcn_Path_calcSingleTraversalStandardDeviation(reference_traversal,fig_num);