% script_test_fcn_DebugTools_addStringToEnd.m
% This is a script to exercise the function: fcn_DebugTools_addStringToEnd.m
% This function was written on 2021_12_12 by S. Brennan
% Questions or comments? sbrennan@psu.edu


% Revision history:
% 2023_01_17:
% -- first write of the code

close all;
clc;



%% Basic case - numeric (adds a space)
input_string = 'test';
value_to_add = 2;
output_string = fcn_DebugTools_addStringToEnd(input_string,value_to_add);
assert(isequal(output_string,'test 2'));

%% Basic case - cell (adds a space)
input_string = 'test';
value_to_add = {'2'};
output_string = fcn_DebugTools_addStringToEnd(input_string,value_to_add);
assert(isequal(output_string,'test 2'));

%% Basic case - string (adds a space)
input_string = 'test';
value_to_add = '2';
output_string = fcn_DebugTools_addStringToEnd(input_string,value_to_add);
assert(isequal(output_string,'test 2'));

%% Fail conditions
if 1==0
    %% Bad input
    output_string = fcn_DebugTools_addStringToEnd(input_string,value_to_add);
end
    