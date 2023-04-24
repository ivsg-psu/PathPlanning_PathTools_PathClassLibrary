
% script_test_fcn_DebugTools_number2string
% This is a script to exercise the function: fcn_DebugTools_number2string.m
% This function was written on 2021_12_12 by S. Brennan
% Questions or comments? sbrennan@psu.edu


% Revision history:
% 2022_11_14:
% -- wrote the code originally by copying out of old Exam2 code
% 2021_12_12:
% -- first write of the script
% 2023_02_17
% -- copied code out of AutoExam and into DebugTools

close all;
clc;


%% Basic case - example
stringNumber = fcn_DebugTools_number2string(2.333333333); % Empty result
assert(isequal(stringNumber,'2.33'));
fprintf(1,'%s\n',stringNumber);

%% Basic case - empty number
stringNumber = fcn_DebugTools_number2string([]); % Empty result
assert(isequal(stringNumber,' '));

%% Basic case - many forms of zero
stringNumber = fcn_DebugTools_number2string(0); % Zero
assert(isequal(stringNumber,'0'));

stringNumber = fcn_DebugTools_number2string(0.0000); % Zero
assert(isequal(stringNumber,'0'));

stringNumber = fcn_DebugTools_number2string(00.0000); % Zero
assert(isequal(stringNumber,'0'));

stringNumber = fcn_DebugTools_number2string(0.0); % Zero
assert(isequal(stringNumber,'0'));

stringNumber = fcn_DebugTools_number2string(0.); % Zero
assert(isequal(stringNumber,'0'));

%% Basic case - "round" numbers
stringNumber = fcn_DebugTools_number2string(12.000); 
assert(isequal(stringNumber,'12'));

stringNumber = fcn_DebugTools_number2string(985.000); 
assert(isequal(stringNumber,'985'));

stringNumber = fcn_DebugTools_number2string(2.0000);
assert(isequal(stringNumber,'2'));

stringNumber = fcn_DebugTools_number2string(555.0);
assert(isequal(stringNumber,'555'));

stringNumber = fcn_DebugTools_number2string(-32); 
assert(isequal(stringNumber,'-32'));

%% Basic case - "small" numbers
stringNumber = fcn_DebugTools_number2string(0.234); 
assert(isequal(stringNumber,'0.234000'));

%% Basic case - decimals greater than 1 or less than 10, show 2 decimal places
stringNumber = fcn_DebugTools_number2string(2.234); 
assert(isequal(stringNumber,'2.23'));

stringNumber = fcn_DebugTools_number2string(-2.237); 
assert(isequal(stringNumber,'-2.24'));

%% Basic case - decimals greater than 10 or less than 100, show 1 decimal places
stringNumber = fcn_DebugTools_number2string(12.234); 
assert(isequal(stringNumber,'12.2'));

stringNumber = fcn_DebugTools_number2string(-98.237); 
assert(isequal(stringNumber,'-98.2'));

%% Basic case - decimals greater than 100 or less than 100, show no decimal places
stringNumber = fcn_DebugTools_number2string(122.234); 
assert(isequal(stringNumber,'122'));

stringNumber = fcn_DebugTools_number2string(-988.237); 
assert(isequal(stringNumber,'-988'));

%% Fail conditions
if 1==0
    %% Bad input
    stringNumber = fcn_DebugTools_number2string('abc');
end
    



