% script_test_fcn_DebugTools_extractNumberFromStringCell
% Tests fcn_DebugTools_extractNumberFromStringCell
% Written in 2022_11_15 by S.Brennan


% Test cases
%% Single character in cell array
result = fcn_DebugTools_extractNumberFromStringCell({'1 turtle'});
assert(isequal(result,{'1'}));

%% Decimal number in cell array
result = fcn_DebugTools_extractNumberFromStringCell({'0.4'});
assert(isequal(result,{'0.4'}));

%% Empty character in cell array
result = fcn_DebugTools_extractNumberFromStringCell({''});
assert(isempty(result));

%% No number characters in cell array
result = fcn_DebugTools_extractNumberFromStringCell({'   '});
assert(isempty(result));

%% Negative value in cell array
result = fcn_DebugTools_extractNumberFromStringCell({'-1 turtles'});
assert(isequal(result,{'-1'}));

%% Two numbers in a cell array 
input = {'1 if by land and 2 if by sea'};
result = fcn_DebugTools_extractNumberFromStringCell(input);
assert(isequal(result,{'1'}));

%% No number in cell array - should be empty
result = fcn_DebugTools_extractNumberFromStringCell({'turtles'});
assert(isempty(result));

%% Many numbers in a cell array, and the one we seek repeats. Should grab first one only.
input = {'4 score and 7 years ago, 4 and 20 blackbirds'};
result = fcn_DebugTools_extractNumberFromStringCell(input);
assert(isequal(result,{'4'}));

%% PSU ID number embeded in text
input = {'abc123@psu.edu'};
result = fcn_DebugTools_extractNumberFromStringCell(input);
assert(isequal(result,{'123'}));

%% Number with a leading zero, need to remove this
input = {'016'};
result = fcn_DebugTools_extractNumberFromStringCell(input);
assert(isequal(result,{'16'}));

%% Number with lots of leading zeros, need to remove this
input = {'0000016'};
result = fcn_DebugTools_extractNumberFromStringCell(input);
assert(isequal(result,{'16'}));

%% Number with lots of leading zeros, but only zero, need to simplify
input = {'00000'};
result = fcn_DebugTools_extractNumberFromStringCell(input);
assert(isequal(result,{'0'}));

%% Decimal number in cell array with leading zeros
result = fcn_DebugTools_extractNumberFromStringCell({'0000.4'});
assert(isequal(result,{'0.4'}));

%% Decimal number in cell array with no zeros
result = fcn_DebugTools_extractNumberFromStringCell({'.4'});
assert(isequal(result,{'0.4'}));

%% Negative decimal number in cell array, with leading zero
result = fcn_DebugTools_extractNumberFromStringCell({'-0.4'});
assert(isequal(result,{'-0.4'}));

%% Negative decimal number in cell array, no leadng zero
result = fcn_DebugTools_extractNumberFromStringCell({'-.4'});
assert(isequal(result,{'-0.4'}));

%% Decimal number, negative, in cell array with leading zeros
result = fcn_DebugTools_extractNumberFromStringCell({'-0000.4'});
assert(isequal(result,{'-0.4'}));

%% Decimal number, negative, in cell array with leading zeros and text
result = fcn_DebugTools_extractNumberFromStringCell({'My number is -0000.4'});
assert(isequal(result,{'-0.4'}));

%% Show that cannot catch a leading negative if separated from the number!
result = fcn_DebugTools_extractNumberFromStringCell({'My number is - 0000.4'});
assert(isequal(result,{'0.4'}));

%% Number with lots of leading zeros, need to remove this
input = {'001600'};
result = fcn_DebugTools_extractNumberFromStringCell(input);
assert(isequal(result,{'1600'}));

%% Negative number with lots of leading zeros, need to remove this
input = {'-0000016'};
result = fcn_DebugTools_extractNumberFromStringCell(input);
assert(isequal(result,{'-16'}));

