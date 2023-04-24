% script_test_fcn_DebugTools_parseStringIntoCells.m
% Tests fcn_DebugTools_parseStringIntoCells
% Written in 2022_11_15 by S.Brennan


% Test cases
%% Single character in cell array
result = fcn_DebugTools_parseStringIntoCells({'D'});
assert(isequal(result,{'D'}));

%% Empty character in cell array
result = fcn_DebugTools_parseStringIntoCells({''});
assert(isempty(result));

%% Complex input
inputString = 'This,isatest,of';
result = fcn_DebugTools_parseStringIntoCells(inputString);
assert(isequal(result,[{'This'},{'isatest'},{'of'}]));

%% Complex input with white space and commas that should be ignored
inputString = 'This,   ,  ,  isatest,';
result = fcn_DebugTools_parseStringIntoCells(inputString);
assert(isequal(result,[{'This'},{'isatest'}]));

%% Check the merging of strings
inputString = 'This,   ,  ,  isatest,';
[result,result2] = fcn_DebugTools_parseStringIntoCells(inputString);
assert(isequal(result,[{'This'},{'isatest'}]));
assert(isequal(result2,'thisisatest'));
