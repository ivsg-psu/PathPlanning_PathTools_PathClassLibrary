% script_test_fcn_DebugTools_convertVariableToCellString.m
% Tests fcn_DebugTools_convertVariableToCellString
% Written in 2023_01_18 by S.Brennan


% Test cases
%% Single character in cell array
result = fcn_DebugTools_convertVariableToCellString({'D'});
assert(isequal(result,{'D'}));

%% Empty single input
result = fcn_DebugTools_convertVariableToCellString({''});
assert(isempty(result{1}));

%% Numeric single input
result = fcn_DebugTools_convertVariableToCellString(2.34);
assert(isequal(result,{'2.34'}));

%% Multiple character in cell array
result = fcn_DebugTools_convertVariableToCellString([{'D'},{'abc'}]);
assert(isequal(result,{'D, abc'}));

%% Multiple mixed character, numeric in cell array
result = fcn_DebugTools_convertVariableToCellString([{'D'},{2}]);
assert(isequal(result,{'D, 2'}));

%% Multiple mixed character, numeric in cell array ending in space
result = fcn_DebugTools_convertVariableToCellString([{'D'},{2},'abc ']);
assert(isequal(result,{'D, 2, abc '}));

%% Multiple mixed character, numeric in cell array ending in string with commas
result = fcn_DebugTools_convertVariableToCellString([{'D'},{2},'abc , 123']);
assert(isequal(result,{'D, 2, abc , 123'}));

%% Character input
result = fcn_DebugTools_convertVariableToCellString('abc');
assert(isequal(result,{'abc'}));
