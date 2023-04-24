% script_test_fcn_DebugTools_doStringsMatch
% Tests fcn_DebugTools_doStringsMatch
% Written in 2022_12_09 by S.Brennan

% Revision history:
% 2022_12_09:
% -- wrote the code originally
% 2023_01_17:
% -- moved out of the AutoExam codeset, into DebugTools


% Test cases
%% simple string comparisons, identical characters that are same
result = fcn_DebugTools_doStringsMatch('a','a');
assert(result);

%% simple string comparisons, identical characters that are same, but different cases
student_answer = 'A';
correct_answers = 'a';
result = fcn_DebugTools_doStringsMatch(student_answer,correct_answers);
assert(result);

student_answer = 'a';
correct_answers = 'A';
result = fcn_DebugTools_doStringsMatch(student_answer,correct_answers);
assert(result);

%% simple string comparisons, student answer is part of correct answer so true
student_answer = 'a';
correct_answers = 'abc';
result = fcn_DebugTools_doStringsMatch(student_answer,correct_answers);
assert(result);

%% simple string comparisons, student answer is part of correct answer so true, checking case
student_answer = 'A';
correct_answers = 'abc';
result = fcn_DebugTools_doStringsMatch(student_answer,correct_answers);
assert(result);

%% simple string comparisons, student answer is part of correct answer so true, checking to produce false result if student repeats (FALSE)
student_answer = 'aa';
correct_answers = 'abc';
result = fcn_DebugTools_doStringsMatch(student_answer,correct_answers);
assert(result==false);

%% simple string comparisons, student answer is part of correct answer so true, checking repeats (FALSE)
student_answer = 'aaa';
correct_answers = 'abc';
result = fcn_DebugTools_doStringsMatch(student_answer,correct_answers);
assert(result==false);

%% simple string comparisons, student answer is part of correct answer but too many guesses (FALSE)
student_answer = 'ab';
correct_answers = 'abc';
result = fcn_DebugTools_doStringsMatch(student_answer,correct_answers);
assert(result==false);

%% simple string comparisons, student answer and correct answer mixed up (FALSE)
student_answer = 'abc';
correct_answers = 'a';
result = fcn_DebugTools_doStringsMatch(student_answer,correct_answers);
assert(result==false);

%% simple string comparisons, student answer and correct answer are complex strings
student_answer = 'yes';
correct_answers = 'yes';
result = fcn_DebugTools_doStringsMatch(student_answer,correct_answers);
assert(result);

%% simple string comparisons, student answer and correct answer are complex strings, case insensitive
student_answer = 'YeS';
correct_answers = 'yes';
result = fcn_DebugTools_doStringsMatch(student_answer,correct_answers);
assert(result);

student_answer = 'YeS';
correct_answers = 'YES';
result = fcn_DebugTools_doStringsMatch(student_answer,correct_answers);
assert(result);

student_answer = 'Y';
correct_answers = 'YES';
result = fcn_DebugTools_doStringsMatch(student_answer,correct_answers);
assert(result);

student_answer = 'N';
correct_answers = 'No';
result = fcn_DebugTools_doStringsMatch(student_answer,correct_answers);
assert(result);

student_answer = 'NO';
correct_answers = 'no';
result = fcn_DebugTools_doStringsMatch(student_answer,correct_answers);
assert(result);

%% Checking N/A entries, that they don't trigger correct answers (if "A" is an answer)
student_answer = 'N/A';
correct_answers = 'abc';
result = fcn_DebugTools_doStringsMatch(student_answer,correct_answers);
assert(result==false);