% This is a wrapper script to run all the test scripts in the path
% class library for the purpose of evaluating every assertion test in these
% files
clear all; close all; clc;
all_scripts = dir('script_test_fcn_*.m');
diary 'script_test_fcn_path_all_stdout.txt';
for i_script = 1:length(all_scripts)
    file_name_extended = all_scripts(i_script).name;
    file_name = erase(file_name_extended,'.m');
    if ~strcmp(mfilename,file_name)
        file_name_trunc = erase(file_name,'script_');
        fprintf(1,'Testing script: %.0d of %.0d, %s\n\n',i_script,length(all_scripts),file_name_trunc);
        suite = testsuite(file_name);
        results = run(suite);
    end
end
