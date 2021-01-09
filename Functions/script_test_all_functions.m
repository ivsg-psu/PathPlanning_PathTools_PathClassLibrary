% script_test_all_functions
% This is a script to exercise all functions in the path class
% This function was written on 2020_10_10 by S. Brennan, sbrennan@psu.edu

% Revision history:
%     2021_01_09
%     -- first write of the code

flag_file_log = 0;

all_scripts = dir('script_test_fcn_*');
for i_script = 1:length(all_scripts)
    script_name = all_scripts(i_script).name;
    
    clc
    close all
    
    if flag_file_log
        fid = fopen('test_log.txt','a');  %#ok<*UNRCH> % Open file in append mode
        fprintf(fid,'Testing script: %s\n',script_name);
        fclose(fid);
    end
    
    % Evaluate each script name, dropping the .m at the end
    eval(script_name(1:end-2));

end

close all;
clc;
fprintf(1,'All files tested - none threw errors.\n');