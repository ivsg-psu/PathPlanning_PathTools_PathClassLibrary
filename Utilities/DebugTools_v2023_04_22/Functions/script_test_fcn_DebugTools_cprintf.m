% script_test_fcn_DebugTools_cprintf.m
% This is a script to exercise the function: fcn_DebugTools_cprintf.m
% This function was written on 2023_01_29 by S. Brennan
% Questions or comments? sbrennan@psu.edu


% Revision history:
% 2023_01_29:
% -- first write of the code (this is just a wrapper for cprint)

close all;
clc;



%% Built-in test cases
fcn_DebugTools_cprintf('text',    'regular black text ');
fcn_DebugTools_cprintf('hyper',   'followed %s','by ');
fcn_DebugTools_cprintf('key',     '%d colored ',5);
fcn_DebugTools_cprintf('-comment','& underlined ');
fcn_DebugTools_cprintf('err',     'elements:\n');
fcn_DebugTools_cprintf('cyan',    'cyan ');
fcn_DebugTools_cprintf('_green',  'underlined green ');
fcn_DebugTools_cprintf(-[1,0,1],  'underlined magenta ');
fcn_DebugTools_cprintf([1,0.5,0], 'and multi-\nline orange\n');
fcn_DebugTools_cprintf('*blue',   'and *bold* (R2011b+ only)\n');

%% Comprehensive list
fprintf(1,'\n');
fprintf(1,'Comprehensive list of fcn_DebugTools_cprintf options:\n');
fprintf(1,'\n');
fprintf(1,'Possible pre-defined STYLE names. NOTE: the STYLE entries are not case sensitive:\n');
fcn_DebugTools_cprintf('Text',                   '\t ''Text'' - default: black \n');
fcn_DebugTools_cprintf('Keywords',               '\t ''Keywords'' - default: blue \n');
fcn_DebugTools_cprintf('Comments',               '\t ''Comments'' - default: green \n');
fcn_DebugTools_cprintf('Strings',                '\t ''Strings'' - default: purple \n');
fcn_DebugTools_cprintf('UnterminatedStrings',    '\t ''UnterminatedStrings'' - default: dark red \n');
fcn_DebugTools_cprintf('SystemCommands',         '\t ''SystemCommands'' - default: orange \n');
fcn_DebugTools_cprintf('Errors',                 '\t ''Errors'' - default: light red \n');
fcn_DebugTools_cprintf('Hyperlinks',             '\t ''Hyperlinks'' - default: underlined blue \n');
fprintf(1,'\n');
fprintf(1,'Possible pre-defined COLOR names. NOTE: the COLOR entries are not case sensitive:\n')
fcn_DebugTools_cprintf('Black',                  '\t ''Black'' - default: black \n');
fcn_DebugTools_cprintf('Cyan',                   '\t ''Cyan'' - default: cyan \n');
fcn_DebugTools_cprintf('Magenta',                '\t ''Magenta'' - default: magenta \n');
fcn_DebugTools_cprintf('Blue',                   '\t ''Blue'' - default: blue \n');
fcn_DebugTools_cprintf('Green',                  '\t ''Green'' - default: green \n');
fcn_DebugTools_cprintf('Red',                    '\t ''Red'' - default: red \n');
fcn_DebugTools_cprintf('Yellow',                 '\t ''Yellow'' - default: yellow \n');
fcn_DebugTools_cprintf('White',                  '\t ''White'''); fcn_DebugTools_cprintf('Black',' - default: white \n');
fprintf(1,'\n');
fprintf(1,'Possible UNDERLINED (-) or (_) names. NOTE: not case sensitive:\n')
fcn_DebugTools_cprintf('-Text',                   '\t ''Text'' - default: black \n');
fcn_DebugTools_cprintf('-Keywords',               '\t ''Keywords'' - default: blue \n');
fcn_DebugTools_cprintf('-Comments',               '\t ''Comments'' - default: green \n');
fcn_DebugTools_cprintf('-Strings',                '\t ''Strings'' - default: purple \n');
fcn_DebugTools_cprintf('-UnterminatedStrings',    '\t ''UnterminatedStrings'' - default: dark red \n');
fcn_DebugTools_cprintf('-SystemCommands',         '\t ''SystemCommands'' - default: orange \n');
fcn_DebugTools_cprintf('-Errors',                 '\t ''Errors'' - default: light red \n');
fcn_DebugTools_cprintf('-Hyperlinks',             '\t ''Hyperlinks'' - default: underlined blue \n');
fcn_DebugTools_cprintf('-Black',                  '\t ''Black'' - default: black \n');
fcn_DebugTools_cprintf('-Cyan',                   '\t ''Cyan'' - default: cyan \n');
fcn_DebugTools_cprintf('-Magenta',                '\t ''Magenta'' - default: magenta \n');
fcn_DebugTools_cprintf('-Blue',                   '\t ''Blue'' - default: blue \n');
fcn_DebugTools_cprintf('-Green',                  '\t ''Green'' - default: green \n');
fcn_DebugTools_cprintf('-Red',                    '\t ''Red'' - default: red \n');
fcn_DebugTools_cprintf('-Yellow',                 '\t ''Yellow'' - default: yellow \n');
fcn_DebugTools_cprintf('-White',                  '\t ''White'''); fcn_DebugTools_cprintf('Black',' - default: white \n');
fprintf(1,'\n');
fprintf(1,'Possible BOLD (*) names. NOTE: not case sensitive:\n')
fcn_DebugTools_cprintf('*Text',                   '\t ''Text'' - default: black \n');
fcn_DebugTools_cprintf('*Keywords',               '\t ''Keywords'' - default: blue \n');
fcn_DebugTools_cprintf('*Comments',               '\t ''Comments'' - default: green \n');
fcn_DebugTools_cprintf('*Strings',                '\t ''Strings'' - default: purple \n');
fcn_DebugTools_cprintf('*UnterminatedStrings',    '\t ''UnterminatedStrings'' - default: dark red \n');
fcn_DebugTools_cprintf('*SystemCommands',         '\t ''SystemCommands'' - default: orange \n');
fcn_DebugTools_cprintf('*Errors',                 '\t ''Errors'' - default: light red \n');
fcn_DebugTools_cprintf('Hyperlinks',             '\t ''Hyperlinks'' - DOES NOT WORK!\n')
fcn_DebugTools_cprintf('*Black',                  '\t ''Black'' - default: black \n');
fcn_DebugTools_cprintf('*Cyan',                   '\t ''Cyan'' - default: cyan \n');
fcn_DebugTools_cprintf('*Magenta',                '\t ''Magenta'' - default: magenta \n');
fcn_DebugTools_cprintf('*Blue',                   '\t ''Blue'' - default: blue \n');
fcn_DebugTools_cprintf('*Green',                  '\t ''Green'' - default: green \n');
fcn_DebugTools_cprintf('*Red',                    '\t ''Red'' - default: red \n');
fcn_DebugTools_cprintf('*Yellow',                 '\t ''Yellow'' - default: yellow \n');
fcn_DebugTools_cprintf('*White',                  '\t ''White'''); fcn_DebugTools_cprintf('Black',' - default: white \n');
fprintf(1,'\n');
fprintf(1,'Possible BOLD (*) names. NOTE: not case sensitive:\n')
fcn_DebugTools_cprintf('*Text',                   '\t ''Text'' - default: black \n');
fcn_DebugTools_cprintf('*Keywords',               '\t ''Keywords'' - default: blue \n');
fcn_DebugTools_cprintf('*Comments',               '\t ''Comments'' - default: green \n');
fcn_DebugTools_cprintf('*Strings',                '\t ''Strings'' - default: purple \n');
fcn_DebugTools_cprintf('*UnterminatedStrings',    '\t ''UnterminatedStrings'' - default: dark red \n');
fcn_DebugTools_cprintf('*SystemCommands',         '\t ''SystemCommands'' - default: orange \n');
fcn_DebugTools_cprintf('*Errors',                 '\t ''Errors'' - default: light red \n');
fcn_DebugTools_cprintf('Hyperlinks',             '\t ''Hyperlinks'' - DOES NOT WORK!\n')
fcn_DebugTools_cprintf('*Black',                  '\t ''Black'' - default: black \n');
fcn_DebugTools_cprintf('*Cyan',                   '\t ''Cyan'' - default: cyan \n');
fcn_DebugTools_cprintf('*Magenta',                '\t ''Magenta'' - default: magenta \n');
fcn_DebugTools_cprintf('*Blue',                   '\t ''Blue'' - default: blue \n');
fcn_DebugTools_cprintf('*Green',                  '\t ''Green'' - default: green \n');
fcn_DebugTools_cprintf('*Red',                    '\t ''Red'' - default: red \n');
fcn_DebugTools_cprintf('*Yellow',                 '\t ''Yellow'' - default: yellow \n');
fcn_DebugTools_cprintf('*White',                  '\t ''White'''); fcn_DebugTools_cprintf('Black',' - default: white \n');
fprintf(1,'\n');
fprintf(1,'Color range listing examples: G are rows, B are columns\n')
for ith_R = 0:0.25:1
    fcn_DebugTools_cprintf([ith_R,0,0],'RGB setting: [%.1f G B]\n',ith_R);
    for ith_G = 0:0.25:1
        for ith_B = 0:0.25:1
            fcn_DebugTools_cprintf([ith_R,ith_G,ith_B],'[%.1f %.1f]',ith_G, ith_B);
        end
        fprintf(1,'\n');
    end
    fprintf(1,'\n');
end


%% Fail conditions
if 1==0
    %% Bad input
    output_string = fcn_DebugTools_addStringToEnd(input_string,value_to_add);
end
    