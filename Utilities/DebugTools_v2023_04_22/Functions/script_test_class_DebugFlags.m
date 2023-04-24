% script_test_class_DebugFlags
% Tests class: class_DebugFlags

% REVISION HISTORY:
% 2022_01_13
% -- first written by S. Harnett

%% startup
close all; clear all; clc;
%% clear environment variables
setenv('ENV_FLAG_CHECK_INPUTS','')
setenv('ENV_FLAG__DO_PLOT','')
setenv('ENV_FLAG_DO_DEBUG','')

%% Test heirarchy in contructor

% set nothing, assert defaults (true flags)
my_debugger = class_DebugFlags();
assert(my_debugger.get_local_flag_check_inputs);
assert(my_debugger.get_local_flag_do_plot);
assert(my_debugger.get_local_flag_do_debug);
% set global, don't propogate, assert flag is unchanged
my_debugger.set_global_flag_check_inputs('0');
assert(my_debugger.get_local_flag_check_inputs);
my_debugger.set_global_flag_do_plot('0');
assert(my_debugger.get_local_flag_do_plot);
my_debugger.set_global_flag_do_debug('0');
assert(my_debugger.get_local_flag_do_debug);
% propogate global flag to local level, assert flag has changed
my_debugger.set_local_from_global_flag_check_inputs();
assert(~my_debugger.get_local_flag_check_inputs());
my_debugger.set_local_from_global_flag_do_plot();
assert(~my_debugger.get_local_flag_do_plot());
my_debugger.set_local_from_global_flag_do_debug();
assert(~my_debugger.get_local_flag_do_debug());
% delete test fixture
delete(my_debugger);

% with global still set to 0, construct and assert that defaults aren't
% used
my_other_debugger = class_DebugFlags();
assert(~my_other_debugger.get_local_flag_check_inputs);
assert(~my_other_debugger.get_local_flag_do_plot);
assert(~my_other_debugger.get_local_flag_do_debug);
% delete test fixture
delete(my_other_debugger);

% with global set to 0, construct with overridden local variables and
% assert that these are used
my_last_debugger = class_DebugFlags(1,1,1);
assert(my_last_debugger.get_local_flag_check_inputs);
assert(my_last_debugger.get_local_flag_do_plot);
assert(my_last_debugger.get_local_flag_do_debug);


%% Test setters
% construct, set local to not the default, assert local is not default
my_last_debugger = class_DebugFlags();
my_last_debugger.set_local_flag_check_inputs(0);
assert(~my_last_debugger.get_local_flag_check_inputs);
my_last_debugger.set_local_flag_do_plot(0);
assert(~my_last_debugger.get_local_flag_do_plot);
my_last_debugger.set_local_flag_do_debug(0);
assert(~my_last_debugger.get_local_flag_do_debug);

%% Test loggers
% test log start
my_last_debugger.start_logging();
% test log end
my_last_debugger.end_logging();

%% Test getters
% test get_debug_stack
st = my_last_debugger.get_debug_stack();
% test get_fig_number
fig_num = my_last_debugger.get_fig_number();
% assert that fig_num is greater than the epoch time at which this test was
% originally written
assert(fig_num >= 1642251620957);

% delete test fixture
delete(my_last_debugger);
