% This is an example of testing and untilizing the class defined in
% class_MapGen_DebugUtils

% script_test_class_MapGen_DebugUtils
% Tests class: class_MapGen_DebugUtils

% REVISION HISTORY:
% 2022_01_13
% -- first written by S. Harnett
% 2022_04_26
% -- explanations and demonstrations added by S. Harnett

%% startup
close all; clear all; clc;
%% clear environment variables that the class uses
setenv('ENV_FLAG_CHECK_INPUTS','')
setenv('ENV_FLAG__DO_PLOT','')
setenv('ENV_FLAG_DO_DEBUG','')

%% Purpose of this class
% This class is meant to be used to manage flags for turning on and off blocks
% of code for plotting, debugging, and input sanitation.  These flags can be
% set globally in env vars or locally on a per instance basis.  The class will
% manage conflict resoliution in the getter when local and global values differ.
% The main design principles are (1) explicit definition takes precident over implicit and
% (2) default behavior should show more information rather than less to avoid accidentally masking
% This class also includes features for logging and debugging.

%% Constructor default behavior
% If the user sets nothing at construciton, we can assert the flags have their default values
% (i.e. true flags) because we do not want to hide information from the user as a default

% construct without arguments
my_debugger = class_MapGen_DebugUtils();
% assert all flags defaulted to true
assert(my_debugger.get_local_flag_check_inputs);
assert(my_debugger.get_local_flag_do_plot);
assert(my_debugger.get_local_flag_do_debug);

%% Global environment variable behavior

% If a user sets global flag values, the global flag values will not be
% propogated to the local flags so as not to overwrite explicit local flag values,
% unless the user wants to overwrite the local instances with global values

% set a global flag to false
my_debugger.set_global_flag_check_inputs('0');
% assert that the local flag is still true despite the global value being false
assert(my_debugger.get_local_flag_check_inputs);
my_debugger.set_global_flag_do_plot('0');
assert(my_debugger.get_local_flag_do_plot);
my_debugger.set_global_flag_do_debug('0');
assert(my_debugger.get_local_flag_do_debug);

% Now we can propogate the global value to the local value, overwriting it

% reset the local value (which is currently true) from the global value (which is currently false)
my_debugger.set_local_from_global_flag_check_inputs();
% assert that now the local flag is false, per the global flag
assert(~my_debugger.get_local_flag_check_inputs());
my_debugger.set_local_from_global_flag_do_plot();
assert(~my_debugger.get_local_flag_do_plot());
my_debugger.set_local_from_global_flag_do_debug();
assert(~my_debugger.get_local_flag_do_debug());

% delete test fixture so we can test more constructor features
delete(my_debugger);

% The constructor will inherit global flag values if they are set prior to construction and no
% explicit constructor arguments are given (going with the design principle of explicit
% specification taking precident over implicit)
% So now that the global flags are still set to 0, construct a new instance of the object
% and assert that defaults aren't used, and all flags are false.

% construct a new debugger
my_other_debugger = class_MapGen_DebugUtils();
% assert that local flag values are currently false, despite no argument being given at construction
assert(~my_other_debugger.get_local_flag_check_inputs);
assert(~my_other_debugger.get_local_flag_do_plot);
assert(~my_other_debugger.get_local_flag_do_debug);

% delete test fixture so we can test more constructor features
delete(my_other_debugger);

%% Constructor with arguments

% With the global environment variables set, if we construct a new instance with input arguments,
% we would expect the new input arguments to take precidence over the global default values

% with global set to 0, construct with arguments for true flags
my_last_debugger = class_MapGen_DebugUtils(1,1,1);

% assert that the local values are now true in this instance despite the global default being false
assert(my_last_debugger.get_local_flag_check_inputs);
assert(my_last_debugger.get_local_flag_do_plot);
assert(my_last_debugger.get_local_flag_do_debug);


%% Test setters
% From before, the flags in this instance are now all true.  But we can use setters to reset
% individual flags in this instance to false.

% Call setter to reset a true flag to false
my_last_debugger.set_local_flag_check_inputs(0);

% assert that now the local flag is in fact false
assert(~my_last_debugger.get_local_flag_check_inputs);
my_last_debugger.set_local_flag_do_plot(0);
assert(~my_last_debugger.get_local_flag_do_plot);
my_last_debugger.set_local_flag_do_debug(0);
assert(~my_last_debugger.get_local_flag_do_debug);

%% Test other debugging features
%% Logging
% This class also allows for logging std out.  The following method starts logging.
my_last_debugger.start_logging();
% And this method stops the logging
my_last_debugger.end_logging();
% This method will return the stack trace that we just captured with the logger methods
st = my_last_debugger.get_debug_stack();

%% Non-repeating pseudo-random number generation
% The following method returns a number based on the current epoch time.  This is
% particularly useful for generating new known figure numbers whilst ensuring figure numbers
% will not repeat.
fig_num = my_last_debugger.get_fig_number();
% assert that fig_num is greater than the epoch time at which this test was
% originally written
assert(fig_num >= 1642251620957);

% delete test fixture for garbage collection
delete(my_last_debugger);
