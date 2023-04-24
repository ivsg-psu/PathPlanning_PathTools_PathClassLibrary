% This is an example of class definition and the use of
% environment variables taken from the
% PathPlanning_MapTools_MapGenClassLibrary repo

classdef class_MapGen_DebugUtils < handle

    % debugging flags and their default values
    % flags default to showing the user optional information to avoid masking...
    % potential issues
    % these are private and should be accessed via 'get*' and 'set*' methods
    properties (Access = private)
        flag_check_inputs
        flag_do_plot
        flag_do_debug
        st
    end

    methods (Access = public)

        % constructor can be used to set all flags at once in this order:...
        % flag_check_inputs, flag_do_plot, flag_do_debug
        function obj = class_MapGen_DebugUtils(varargin)
            env_flag_check_inputs = getenv('ENV_FLAG_CHECK_INPUTS');
            env_flag_do_plot = getenv('ENV_FLAG_DO_PLOT');
            env_flag_do_debug = getenv('ENV_FLAG_DO_DEBUG');
            % if the constructor is given args, use these values
            if nargin == 3
                obj.set_local_flag_check_inputs(varargin{1});
                obj.set_local_flag_do_plot(varargin{2});
                obj.set_local_flag_do_debug(varargin{3});
            % otherwise use the environment variables
            elseif ~isempty(env_flag_check_inputs) && ...
                (env_flag_check_inputs == '1' || env_flag_check_inputs == '0') && ...
                ~isempty(env_flag_do_plot) && ...
                (env_flag_do_plot == '1' || env_flag_do_plot == '0') && ...
                ~isempty(env_flag_do_debug) && ...
                (env_flag_do_debug == '1' || env_flag_do_debug == '0')
                obj.flag_check_inputs = env_flag_check_inputs == '1';
                obj.flag_do_plot = env_flag_do_plot == '1';
                obj.flag_do_debug = env_flag_do_debug == '1';
            % otherwise default to flagging things on
            else
                obj.flag_check_inputs = 1;
                obj.flag_do_plot = 1;
                obj.flag_do_debug = 1;
                fprintf(['Using default value of 1 for flag_check_inputs, ',...
                'flag_do_plot, and flag_do_debug\n'])
            end
        end

        function set_local_flag_check_inputs(obj, flag_value)
            if flag_value == 1 || flag_value == 0 || flag_value == true || flag_value == false
                obj.flag_check_inputs = flag_value;
            else
                fprintf(['Desired value for flag_check_inputs was %s.\n',...
                'It should be 0 or 1 (or true or false)\n'],...
                flag_value)
            end
        end

        function set_local_flag_do_plot(obj, flag_value)
            if flag_value == 1 || flag_value == 0 || flag_value == true || flag_value == false
                obj.flag_do_plot = flag_value;
            else
                fprintf(['Desired value for flag_do_plot was %s.\n',...
                'It should be 0 or 1 (or true or false)\n'],...
                flag_value)
            end
        end

        function set_local_flag_do_debug(obj, flag_value)
            if flag_value == 1 || flag_value == 0 || flag_value == true || flag_value == false
                obj.flag_do_debug = flag_value;
            else
                fprintf(['Desired value for flag_do_debug was %s.\n',...
                'It should be 0 or 1 (or true or false)\n'],...
                flag_value)
            end
        end

        function set_global_flag_check_inputs(obj,flag_value)
            if flag_value == '1' || flag_value == 1 || flag_value == true
                setenv('ENV_FLAG_CHECK_INPUTS','1');
            elseif flag_value == '0' || flag_value == 0 || flag_value == false
                setenv('ENV_FLAG_CHECK_INPUTS','0');
            else
                fprintf(['Environment variable ENV_FLAG_CHECK_INPUTS was %s.\n',...
                'It should be a boolean; an integer of 1 or 0,', ...
                'or a one character string of 1 or 0\n'],...
                flag_value)
            end
        end

        function set_global_flag_do_plot(obj,flag_value)
            if flag_value == '1' || flag_value == 1 || flag_value == true
                setenv('ENV_FLAG_DO_PLOT','1');
            elseif flag_value == '0' || flag_value == 0 || flag_value == false
                setenv('ENV_FLAG_DO_PLOT','0');
            else
                fprintf(['Environment variable ENV_FLAG_DO_PLOT was %s.\n',...
                'It should be a boolean; an integer of 1 or 0,',...
                'or a one character string of 1 or 0\n'],...
                flag_value)
            end
        end

        function set_global_flag_do_debug(obj,flag_value)
            if flag_value == '1' || flag_value == 1 || flag_value == true
                setenv('ENV_FLAG_DO_DEBUG','1');
            elseif flag_value == '0' || flag_value == 0 || flag_value == false
                setenv('ENV_FLAG_DO_DEBUG','0');
            else
                fprintf(['Environment variable ENV_FLAG_DO_DEBUG was %s.\n',...
                'It should be a boolean; an integer of 1 or 0,',...
                'or a one character string of 1 or 0\n'],...
                flag_value)
            end
        end

        function set_local_from_global_flag_check_inputs(obj)
            env_flag_check_inputs = getenv('ENV_FLAG_CHECK_INPUTS');
            if env_flag_check_inputs == '1'
                obj.flag_check_inputs = 1;
            elseif env_flag_check_inputs == '0'
                obj.flag_check_inputs = 0;
            else
                fprintf(['Environment variable ENV_FLAG_CHECK_INPUTS was %s.\n',...
                'It should be a one character string of 1 or 0\n'],...
                env_flag_check_inputs)
            end
        end

        function set_local_from_global_flag_do_plot(obj)
            env_flag_do_plot = getenv('ENV_FLAG_DO_PLOT');
            if env_flag_do_plot == '1'
                obj.flag_do_plot = 1;
            elseif env_flag_do_plot == '0'
                obj.flag_do_plot = 0;
            else
                fprintf(['Environment variable ENV_FLAG_DO_PLOT was %s.\n',...
                'It should be a one character string of 1 or 0\n'],...
                env_flag_do_plot)
            end
        end

        function set_local_from_global_flag_do_debug(obj)
            env_flag_do_debug = getenv('ENV_FLAG_DO_DEBUG');
            if env_flag_do_debug == '1'
                obj.flag_do_debug = 1;
            elseif env_flag_do_debug == '0'
                obj.flag_do_debug = 0;
            else
                fprintf(['Environment variable ENV_FLAG_DO_DEBUG was %s.\n',...
                'It should be a one character string of 1 or 0\ncl'],...
                env_flag_do_debug)
            end
        end

        % get methods will use local variables
        function flag_value = get_local_flag_check_inputs(obj)
            flag_value = obj.flag_check_inputs;
        end

        function flag_value = get_local_flag_do_plot(obj)
            flag_value = obj.flag_do_plot;
        end

        function flag_value = get_local_flag_do_debug(obj)
            flag_value = obj.flag_do_debug;
        end

        function start_logging(obj)
            obj.st = dbstack; %#ok<*UNRCH>
            fprintf(1,'STARTING function: %s, in file: %s\n',obj.st(2).name,obj.st(2).file)
        end

        function end_logging(obj)
            fprintf(1,'ENDING function: %s, in file: %s\n\n',obj.st(2).name,obj.st(2).file)
        end

        function st = get_debug_stack(obj)
            st = obj.st;
        end

        function fig_num = get_fig_number(obj)
        % gives epoch time in thou of seconds providing a seemingly random...
        % figure number that won't be reproduced unless this method is called...
        % twice in a thousandth of a second
            fig_num = convertTo(datetime,'epochtime','TicksPerSecond',1000);
        end

    end
end
