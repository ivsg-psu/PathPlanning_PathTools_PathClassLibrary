function fcn_GPS_checkInputsToFunctions(...
    variable,variable_type_string)

% fcn_GPS_checkInputsToFunctions
% Checks the variable types commonly used in the Path class to ensure that
% they are correctly formed. This function is typically called at the top
% of most Path class functions. The input is a variable and a string
% defining the "type" of the variable. This function checks to see that
% they are compatible. For example, the 'path' type variables in the class
% function are N x 2 arrays of [x y] pairs; if someone had a path variable
% called "path_example", they could check that this fit the path type by
% calling fcn_Path_checkInputsToFunctions(path_example,'path'). This
% function would then check that the array was N x 2, and if it was not, it
% would send out an error warning.
%
% FORMAT:
%
%      fcn_GPS_checkInputsToFunctions(...
%      variable,variable_type_string)
%
% INPUTS:
%
%      variable: the variable to check
%
%      variable_type_string: a string representing the variable type to
%      check. The current strings include:
%
%            'station' - checks that the station type is N x 1 and is a
%            number.
%
%            'stations' - checks that the station type is N x 1 and is a
%            number, with N >= 2
%
%            'path'  - checks that the path type is N x 2 with N>=2
%
%            'paths'  - checks that the path type is N x 2 with N>=3
%
%            'traversal' - checks if a structure with X, Y, and Station,
%            and that each has an N x 1 vector within all of same length.
%            Further, the Station field must be strictly increasing.
%
%
%      Note that the variable_type_string is not case sensitive: for
%      example, 'station' and 'Station' or 'STAtion' all give the same
%      result.
%
% OUTPUTS:
%
%      No explicit outputs, but produces MATLAB error outputs if conditions
%      not met, with explanation within the error outputs of the problem.
%
% EXAMPLES:
%
% See the script: script_test_fcn_Path_checkInputsToFunctions
% for a full test suite.
%
% DEPENDENCIES:
%
%      Uses MATLABs dbstack feature to trace dependencies 
%
% This function was written on 2021_01_25 by Satya Prasad
% Questions or comments? szm888@psu.edu

% Revision history:
%      2021_01_25:
%      -- first write of the code moving this functionality out of

flag_do_debug = 0; % Flag to debug the results
flag_do_plot = 0; % Flag to plot the results
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
end


%% check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag_check_inputs == 1   
    % Are there the right number of inputs?
    if 2 ~= nargin
        error('Incorrect number of input arguments')
    end
    
    % Check the string input, make sure it is characters
    if ~ischar(variable_type_string)
        error('The variable_type_string input must be a string type, for example: ''Path'' ');
    end
    
end

%% Start of main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grab the variable name
variable_name = inputname(1);

% Make the variable lower case
variable_type_string = lower(variable_type_string);

%% A path in Geodetic coordinate system
if strcmpi(variable_type_string,'path_LLA')
    % Check the size of inputs
    if 1 > size(variable,1) || 3 ~= size(variable,2)
        error('The %s input must be a path type (Geodetic), namely a N x 3 vector with N >= 1', variable_name);
    end
    
    % Check the type and validity of inputs
    if ~isnumeric(variable) || any(isnan(variable),'all')
        error('The %s input must be a path type (Geodetic), namely a N x 3 vector with N >= 1', variable_name);
    end
    
    % Check the domain of inputs (latitude and longitude)
    if (any(-90.0 > variable(:,1)) || any(90.0 < variable(:,1)) || ...
            any(-180.0 > variable(:,2)) || any(180.0 < variable(:,2)))
        error('The %s input is out of range', variable_name);
    end
end

%% A path in ECEF coordinate system
if strcmpi(variable_type_string,'path_ECEF')
    % Check the size of inputs
    if 1 > size(variable,1) || 3 ~= size(variable,2)
        error('The %s input must be a path type (ECEF), namely a N x 3 vector with N >= 1', variable_name);
    end
    
    % Check the type and validity of inputs
    if ~isnumeric(variable) || any(isnan(variable),'all')
        error('The %s input must be a path type (ECEF), namely a N x 3 vector with N >= 1', variable_name);
    end
end

%% A path in ENU coordinate system
if strcmpi(variable_type_string,'path_ENU')
    % Check the size of inputs
    if 1 > size(variable,1) || 3 ~= size(variable,2)
        error('The %s input must be a path type (ENU), namely a N x 3 vector with N >= 1', variable_name);
    end
    
    % Check the type and validity of inputs
    if ~isnumeric(variable) || any(isnan(variable),'all')
        error('The %s input must be a path type (ENU), namely a N x 3 vector with N >= 1', variable_name);
    end
end

%% A point in Geodetic coordinate system
if strcmpi(variable_type_string,'point_LLA')
    % Check the size of inputs
    if 1 ~= size(variable,1) || 3 ~= size(variable,2)
        error('The %s input must be a point type (Geodetic), namely an 1 x 3 vector', variable_name);
    end
    
    % Check the type and validity of inputs
    if ~isnumeric(variable) || any(isnan(variable))
        error('The %s input must be a point type (Geodetic), namely an 1 x 3 vector', variable_name);
    end
    
    % Check the domain of inputs (latitude and longitude)
    if ((-90.0 > variable(1,1)) || (90.0 < variable(1,1)) || ...
            (-180.0 > variable(1,2)) || (180.0 < variable(1,2)))
        error('The %s input is out of range', variable_name);
    end
end

%% A point in ECEF coordinate system
if strcmpi(variable_type_string,'point_ECEF')
    % Check the size of inputs
    if 1 ~= size(variable,1) || 3 ~= size(variable,2)
        error('The %s input must be a point type (ECEF), namely an 1 x 3 vector', variable_name);
    end
    
    % Check the type and validity of inputs
    if ~isnumeric(variable) || any(isnan(variable))
        error('The %s input must be a point type (ECEF), namely an 1 x 3 vector', variable_name);
    end
end

%% A point in ENU coordinate system
if strcmpi(variable_type_string,'point_ENU')
    % Check the size of inputs
    if 1 ~= size(variable,1) || 3 ~= size(variable,2)
        error('The %s input must be a point type (ENU), namely an 1 x 3 vector', variable_name);
    end
    
    % Check the type and validity of inputs
    if ~isnumeric(variable) || any(isnan(variable))
        error('The %s input must be a point type (ENU), namely an 1 x 3 vector', variable_name);
    end
end

%% Plot the results (for debugging)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _
%  |  __ \     | |
%  | |  | | ___| |__  _   _  __ _
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_do_plot
    % Nothing to plot here
end % Ends the flag_do_plot if statement

if flag_do_debug
    fprintf(1,'The variable: %s was checked that it meets type: %s, and no errors were detected.\n',variable_name,variable_type_string);
end
if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends the function
