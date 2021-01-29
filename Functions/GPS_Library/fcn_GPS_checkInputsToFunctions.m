function fcn_GPS_checkInputsToFunctions(...
    variable,variable_type_string)

% fcn_GPS_checkInputsToFunctions
% Checks the variable types commonly used in the GPS class to ensure that
% they are correctly formed. This function is typically called at the top
% of most GPS class functions. The input is a variable and a string
% defining the "type" of the variable. This function checks to see that
% they are compatible. For example, the 'path_LLA' type variables in the 
% class function are N x 3 arrays of [latitude longitude altitude] pairs; 
% if someone had a path_LLA variable called "path_LLA_example", they could 
% check that this fit the path type by calling 
% fcn_Path_checkInputsToFunctions(path_LLA_example,'path_LLA'). This
% function would then check that the array was N x 3, and if it was not, it
% would send out an error warning. In this case, it also checks for domain
% of latitude and longitude.
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
%            'path_LLA' - checks that the path_LLA type is N x 3 and is a
%            number. It also checks for the domain of latitude and
%            longitude.
%
%            'path_ECEF' - checks that the path_ECEF type is N x 3 and is a
%            number.
%
%            'path_ENU' - checks that the path_ENU type is N x 3 and is a
%            number.
%
%            'point_LLA' - checks that the point_LLA type is 1 x 3 and is a
%            number. It also checks for the domain of latitude and
%            longitude.
%
%            'point_ECEF' - checks that the point_ECEF type is 1 x 3 and is
%            a number.
%
%            'point_ENU' - checks that the point_ENU type is 1 x 3 and is a
%            number.
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
%
% DEPENDENCIES:
%
%      Uses MATLABs dbstack feature to trace dependencies 
%
% This function was written on 2021_01_25 by Satya Prasad
% Questions or comments? szm888@psu.edu

% Revision history:
%      2021_01_25:
%      -- first write of the code moving this functionality out of all
%      functions

flag_do_debug = 0; % Flag to debug the results
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
if flag_do_debug
    fprintf(1,'The variable: %s was checked that it meets type: %s, and no errors were detected.\n',variable_name,variable_type_string);
end
if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends the function
