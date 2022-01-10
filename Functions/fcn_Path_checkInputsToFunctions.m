function fcn_Path_checkInputsToFunctions(...
    variable,variable_type_string)

% fcn_Path_checkInputsToFunctions
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
%      fcn_Path_checkInputsToFunctions(...
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
%            'path2or3D'  - checks that the path type is N x 2 or N x 3, with N>=2
%
%            'elevated_path'  - checks that the elevated path type is N x 3 
%            with N>=2
%
%            'paths'  - checks that the path type is N x 2 with N>=3
%
%            'traversal' - checks if a structure with X, Y, and Station,
%            and that each has an N x 1 vector within all of same length.
%            Further, the Station field must be strictly increasing.
%
%            'traversals' - checks if a structure containing a subfield
%            that is a cell array of traveral{i}, e.g. "data.traversal{3}",
%            with each traversal also meeting traversal requirements.
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
% This function was written on 2021_01_06 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
%      2021_01_06:
%      -- first write of the code moving this functionality out of
%      fcn_Path_FindOrthogonalScatterFromPathToPaths.m
%      2021_01_07:
%      -- fixed errors in comments as 'traversals' description was missing
%      -- created 'path' and 'paths' checks
%      2021_03_06:
%      -- created 'elevated_path' checks
%      2021_03_20:
%      -- created 'path2or3D' checks
%      -- updated traversal type to allow above as type, added comments


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
    if nargin < 2 || nargin > 2
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

%% Station
if strcmpi(variable_type_string,'station')
    % Check the station input
    if (length(variable(1,:))~=1) || ~isnumeric(variable)
        error('The %s input must be a station type, namely an N x 1 vector with N>=1',variable_name);
    end
end

%% Stations
if strcmpi(variable_type_string,'stations')
    % Check the station input type (e.g. stations must be a station)
    if (length(variable(1,:))~=1) || ~isnumeric(variable)
        error('The %s input must be a stations type, namely an N x 1 vector with N>1',variable_name);
    end
    
    % Check the station input row length
    if length(variable(:,1))==1
        error('The %s input must be a stations type, namely an N x 1 vector with N>1',variable_name);
    end
end

%% Path
% Must be an N x 2 with N>=2
if strcmpi(variable_type_string,'path')
    % Check the paths input type (e.g. paths must be a matrix of numbers)
    if (length(variable(1,:))~=2) || ~isnumeric(variable)
        error('The %s input must be a path type, namely an N x 2 vector with N>=2',variable_name);
    end
    
    % Check the paths input row length
    if length(variable(:,1))<2
        error('The %s input must be a path type, namely an N x 2 vector with N>=2',variable_name);
    end
end

%% path2or3D
% Must be an N x 2 or N x 3, with N>=2
if strcmpi(variable_type_string,'path2or3D')
    % Check the paths input type (e.g. paths must be a matrix of numbers)
    if (length(variable(1,:))~=2 && length(variable(1,:))~=3) || ~isnumeric(variable)
        error('The %s input must be a path2or3D type, namely an N x 2 or N x 3 vector with N>=2',variable_name);
    end
    
    % Check the paths input row length
    if length(variable(:,1))<2
        error('The %s input must be a path2or3D type, namely an N x 2 or N x 3 vector with N>=2',variable_name);
    end
end

%% Elevated Path
% Must be a N x 3 with N>=2
if strcmpi(variable_type_string,'elevated_path')
    % Check the elevated paths input type (e.g. elevated paths must be a matrix of numbers)
    if (length(variable(1,:))~=3) || ~isnumeric(variable)
        error('The %s input must be a elevated path type, namely an N x 3 vector with N>=2',variable_name);
    end
    
    % Check the paths input row length
    if length(variable(:,1))<2
        error('The %s input must be a elevated path type, namely an N x 3 vector with N>=2',variable_name);
    end
end

%% Paths
% Must be an N x 2 with N>=3
if strcmpi(variable_type_string,'paths')
    % Check the paths input type (e.g. paths must be a matrix of numbers)
    if (length(variable(1,:))~=2) || ~isnumeric(variable)
        error('The %s input must be a paths type, namely an N x 2 vector with N>=3',variable_name);
    end
    
    % Check the paths input row length
    if length(variable(:,1))<3
        error('The %s input must be a paths type, namely an N x 2 vector with N>=3',variable_name);
    end
end

%% Traversal
if strcmpi(variable_type_string,'traversal')
    % Check the reference_traversal subfields exist
    try
        X_central = variable.X;
        Y_central = variable.Y;
        Z_central = variable.Z;
        Station_central = variable.Station;
    catch
        error('The %s input must be a traversal type, namely being a structure with fields X, Y, Z, and Station, each N x 1 numeric arrays. The fields were not found. ',variable_name);
    end
    
    % Check that all are numeric
    if  ~isnumeric(X_central) ||  ~isnumeric(Y_central) ||  ~isnumeric(Z_central) ||  ~isnumeric(Station_central) 
        error('The %s input must be a traversal type, namely a structure with fields X, Y, Z, and Station, each N x 1 numeric arrays. At least one data field is non-numeric.',variable_name);
    end
    
    % Check that all are 1-dimensional columns
    if (length(X_central(1,:))~=1) || (length(Y_central(1,:))~=1) || (length(Z_central(1,:))~=1) ||  (length(Station_central(1,:))~=1)
        error('The %s input must be a traversal type, namely a structure with fields X, Y, Z, and Station, each N x 1 numeric arrays. At least one data field has multiple columns.',variable_name);
    end
    
    % Check that their lengths are all the same
    if (length(X_central(:,1))~=length(Y_central(:,1))) || ((length(X_central(:,1))~=length(Z_central(:,1))))  || ((length(X_central(:,1))~=length(Station_central(:,1))))
        warning('The %s input has variables whose lengths do not match. See the workspace for variable information.',variable_name);        
        fprintf(1,'\tX length is: %.0d.\n',length(X_central(:,1)));        
        fprintf(1,'\tY length is: %.0d.\n',length(Y_central(:,1)));      
        fprintf(1,'\tZ length is: %.0d.\n',length(Z_central(:,1)));
        fprintf(1,'\tStation length is: %.0d.\n',length(Station_central(:,1)));
        error('The %s input must be a traversal type, namely a structure with fields X, Y, Z, and Station, each N x 1 numeric arrays. The lengths do not match.',variable_name);
    end
       
    % Make sure the station field is sorted
    if ~issorted(Station_central,'strictascend')
        error('The Station field on the %s input must be strictly increasing',variable_name);
    end
end
    
%% Traversals
if strcmpi(variable_type_string,'traversals')
    % Check the all_traversals variables
    try
        test = variable.traversal{1}; %#ok<NASGU>
    catch
        error('The variable: %s is expected to be a traversals type, and must a subfield called traversal which is a cell array (e.g. variable.traversal{2}). An array was not found.',variable_name);
    end
    
    for i_traversal = 1:length(variable.traversal)
        try
        fcn_Path_checkInputsToFunctions(...
            variable.traversal{i_traversal},'traversal');
        catch ME
            error('The variable: %s is expected to be a traversals type, and must have a subfield called traversal which is a cell array; for eexample: variable.traversal{2}. An error in structure type was found in the %.0d index. The detail is: %s',variable_name,i_traversal,ME.message);
        end
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

