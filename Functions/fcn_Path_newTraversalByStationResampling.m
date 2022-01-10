function new_traversal = ...
    fcn_Path_newTraversalByStationResampling(input_traversal, new_stations, varargin)
% fcn_Path_newTraversalByStationResampling
% creates a new traversal by resampling a given traversal at given station
% points. 
%
% Note: if the stations are intended to align in space between the
% input_traversal and new_traversal traversals, then the first station
% point must be zero.
%
% If the stations are outside the station range of the input traversal,
% then extraploation is used to extend the input_traversal linearly
% outward. This can result in bad data if the path is not approximately
% linear at the endpoints.
%
% FORMAT: 
%
%      [new_traversal] = ...
%      fcn_Path_newTraversalByStationResampling(...
%            input_traversal,
%            new_stations,
%            (fig_num));
%
% INPUTS:
%
%      input_traversal: a traversal type data structure, namely a structure
%      with subfields of X, Y, Station, etc. in the following form
%           traversal.X
%
%      new_stations: an N x 1 column of stations that contain the locations
%      where the input traversal is to be estimated, with N>=1.
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      new_traversal: the resulting traversal representing the resampled
%      result
%
% DEPENDENCIES:
%
%      fcn_Path_checkInputsToFunctions
%      fcn_Path_plotTraversalsXY
%      fcn_Path_convertPathToTraversalStructure
%
% EXAMPLES:
%      
%     See the script: 
%     script_test_fcn_Path_newTraversalByStationResampling
%     for a full test suite.
%
% This function was written on 2022_01_05 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%     
%     2022_01_05: 
%     -- wrote the code originally - in support of
%     fcn_Path_findAverageTraversalViaOrthoProjection


flag_do_debug = 0; % Flag to show the results for debugging
flag_do_plots = 0; % Flag to plot the final results
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
    if nargin < 2 || nargin > 3
        error('Incorrect number of input arguments')
    end

    % Check the input_traversal input
    fcn_Path_checkInputsToFunctions(input_traversal, 'traversal');
    
    % Check the new_stations input
    fcn_Path_checkInputsToFunctions(new_stations, 'station');
          
end
    
% Does user want to show the plots?
if 3 == nargin
    fig_num = varargin{1};
    figure(fig_num);
    flag_do_plots = 1;
else
    if flag_do_debug
        fig = figure; 
        fig_num = fig.Number;
        flag_do_plots = 1;
    end
end

%% Main code starts here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

current_X = input_traversal.X;
current_Y = input_traversal.Y;
current_Station = input_traversal.Station;


interp_X       = interp1(current_Station,current_X,new_stations,'linear','extrap');
interp_Y       = interp1(current_Station,current_Y,new_stations,'linear','extrap');
new_traversal = fcn_Path_convertPathToTraversalStructure([interp_X, interp_Y]);

% The calculation of station is a bit off in the conversion, so fix it here
% in the cases where the user explicitly starts the station count at zero.
% For long station lists, this prevents round-off errors from accumulating.
if 0==new_stations(1,1)
    new_traversal.Station = new_stations; 
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
if flag_do_plots
    
    % plot the final XY result
    figure(fig_num);
    clf;
    
    hold on;
    plot(input_traversal.X,input_traversal.Y,'b.-','Linewidth',4,'Markersize',30);
    plot(new_traversal.X,new_traversal.Y,'r.-','Linewidth',2,'Markersize',20);
    legend('Input traversal','New traversal');
    title('original paths and final average path');
    xlabel('x [m]');
    ylabel('y [m]');
    
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end
end % End of function

