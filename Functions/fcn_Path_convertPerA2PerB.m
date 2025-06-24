function B_values = ...
    fcn_Path_convertPerA2PerB(...
    start_A, start_B, ...
    vector_A, vector_B, A_values,...
    varargin)   
%% fcn_Path_convertPerA2PerB 
% give the starting points and vectors for two vectors A and B, and the
% length(s) along vector A, calculates lengths along vector B
%
% FORMAT: 
%
%      percentage_B = ...
%         fcn_Path_convertPerA2PerB(...
%         start_A, start_B, ...
%         vector_A, vector_B, A_values,...
%         (fig_num))
%
% INPUTS:
%
%      start_A: an K x 2 vector containing the X,Y point of the start of
%      the A vector. K can either be 1, or must be equal to the J in vector B.
%      If K equals J, then the solution is found for each combination of
%      A, B, A_values. If K equals 1, the solution is found by fixing the A
%      vector and associating each A_value to each B value.
%
%      start_B: an J x 2 vector containing the X,Y point of the start of
%      the B vector. If J is not equal to K, then K must be 1. J must be
%      equal to L, the number of inputs in A_values
%
%      vector_A: an K x 2 vector containing the A vector
%
%      vector_B: an J x 2 vector containing the B vector
%
%      A_values: an L x 1 vector containing the scalars to use to find the
%      B values
%
%     (OPTIONAL INPUTS)
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%      percentage_B: a N x 1 vector representing the A distances along the
%      B vector
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%      
%       See the script: script_test_fcn_Path_convertPerA2PerB.m
%       for a full test suite. 
%
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%      2025_06_19 - S. Brennan
%      - wrote the code


%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==6 && isequal(varargin{end},-1))
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_PATHCLASS_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_PATHCLASS_FLAG_CHECK_INPUTS");
    MATLABFLAG_PATHCLASS_FLAG_DO_DEBUG = getenv("MATLABFLAG_PATHCLASS_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_PATHCLASS_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_PATHCLASS_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_PATHCLASS_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_PATHCLASS_FLAG_CHECK_INPUTS);
    end
end

% flag_do_debug = 1;

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 999978; %#ok<NASGU>
else
    debug_fig_num = []; %#ok<NASGU>
end

%% check input arguments?
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

if 0==flag_max_speed
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(5,6);

        % Check the start_A input
        fcn_DebugTools_checkInputsToFunctions(start_A, '2column_of_numbers');

        % Check the start_B input
        fcn_DebugTools_checkInputsToFunctions(start_B, '2column_of_numbers');

        % Check the vector_A input
        NinA = length(start_A(:,1));
        NinB = length(start_B(:,1));
        fcn_DebugTools_checkInputsToFunctions(vector_A, '2column_of_numbers', [2 NinA]);

        % Check the vector_B input
        fcn_DebugTools_checkInputsToFunctions(vector_B, '2column_of_numbers', [1 NinB]);

        % Check the A_values input
        fcn_DebugTools_checkInputsToFunctions(A_values, '1column_of_numbers');
    end
end

NinA = length(start_A(:,1));
NinB = length(start_B(:,1));
Ninputs = length(A_values(:,1));

% Figure out which situation we have
if NinA==NinB
    assert(Ninputs==Ninputs);    
    if Ninputs~=NinA
        start_A = ones(Ninputs,1)*start_A;
        vector_A = ones(Ninputs,1)*vector_A;
        start_B = ones(Ninputs,1)*start_B;
        vector_B = ones(Ninputs,1)*vector_B;
    end
else
    if NinA ==1
        % A and B have different sizes. Repeat the A vector values for
        % every value needed in the calculation.
        assert(NinB==Ninputs);
        start_A = ones(Ninputs,1)*start_A;
        vector_A = ones(Ninputs,1)*vector_A;
    else
        warning('on','backtrace');
        warning('Expecting a A and B vectors to have same length, or that A vector has length 1. Length of A and B are (respectively): %.0f %.0f',NinA, NinB);
        error('Cross product must be zero!');
    end
end


% Does user want to show the plots?
flag_do_plot = 0; % Default is to NOT show plots
if (0==flag_max_speed) && (6 == nargin) 
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp; 
        flag_do_plot = 1;
    end
else
    if flag_do_debug
        fig = figure; 
        fig_num = fig.Number; 
        flag_do_plot = 1;
    end
end


%% Calculations begin here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check the cross product?
if (0==flag_max_speed)
    A_cross_B = fcn_INTERNAL_crossProduct(vector_A,vector_B);
    if any(A_cross_B>(eps*1000))
        warning('on','backtrace');
        warning('Expecting a A and B vectors to be collinear, but their cross procudt is: %.5f',A_cross_B);
        error('Cross product must be zero!');
    end
end

% Calculate the result. See the PPT presentation for derivation of this
% formula for the B_values.
A_dot_B = sum(vector_A.*vector_B,2);
B_dot_B = sum(vector_B.*vector_B, 2);
As_minus_Bs_dot_B = sum((start_A - start_B).*vector_B,2);

B_values = A_values.*(A_dot_B)./B_dot_B + As_minus_Bs_dot_B./B_dot_B;

%% Any debugging?
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

    % For debugging
    predicted_A_points = A_values.*vector_A + start_A;
    predicted_B_points = B_values.*vector_B + start_B;

    % check whether the figure already has data
    h_fig = figure(fig_num);
    flag_rescale_axis = 0; 
    if isempty(get(h_fig,'Children')) 
        flag_rescale_axis = 1; 
    else
        child_handle = get(h_fig,'Children');
        if isfield(child_handle,'TileArrangement') && strcmp(get(child_handle,'TileArrangement'),'flow')
            flag_rescale_axis = 1;
        end
    end

    hold on;
    axis equal;
    grid on; 

    % Find size of plotting domain
    allPoints = [start_A; start_B; start_A+vector_A; start_B+vector_B; predicted_A_points; predicted_B_points];

    max_plotValues = max(allPoints);
    min_plotValues = min(allPoints);
    sizePlot = max(max_plotValues) - min(min_plotValues);
    nudge = sizePlot*0.006; %#ok<NASGU>

    % Set size of plotting domain
    if flag_rescale_axis

        percent_larger = 0.3;
        axis_range = max_plotValues - min_plotValues;
        if (0==axis_range(1,1))
            axis_range(1,1) = 2/percent_larger;
        end
        if (0==axis_range(1,2))
            axis_range(1,2) = 2/percent_larger;
        end

        % Force the axis to be equal
        min_vertexValuesInPlot = min(min_plotValues);
        max_vertexValuesInPlot = max(max_plotValues);

        % Stretch the axes
        stretched_min_vertexValues = min_vertexValuesInPlot - percent_larger.*axis_range;
        stretched_max_vertexValues = max_vertexValuesInPlot + percent_larger.*axis_range;
        axesTogether = [stretched_min_vertexValues; stretched_max_vertexValues];
        newAxis = reshape(axesTogether, 1, []);
        axis(newAxis);

    end
    % goodAxis = axis;

    % Plot the A and B vector

    % Plot the A points
    quiver(start_A(:,1),start_A(:,2),vector_A(:,1),vector_A(:,2),0, 'r','Linewidth',2,'MaxHeadSize',1);
    colorsUsed = zeros(Ninputs,3);
    for ith_input = 1:Ninputs
        h_plot = plot(predicted_A_points(ith_input,1),predicted_A_points(ith_input,2),'.','Linewidth',5,'Markersize',30);
        colorsUsed(ith_input,:) = get(h_plot,'Color');
        if NinA==1 && NinB==1
            plot(predicted_A_points(ith_input,1),predicted_A_points(ith_input,2),'k.','Linewidth',5,'Markersize',30);
        end
    end

    if NinB==1
        % Use the same color for all points, if there's only one B vector
        colorsUsed = ones(Ninputs,1)*colorsUsed(1,:);
    end

    % Plot the B inputs
    for ith_input = 1:Ninputs
        quiver(start_B(ith_input,1),start_B(ith_input,2),vector_B(ith_input,1),vector_B(ith_input,2),0, 'Linewidth',2,'MaxHeadSize',1,'Color',colorsUsed(ith_input,:));
        plot(predicted_B_points(ith_input,1),predicted_B_points(ith_input,2),'.','Linewidth',5,'Markersize',15,'Color',colorsUsed(ith_input,:));
    end
    
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end
end

%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§

%% fcn_INTERNAL_crossProduct
% Calculate cross products
function result = fcn_INTERNAL_crossProduct(v,w)
result = v(:,1).*w(:,2)-v(:,2).*w(:,1);
end % Ends fcn_INTERNAL_crossProduct
