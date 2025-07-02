function clean_path = fcn_Path_cleanPathFromForwardBackwardJogs...
    (path_with_jogs,varargin)
% Finds and removes situations where the path is jumping forward and
% backward. This is detected by finding situations where the angle between
% segments is more than a threshold (currently pi/4), and then taking these
% segments, and the one before and after, and removing them. It then
% re-scans the path (up to 3 times) to again check for these situations.
%
% FORMAT: 
%
%     clean_path = fcn_Path_cleanPathFromForwardBackwardJogs...
%     (path_with_jogs, (fig_num));
%
% INPUTS:
%
%      path_with_jogs: a paths type consisting of (N x 2) array with N>=3
%
%      (OPTIONAL INPUTS)
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%      clean_path: the resulting path
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_Path_calcDiffAnglesBetweenPathSegments
%
% EXAMPLES:
%      
%     See the script: 
%     script_test_fcn_Path_cleanPathFromForwardBackwardJogs
%     for a full test suite.
%
% This function was written on 2021_01_09 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2021_01_09:
% -- wrote the code originally
% 2025_06_23 - S. Brennan
% -- Updated debugging and input checks

% TO-DO
% (none)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 2; % The largest Number of argument inputs to the function
flag_max_speed = 0;
if (nargin==MAX_NARGIN && isequal(varargin{end},-1))
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
        narginchk(1,MAX_NARGIN);

        % Check the Path variables
        fcn_DebugTools_checkInputsToFunctions(path_with_jogs, 'paths');

    end
end

% Does user want to show the plots?
flag_do_plots = 0; % Default is to NOT show plots
if (0==flag_max_speed) && (MAX_NARGIN == nargin) 
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        fig_num = temp;
        figure(fig_num);
        flag_do_plots = 1;
    end
end

if flag_do_debug
    fig_debug = 9585; 
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

Npath = length(path_with_jogs(:,1));
iteration_count = 1;
flag_average_is_good = 0;


% For debugging
if flag_do_debug
    figure(fig_debug);
    clf;

    plot(path_with_jogs(:,1),path_with_jogs(:,2),'k.-');
    hold on;
end


working_path_with_jogs = path_with_jogs;
points_removed = [];
while (0==flag_average_is_good)  && (iteration_count<=3)
    % Calculate angle changes between points
    diff_angles = fcn_Path_calcDiffAnglesBetweenPathSegments(working_path_with_jogs, -1);
    diff_angles_fullLength = [0; diff_angles];
    
    % Find outliers
    outliers = find(abs(diff_angles_fullLength)>pi/4);

    if flag_do_debug
        figure(fig_debug);
        plot(working_path_with_jogs(outliers,1),working_path_with_jogs(outliers,2),'ro','MarkerSize',10);
    end

    % For debugging
    if flag_do_debug
        figure(888);
        clf;

        x_indices = 1:length(diff_angles_fullLength);
        plot(x_indices,diff_angles_fullLength*180/pi,'k.-','MarkerSize',10);
        hold on;
        plot(x_indices(outliers), diff_angles_fullLength(outliers)*180/pi,'ro');
        xlabel('Indices');
        ylabel('Angle change (deg)');
    end

    % Are there any back/forth jogs?
    if ~isempty(outliers)

        % If there are, find pairs, e.g. outliers in sequence. These occur
        % where the differences in outliers is 1
        outlierIndexDifferences = diff(outliers);
        pairStarts = outliers(outlierIndexDifferences==1);
        realOutliers = nan(length(pairStarts),1);
        if ~isempty(realOutliers)
            for ith_pair = 1:length(pairStarts)
                startingIndex = pairStarts(ith_pair);
                if startingIndex==1 || startingIndex>=Npath-1
                    error('jogs found at start/end - this will likely cause errors.');
                end

                % Find the segment that the jogs occur "within". We will use
                % this segment to find which one is outlier is most away from the
                % segment, and just eliminate that one. To do this, we convert
                % the segment into a vector, rotate the vector by 90 degrees,
                % take the dot product of the ortho vector with the vector to
                % each point, and find the largest magnitude value.
                segmentStartIndex = startingIndex-1;
                segmentEndIndex   = startingIndex+2;
                segmentStart = path_with_jogs(segmentStartIndex,:);
                segmentEnd   = path_with_jogs(segmentEndIndex,:);
                segmentVector = segmentEnd-segmentStart;
                segmentOrthoVector = segmentVector*[0 1; -1 0];
                outlier1Vector = path_with_jogs(startingIndex,:) - segmentStart;
                outlier2Vector = path_with_jogs(startingIndex+1,:) - segmentStart;
                % Take dot products
                orthoDistance1 = abs(sum(outlier1Vector.*segmentOrthoVector,2));
                orthoDistance2 = abs(sum(outlier2Vector.*segmentOrthoVector,2));
                if orthoDistance1>orthoDistance2
                    realOutliers(ith_pair) = startingIndex;
                else
                    realOutliers(ith_pair) = startingIndex+1;
                end
            end % Ends loop through pairs


            if flag_do_debug
                figure(fig_debug);
                plot(working_path_with_jogs(realOutliers,1),working_path_with_jogs(realOutliers,2),'rx','MarkerSize',10);
            end
        else % No back/forth jogs left!
            flag_average_is_good = 1;
        end
        % % ID adjacent points
        % outliers = unique([outliers;outliers+1;outliers-1]);
        % outliers = min(outliers,length(working_path_with_jogs)-1);
        % outliers = max(1,outliers);
        
        % Create a set of indices we will save
        indices = (1:length(working_path_with_jogs(:,1)))';
        indices(realOutliers) = 0;
        
        % Save the clean path
        clean_path = working_path_with_jogs(indices~=0,:);        
        
        % Save the points that were removed
        points_removed = [points_removed; working_path_with_jogs(indices==0,:)]; %#ok<AGROW>

        
    else % No back/forth jogs left!
        flag_average_is_good = 1;
        clean_path = working_path_with_jogs;
    end
    
    % Show results for debugging?
    if flag_do_debug
        figure(fig_debug);
        plot(clean_path(:,1),clean_path(:,2),'b-');
    end
    
    % Increment the iteration count
    iteration_count = iteration_count + 1;
    
    % Reset the path average for the next round
    working_path_with_jogs = clean_path;
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
    
    % Prep the figure for plotting
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end      
    
    % Is this 2D or 3D?
    dimension_of_points = length(path_with_jogs(1,:));

    % Find size of plotting domain
    max_plotValues = max(path_with_jogs);
    min_plotValues = min(path_with_jogs);
    sizePlot = max(max_plotValues) - min(min_plotValues);
    nudge = sizePlot*0.006; %#ok<NASGU>


    % Find size of plotting domain
    if flag_rescale_axis
        percent_larger = 0.3;
        axis_range = max_plotValues - min_plotValues;
        if (0==axis_range(1,1))
            axis_range(1,1) = 2/percent_larger;
        end
        if (0==axis_range(1,2))
            axis_range(1,2) = 2/percent_larger;
        end
        if dimension_of_points==3 && (0==axis_range(1,3))
            axis_range(1,3) = 2/percent_larger;
        end

        
        % Force the axis to be equal
        min_valuesInPlot = min(min_plotValues);
        max_valuesInPlot = max(max_plotValues);

        % Stretch the axes
        stretched_min_vertexValues = min_valuesInPlot - percent_larger.*axis_range;
        stretched_max_vertexValues = max_valuesInPlot + percent_larger.*axis_range;
        axesTogether = [stretched_min_vertexValues; stretched_max_vertexValues];
        newAxis = reshape(axesTogether, 1, []);
        axis(newAxis);

    end
    % goodAxis = axis;

    hold on;
    grid on;
    grid minor;
    axis equal;
    
    plot(path_with_jogs(:,1),path_with_jogs(:,2),'k.-','Linewidth',5,'Markersize',40,'DisplayName','Path with jogs');
    plot(points_removed(:,1),points_removed(:,2),'ro','DisplayName','Outliers','MarkerSize',10,'LineWidth',3);
    plot(clean_path(:,1),clean_path(:,2),'c.-','Linewidth',2,'Markersize',25,'DisplayName','Fixed path, no jogs');
    legend;

    xlabel('X [m]')
    ylabel('Y [m]')

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end
end % End of function



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


