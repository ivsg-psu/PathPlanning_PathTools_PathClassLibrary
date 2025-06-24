% script_test_fcn_Path_fillPathViaUserInputs
% This is a script to exercise the function:
% fcn_Path_fillPathViaUserInputs.m
% This function was written on 2020_10_15 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revisions
% 2020_10_15
% -- First write of the code
% 2021_01_11
% -- Modified to make the code reentrant for stand-alone operation
% 2021_01_20
% -- Fixed minor error in character being used as a number, and missing
% argument on fprintf function

close all;



%% BASIC example 1
if 1==0  % Intentionally comment this out so that doesn't autorun. Forces user to read the script!
    
    %% Code starts here
    fig_num = 1;
    h = figure(fig_num);
    hold on;
    
    num_iterations = input('How many paths do you want to draw? [Hit enter for default of 3]:','s');
    if isempty(num_iterations)
        num_iterations = 3;
    else
        num_iterations = str2double(num_iterations);
    end
    fprintf(1,'\n Filling in %.0d paths.\n',num_iterations);
    fprintf(1,'Instructions: \n');
    fprintf(1,'Left click on the plot to create points. \n');
    fprintf(1,'Right click on the plot to remove points \n');
    fprintf(1,'Double click on the plot to end the path creation. \n');
    fprintf(1,'When the last path is completed, another plot will be created to show results. \n');
    
    
    % Initialize the paths_array
    clear paths_array
    paths_array{num_iterations} = [0 0];
    for i_path = 1:num_iterations
        
        % Set the title header
        UserData.title_header = sprintf('Path %.0d of %.0d',i_path,num_iterations);
        
        % Save the results
        set(gcf,'UserData',UserData);
        
        pathXY = fcn_Path_fillPathViaUserInputs(fig_num);
        paths_array{i_path} = pathXY;
    end
    
    % Plot the results
    clear data;
    % Convert paths to traversal structures
    for i_Path = 1:length(paths_array)
        traversal = fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
        data.traversal{i_Path} = traversal;
    end
    
    % Plot the results
    fig_num = 12;
    fcn_Path_plotTraversalsYaw(data,fig_num);
    fig_num = 13;
    fcn_Path_plotTraversalsXY(data,fig_num);
    
end