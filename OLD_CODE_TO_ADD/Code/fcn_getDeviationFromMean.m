function deviation = fcn_getDeviationFromMean(mean, data)


mean_projection_vector = [diff(mean(:,1)), diff(mean(:,2)), diff(mean(:,3))]; % [x(1)-x(2), y(1)-y(2); x(2)-x(3), y(2)-y(3);...]

% repeat the first row because the diff function produces incorrect length (1 less)
mean_projection_vector = [mean_projection_vector(1,:); mean_projection_vector]; 

%mean_projection_vector = [mean_projection_vector zeros(length(mean_projection_vector(:,1)), 1)]; % Turn the vector into an xyz vector by adding a col of zeros

% Need to make all trajectories the same length as mean


for i_traversal = 1:length(data.traversal)
    
    % Grab the deviation vector:[x-x_mean, y-y_mean, 0]
    deviation_vector = [data.traversal{i_traversal}.X'-mean_Data.mean_xEast, data.traversal{i_traversal}.Y'-mean_Data.mean_yNorth, 0*mean_Data.mean_yNorth];
   
    % Grab the yaw angle for this run
    yaw_this_traversal = aligned_Data_ByStation.traversal{i_traversal}.yaw';

    
    % Calculate the yaw angle errors for this run
    yaw_error_this_traversal = aligned_Data_ByStation.traversal{i_traversal}.yaw'-mean_Data.mean_yaw;
    yaw_error_this_traversal_unwrapped = mod(yaw_error_this_traversal+180,360)-180;
    
    % Calculate the distance from the mean
    dist_from_mean = sum(deviation_vector.^2,2).^0.5;
    all_dist = [all_dist, dist_from_mean];
    % Do the cross product so we can check which direction the error is in
    cross_product_projection_to_deviation = cross(mean_projection_vector,deviation_vector);
    
    % Keep just the z-component - the 3rd column
    sign_of_deviation = 2.0*(cross_product_projection_to_deviation(:,3)>0)-1;
    
    % Grab the error
    error_relative_to_path = dist_from_mean.*sign_of_deviation;
    
    % get the closest point from the mean on each lap?
    
    
    all_error_lap = [all_error_lap, error_relative_to_path];
    all_error = [all_error; error_relative_to_path]; %#ok<AGROW>
    all_yaw_error = [all_yaw_error; yaw_error_this_traversal]; %#ok<AGROW>
    all_yaw = [all_yaw, yaw_this_traversal]; %#ok<AGROW>
    all_yaw_error_unwrapped = [all_yaw_error_unwrapped, yaw_error_this_traversal_unwrapped]; %#ok<AGROW>
    
end

end

