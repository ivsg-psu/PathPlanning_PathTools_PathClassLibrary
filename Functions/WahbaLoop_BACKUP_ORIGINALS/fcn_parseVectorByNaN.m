function [data,indices] = fcn_parseVectorByNaN(x)
%fcn_parseVectorByNaN Parses the data into segments that do not contain a
%NaN value.
% For example
%    temp = 1:20; temp(8) = nan;
%    x = temp'
%    [data,indices] = fcn_parseVectorByNaN(x);
%    data{1}
%    data{2}

flag_nan = isnan(x);
if any(flag_nan) % This only happens if there's NaN within the data array
    x_test = [nan; x; nan]; % Pad start and end of vector with NaN to enable a search
    flag_nan = isnan(x_test); % Tag all the points that have NaNs within
    changes = diff(flag_nan); % Mark the points as transitioning into or out of NaN

    start_points = find(changes==-1); % Mark points transitioning out of NaN, these are starting points
    for i=1:length(start_points) % Loop through starting points
        start_index = start_points(i); % Grab current start point
        changes(1:start_index)=0; % Clear entries in changes up to this start point
        end_index = find(changes==1,1)-1; % Find the next end point
        data{i} = x(start_index:end_index,1); %#ok<*AGROW> % Grab the data between this start/end point combination
        indices{i} = (start_index:end_index)';
    end
else
    data{1} = x;
    indices{1} = (1:length(x(:,1)))';
end
return

