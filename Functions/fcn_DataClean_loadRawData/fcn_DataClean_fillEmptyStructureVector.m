function EmptyVector = fcn_DataClean_fillEmptyStructureVector(structure)

% This function is used to create a zeros vector
% Input Variables:
%      structure = a data structure with ROS_Time field(format:struct)

% Returned Results:
%      EmptyVector
% Author: Liming Gao
% Created Date: 2020_11_15
% Modify Date: 2019_11_22
%
% Updates:
%
% To do lists:
% 1. 

%%
EmptyVector = 0*structure.ROS_Time;

return
