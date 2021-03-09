%%
clear, clc

%% Load Data

load('MappingVan_DecisionMakingTestTrack_02242020.mat');
fn = fieldnames(MappingVan_DecisionMakingTestTrack_02242020);

%% 
% The start and end index are gotten from script_determineSeparationPoints
start_idx = [1,     15600, 42600, 57500, 67500, 91500];
end_idx   = [15360, 38550, 56500, 67500, 85950, 100000];
    
for i = 1:length(start_idx) % from 1 to 6 (6 different loops)
    
    for j = 1:length(fn) % from 1 to 12
        fn_str = string(fn(j));
        %GPS.(fn_str)(:, end_idx(i):end) = [];
        %GPS.(fn_str)(:, 1:start_idx(i)) = [];
        data_i.(fn_str) = MappingVan_DecisionMakingTestTrack_02242020.(fn_str);
        fn2 = fieldnames(data_i.(fn_str)); % second layer field names
        for k = 1:length
        
        
        
        
        %(:, start_idx(i):end_idx(i));
    end
    %file_name = ['loop_' num2str(i)];
    %save(file_name, "data_i")

end