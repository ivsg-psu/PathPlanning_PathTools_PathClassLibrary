inputs.start_A   = [0 0; 0 5;  -5 10; 5 15];
inputs.start_B   = [0 1; -5 6; 0 11; -5 16];
inputs.vector_A   = [5 0; 10 0; 5 0; 20 0];
inputs.vector_B   = [10 0; 5 0; 20 0; 5 0];
inputs.A_values   = [0; 1; 0.5; 0.5];

Niterations = 10000;

% Do calculation with pre-calculation, FAST_MODE on
inputs.fig_num   = -1;
for ith_test = 1:Niterations
    % Call the function
    actual.percentage_B = ...
        fcn_Path_convertPerA2PerB(...
        inputs.start_A, inputs.start_B, ...
        inputs.vector_A, inputs.vector_B, inputs.A_values,...
        (inputs.fig_num));

end