% function [bayes_avg, bayes_sigma] = fcn_bayesianAverage(value1, sigma1, value2, simgma2, varargin)
% varargin can take many arguments, arguments follow the order 
% value1, sigma1, value2, simag2 .......
% where 'value%d' is the parameter and 'sigma%d' is the corresponding std-deviation
% 'bayes_avg' is the weighted average
% 'bayes_sigma' is the weighted std-dev
% https://www.fil.ion.ucl.ac.uk/~wpenny/course/bayes.pdf
% https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
%
% Author: Srivenkata Satya Prasad Maddipatla (e.g. "Satya")
% Contact: szm888@psu.edu
% Date of first writing: 19th Oct, 2019
%
% Edited by Dr. Sean Brennan: sbrennan@psu.edu
%
% Revision history:
%   2019_10_19 - added more details to the example, see script script_check_fcn_bayesianAverage



function [bayes_avg, bayes_sigma] = fcn_bayesianAverage(value1, sigma1, value2, sigma2, varargin)
    
    m = length(varargin);       % Number of extra arguments to the function/more than 2-pairs
    value = [value1, value2, varargin{1:2:m}];    % Values satacked up as a matrix
    sigma = [sigma1, sigma2, varargin{2:2:m}];    % Corresponding std-dev stacked up
    
    variance = sigma.^2;                % Variance
    weights = 1 ./ variance;            % Inverse of variance
    weightValue = value .* weights;
    weighted_sum = sum(weightValue, 2); % Weighted sum
    sum_weights = sum(weights, 2);      % Sum of the inverse of variance
    
    bayes_avg = weighted_sum ./ sum_weights;
    bayes_sigma = 1 ./ sqrt(sum_weights);
    
end