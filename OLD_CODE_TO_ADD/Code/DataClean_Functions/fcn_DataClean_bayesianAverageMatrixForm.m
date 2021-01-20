function [bayes_avg, bayes_sigma] = fcn_DataClean_bayesianAverageMatrixForm(dataMatrix,sigmaMatrix)
% where dataMatrix is a matrix of the data to be averaged, arranged in
% columns where column 1 is data 1, column 2 is data 2, etc. where
% sigmaMatrix is a matrix of the standard deviations of the data, arranged
% in columns where column 1 is the stddev of data 1, column 2 is the stddev
% of data 2, etc.
% For outputs:
% 'bayes_avg' is the weighted average
% 'bayes_sigma' is the weighted std-dev
%
% Based on the algorithms for the weighted aritmetic mean:
% https://www.fil.ion.ucl.ac.uk/~wpenny/course/bayes.pdf
% https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
%
% Authors: Srivenkata Satya Prasad Maddipatla (e.g. "Satya")
% Contact: szm888@psu.edu, sbrennan@psu.edu
%
% Revision history:
%     2019_10_19 - First write based on Satya's version, added more details to the example, see script script_check_fcn_bayesianAverage
%     2019_10_20 - added matrix form support, plus error checking
%     2020_11_10 - changed function names in prep for DataClean class



%% Error checking
% Check to see if data is an integer type - if so, it will cause problems
% when doing floating-point calculations and give odd results.
if isinteger(dataMatrix)
    error(message('MATLAB:var:integerClass'));
end

% Check to see if the number of arguments is correct
if nargin ~= 2
    error('Invalid number of arguments: Expecting two arguments, first is data matrix, second is sigma matrix.');
end

% Check to see if the matrices are of same size
if ~isequal(size(dataMatrix),size(sigmaMatrix))
    error('Inconsistent dimensions on input arguments. Data matrix must have same size as sigma matrix.');
end

% Check to see if the matrices do not have nan
if max(max(isnan(dataMatrix))) || max(max((isnan(sigmaMatrix))))
    error('NaN data detected in arguments. Data and sigma matrices must be numeric values.');
end

% Check to see if the std-dev is non-positive
if (1 ~= min(min((sign(sigmaMatrix)))))
    error('Non-positive detected in sigma matrix. Std-dev must be positive');
end

% NOTE: need to add other argument checks if not using built-in variance
% command, var, below. Look at var.m code for examples

%% Start data processing

value = dataMatrix;    % Values satacked up as a matrix
sigma = sigmaMatrix;   % Corresponding std-dev stacked up

% Seems like the "var" command could do much of this... should look at
% this...
variance     = sigma.^2;            % Variance
weights      = 1 ./ variance;       % Inverse of variance
weightValue  = value .* weights;
weighted_sum = sum(weightValue, 2); % Weighted sum
sum_weights  = sum(weights, 2);     % Sum of the inverse of variance

bayes_avg    = weighted_sum ./ sum_weights;
bayes_sigma  = 1 ./ sqrt(sum_weights);

end