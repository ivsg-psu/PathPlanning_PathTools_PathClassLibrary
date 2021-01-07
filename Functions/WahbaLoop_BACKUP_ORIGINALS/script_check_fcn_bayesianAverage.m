%% This is a script that checks the function: fcn_bayesianAverage
% function [bayes_avg, bayes_sigma] = fcn_bayesianAverage(value1, sigma1, value2, simgma2, varargin)
% varargin can take many arguments, arguments follow the order 
% value1, sigma1, value2, simag2 .......
% where 'value%d' is the parameter and 'sigma%d' is the corresponding std-deviation
% 'bayes_avg' is the weighted average
% 'bayes_sigma' is the weighted std-dev

% Prep the workspace
clear all
close all
clc

% create random values
A = 0.5*randn(50,1);
A_s = 0.5*ones(50,1);
B = [0.02*randn(15,1); 100*randn(20,1); 0.02*randn(15,1)];
B_s = [0.02*ones(15,1); 100*ones(20,1); 0.02*ones(15,1)];
C = [1*randn(15,1); 50*randn(20,1); 1*randn(15,1)];
C_s = [1*ones(15,1); 50*ones(20,1); 1*ones(15,1)];

bayes_avg = fcn_bayesianAverage(A,A_s,B,B_s,C,C_s);

figure()
plot(A(:,1),'b');   % Similar to novatel
hold on;
plot(B(:,1),'g');   % Similar to Hemisphere
plot(C(:,1),'C');   % Similar to Hemisphere
plot(bayes_avg,'r--');  % weighted
legend('0.5','0.02 and 100','1 and 50','avg');

%% Create better example
% First, create test data range, and bad indices
xvector = (1:20)';
bad_parts_of_x = (5:15);
clean_data1 = 4*ones(length(xvector),1);
clean_data2 = 6*ones(length(xvector),1);
clean_data3 = 8*ones(length(xvector),1);

% Next, add noise to each with a given sigma -
sigma1 = 0.2*ones(length(xvector),1);
sigma2 = 0.4*ones(length(xvector),1);
sigma3 = 1*ones(length(xvector),1);

% Show that the weighted average shifts toward variable of least variance
[bayes_avg,sigmabayes] = fcn_bayesianAverage(...
    clean_data1, sigma1,...
    clean_data2, sigma2,...
    clean_data3, sigma3);

% Plot results
figure(111);
clf;
plot(...
    xvector,clean_data1,'r-',...
    xvector,clean_data1+sigma1,'r--',...
    xvector,clean_data1-sigma1,'r--',...
    xvector,clean_data2,'g-',...
    xvector,clean_data2+sigma2,'g--',...
    xvector,clean_data2-sigma2,'g--',...
    xvector,clean_data3,'b-',...
    xvector,clean_data3+sigma3,'b--',...
    xvector,clean_data3-sigma3,'b--');
hold on;
plot(xvector,bayes_avg,'k-','LineWidth',2);
plot(xvector,bayes_avg+sigmabayes,'k--','LineWidth',1);
plot(xvector,bayes_avg-sigmabayes,'k--','LineWidth',1);

%% Now change the variances around to show that this shifts average

% Next, add noise to each with a given sigma -
sigma1 = 1*ones(length(xvector),1);
sigma2 = 0.4*ones(length(xvector),1);
sigma3 = 0.2*ones(length(xvector),1);

% Show that the weighted average shifts toward variable of least variance
[bayes_avg,sigmabayes] = fcn_bayesianAverage(...
    clean_data1, sigma1,...
    clean_data2, sigma2,...
    clean_data3, sigma3);

% Plot results
figure(222);
clf;
plot(...
    xvector,clean_data1,'r-',...
    xvector,clean_data1+sigma1,'r--',...
    xvector,clean_data1-sigma1,'r--',...
    xvector,clean_data2,'g-',...
    xvector,clean_data2+sigma2,'g--',...
    xvector,clean_data2-sigma2,'g--',...
    xvector,clean_data3,'b-',...
    xvector,clean_data3+sigma3,'b--',...
    xvector,clean_data3-sigma3,'b--');
hold on;
plot(xvector,bayes_avg,'k-','LineWidth',2);
plot(xvector,bayes_avg+sigmabayes,'k--','LineWidth',1);
plot(xvector,bayes_avg-sigmabayes,'k--','LineWidth',1);


%% Now change the variances around locally

% Next, add noise to each with a given sigma -
sigma1 = 1*ones(length(xvector),1);
sigma2 = 0.4*ones(length(xvector),1);
sigma3 = 0.2*ones(length(xvector),1);
sigma3(bad_parts_of_x,1) = 0.8;

% Show that the weighted average shifts toward variable of least variance
[bayes_avg,sigmabayes] = fcn_bayesianAverage(...
    clean_data1, sigma1,...
    clean_data2, sigma2,...
    clean_data3, sigma3);

% Plot results
figure(333);
clf;
plot(...
    xvector,clean_data1,'r-',...
    xvector,clean_data1+sigma1,'r--',...
    xvector,clean_data1-sigma1,'r--',...
    xvector,clean_data2,'g-',...
    xvector,clean_data2+sigma2,'g--',...
    xvector,clean_data2-sigma2,'g--',...
    xvector,clean_data3,'b-',...
    xvector,clean_data3+sigma3,'b--',...
    xvector,clean_data3-sigma3,'b--');
hold on;
plot(xvector,bayes_avg,'k-','LineWidth',2);
plot(xvector,bayes_avg+sigmabayes,'k--','LineWidth',1);
plot(xvector,bayes_avg-sigmabayes,'k--','LineWidth',1);
