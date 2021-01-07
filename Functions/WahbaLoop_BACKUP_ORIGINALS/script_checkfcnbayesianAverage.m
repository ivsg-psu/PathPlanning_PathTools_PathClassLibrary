%% To check
% function [bayes_avg, bayes_sigma] = fcn_bayesianAverage(value1, sigma1, value2, simgma2, varargin)
% varargin can take many arguments, arguments follow the order 
% value1, sigma1, value2, simag2 .......
% where 'value%d' is the parameter and 'sigma%d' is the corresponding std-deviation
% 'bayes_avg' is the weighted average
% 'bayes_sigma' is the weighted std-dev

clear all
close all
clc

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
