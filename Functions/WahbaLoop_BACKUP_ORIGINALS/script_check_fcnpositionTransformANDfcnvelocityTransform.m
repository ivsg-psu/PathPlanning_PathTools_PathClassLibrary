clear all
clc
close all

%% Check fcn_positionTransform
x = -250:1:250;
pointArray = [x', (x.^2)', abs(x)'];
angularPos = [pi*randn(length(x),1), pi*randn(length(x),1), pi*randn(length(x),1)];
pos = 0.0254*[21, 0, 0];

pointsTrans = fcn_positionTransform(pointArray, angularPos, pos);
%pointsTrans = fcn_positionTransform(pointArray, angularPos, pos,'deg');

figure(10001)
scatter3(pointArray(:,1), pointArray(:,2), pointArray(:,3));
figure(10002)
scatter3(pointsTrans(:,1), pointsTrans(:,2), pointsTrans(:,3));
%% Check fcn_velocityTransform
v = 0:0.1:100;
velArray = [v', (v.^2)', abs(v)'];
angularPos = [pi*randn(length(v),1), pi*randn(length(v),1), pi*randn(length(v),1)];
angularVel = [randn(length(v),1), randn(length(v),1), randn(length(v),1)];
pos = 0.0254*[21, 0, 0];



velTrans = fcn_velocityTransform(velArray, angularPos, angularVel, pos);
%velTrans = fcn_velocityTransform(velArray, angularPos, angularVel, pos, 'deg');

figure(10003)
scatter3(velArray(:,1), velArray(:,2), velArray(:,3));
figure(10004)
scatter3(velTrans(:,1), velTrans(:,2), velTrans(:,3));