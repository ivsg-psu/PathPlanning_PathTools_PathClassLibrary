% script to test the KF functionality

%% Create test data
dt = 0.01;
t = (0:dt:10)';
x = max(min(sin(t),0.7),-0.4);
xdot = [0;diff(x)]/dt;

% Add noise to state and/or derivative
x1_Sigma = 0.01*ones(length(t),1);
x1 = x + x1_Sigma.*randn(length(t),1);
x1_Sigma = 100*ones(length(t),1);

% Disturb x1 signficantly
bad_indices = find(t>4 & t<5);
x1(bad_indices) = x1(bad_indices)+0.2;



x1dot_Sigma = 0.01*ones(length(t),1);
x1dot = xdot + x1dot_Sigma.*randn(length(t),1);


% Show the signalw we created
figure(1);
clf;
hold on;
grid minor;
plot(t,x1,'r');
plot(t,x,'b','Linewidth',3);

figure(2);
clf;
hold on;
grid minor;
plot(t,x1dot,'r');
plot(t,xdot,'b','Linewidth',1);


t_x1 = t;
t_x1dot = t;

nameString = 'Test';
[x_kf,sigma_x] = fcn_KFmergeStateAndStateDerivative(t_x1,x1,x1_Sigma,t_x1dot,x1dot,x1dot_Sigma,nameString);

figure(1);
plot(t,x_kf,'g','Linewidth',1);