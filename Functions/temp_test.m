
% STEPS:
% Functionalize St_to_XY for each case
% Functionalize XY_to_St for each case

close all;

% Define the central traversal - a right angle downward turn
central_path = [-2 0; 2 0; 2 -2];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
stations = 3; %[1; 2; 2.5];


% Call the Ortho library to get orthogonal vectors
flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

fig_num = 14;  % Define the figure

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalTraversalVectorsAtStations(...
    stations,central_traversal,flag_rounding_type,fig_num);
title('Vertex projection via averaging everywhere (flag=4)');
axis equal;
axis square

unit_vectors = unit_normal_vector_end - unit_normal_vector_start;





% Pmid = [1 1];
% Pend = [3 -1];
% frac = 0.5;
% v_end = [2 0]; % Can be anything
% dP = 2;

Pmid = [0 0];
Pend = [2 0];
frac = 0.8765;
v_end = [1 1]; % Can be anything - converted to a unit vector below
dP = 3.24;


% Preliminary calculations to get unit vectors
mag_v_end = sum(v_end.^2,2).^0.5;
unit_v_end = v_end/mag_v_end;

v_Pmid_to_Pend = Pend-Pmid;
mag_v_Pmid_to_Pend = sum(v_Pmid_to_Pend.^2,2).^0.5;
unit_v_Pmid_to_Pend = v_Pmid_to_Pend/mag_v_Pmid_to_Pend;

v_mid = v_Pmid_to_Pend*[0 1; -1 0];
mag_v_mid = sum(v_mid.^2,2).^0.5;
unit_v_mid = v_mid/mag_v_mid;


% Calculate the position, P
% This is a test to show that averaging arbitrary vectors together, if not
% in unit form, gives a different vector than the average of unit vectors.
% So we have to be careful to use unit vectors
if 1==0
    v_bad = (1-frac)*v_mid + frac*v_end;
    mag_v_bad = sum(v_bad.^2,2).^0.5;
    unit_v_bad = v_bad/v_bad;

    v_unit = (1-frac)*unit_v_mid + frac*unit_v_end;
    mag_v_unit = sum(v_unit.^2,2).^0.5;
    unit_v_unit = v_unit/mag_v_unit;
end
v = (1-frac)*unit_v_mid + frac*unit_v_end;
mag_v = sum(v.^2,2).^0.5;
unit_v = v/mag_v;

fprintf(1,'Comparison of unit normal vectors:');
fprintf(1,'Unit vector from orthogonal projection at s = %0.2f:\n',stations)
disp(unit_vectors)
fprintf(1,'Unit vector calculated by hand: \n');
disp(unit_v);

P1 = unit_normal_vector_start + dP*unit_vectors;

Pbase = Pmid + frac*v_Pmid_to_Pend;
P = dP*unit_v + Pbase;

% More preliminary calculations
% Use cross-product of v_mid and v_end to find theta
cross_product = crossProduct(unit_v_end,unit_v_mid);
theta = asin(cross_product);
L = mag_v_Pmid_to_Pend;
H = L/tan(theta);

Px = P(1,1);
Py = P(1,2);
Vex = unit_v_end(1,1);
Vey = unit_v_end(1,2);
Vmx = unit_v_mid(1,1);
Vmy = unit_v_mid(1,2);

% Check that the intermediate calculations match
Px_calc = frac*L + Py*Vex*(frac/((1-frac)+frac*Vey));

fprintf(1,'Comparison of Px values:');
fprintf(1,'Px from original point:\n');
disp(Px);
fprintf(1,'Px from intermediate calculations: \n');
disp(Px_calc);

% Solve for s using quadratic formula
a = (L-L*Vey);
b = (Px*Vey -Px -L -Py*Vex);
c = Px;
f_calc1 = (-b + (b^2-4*a*c).^0.5)/(2*a);
f_calc2 = (-b - (b^2-4*a*c).^0.5)/(2*a);

fprintf(1,'Comparison of frac values:');
fprintf(1,'frac from original point:\n');
disp(frac);
fprintf(1,'frac from calculations: \n');
disp(f_calc2);

Vy_calc_1 = (1-f_calc2)+f_calc2*Vey;
Vx_calc_1 = frac*Vex; 

mag = (Vy_calc_1.^2 + Vx_calc_1.^2).^0.5;

Vy_calc = Vy_calc_1/mag;

d_calc = Py/Vy_calc;

fprintf(1,'Comparison of d values:');
fprintf(1,'d from original point:\n');
disp(dP);
fprintf(1,'d from calculations: \n');
disp(d_calc);

% Plot the results
figure(232);
clf;
hold on;
grid on;
axis equal
axis([-2 7 -2 7]);

points = [Pend; Pmid];
plot(points(:,1),points(:,2),'r.-','LineWidth',3,'MarkerSize',20);
text(Pmid(:,1),Pmid(:,2),'Pmid');
text(Pend(:,1),Pend(:,2),'Pend');
quiver(Pend(:,1),Pend(:,2),v_end(:,1),v_end(:,2),0,'r','LineWidth',2);
quiver(Pmid(:,1),Pmid(:,2),v_mid(:,1),v_mid(:,2),0,'r','LineWidth',2);
plot(P(:,1),P(:,2),'b.','LineWidth',3,'MarkerSize',20);
text(P(:,1),P(:,2),'P');
% Show the vectors
quiver(Pbase(:,1),Pbase(:,2),unit_v(:,1),unit_v(:,2),0,'b','LineWidth',2);
quiver(Pmid(:,1),Pmid(:,2),unit_v_mid(:,1),unit_v_mid(:,2),0,'b','LineWidth',2);
quiver(Pend(:,1),Pend(:,2),unit_v_end(:,1),unit_v_end(:,2),0,'b','LineWidth',2);

%% Now try with non-normal vectors
close all;

% Define the central traversal - a right angle downward turn
central_path = [4 4; 6 6; 8 4];
central_traversal = fcn_Path_convertPathToTraversalStructure(central_path);
stations = 2.5; %[1; 2; 2.5];


% Call the Ortho library to get orthogonal vectors
flag_rounding_type = 4;  % use average projection of prior and following segments always, with interpolation

fig_num = 15;  % Define the figure

% Calculate the unit normal vectors at given stations and put results into
% the figure.
[unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalTraversalVectorsAtStations(...
    stations,central_traversal,flag_rounding_type,fig_num);
title('Vertex projection via averaging everywhere (flag=4)');
axis equal;
axis square

unit_vectors = unit_normal_vector_end - unit_normal_vector_start;





% Pmid = [1 1];
% Pend = [3 -1];
% frac = 0.5;
% v_end = [2 0]; % Can be anything
% dP = 2;

Pmid = [5 5];
Pend = [6 6];
frac = (2.5-(2)^0.5)/(2^0.5);
v_end = [0 1]; % Can be anything - converted to a unit vector below
dP = 3.24;


% Preliminary calculations to get unit vectors
mag_v_end = sum(v_end.^2,2).^0.5;
unit_v_end = v_end/mag_v_end;

v_Pmid_to_Pend = Pend-Pmid;
mag_v_Pmid_to_Pend = sum(v_Pmid_to_Pend.^2,2).^0.5;
unit_v_Pmid_to_Pend = v_Pmid_to_Pend/mag_v_Pmid_to_Pend;

v_mid = v_Pmid_to_Pend*[0 1; -1 0];
mag_v_mid = sum(v_mid.^2,2).^0.5;
unit_v_mid = v_mid/mag_v_mid;


% Calculate the position, P
% This is a test to show that averaging arbitrary vectors together, if not
% in unit form, gives a different vector than the average of unit vectors.
% So we have to be careful to use unit vectors
if 1==0
    v_bad = (1-frac)*v_mid + frac*v_end;
    mag_v_bad = sum(v_bad.^2,2).^0.5;
    unit_v_bad = v_bad/v_bad;

    v_unit = (1-frac)*unit_v_mid + frac*unit_v_end;
    mag_v_unit = sum(v_unit.^2,2).^0.5;
    unit_v_unit = v_unit/mag_v_unit;
end
v = (1-frac)*unit_v_mid + frac*unit_v_end;
mag_v = sum(v.^2,2).^0.5;
unit_v = v/mag_v;

fprintf(1,'Comparison of unit normal vectors:');
fprintf(1,'Unit vector from orthogonal projection at s = %0.2f:\n',frac)
disp(unit_vectors)
fprintf(1,'Unit vector calculated by hand: \n');
disp(unit_v);

Pbase = Pmid + frac*v_Pmid_to_Pend;
P = dP*unit_v + Pbase;

% More preliminary calculations
% Use cross-product of v_mid and v_end to find theta
cross_product = crossProduct(unit_v_end,unit_v_mid);
theta = asin(cross_product);
L = mag_v_Pmid_to_Pend;
H = L/tan(theta);

Px = P(1,1);
Py = P(1,2);
Pmx = Pmid(1,1);
Pmy = Pmid(1,2);
Vex = unit_v_end(1,1);
Vey = unit_v_end(1,2);
Vmx = unit_v_mid(1,1);
Vmy = unit_v_mid(1,2);

% Check that the intermediate calculations match
Px_calc = Pmid(1,1) + Vmy*frac*L + dP*((1-frac)*Vmx + frac*Vex)/mag_v;
Py_calc = Pmid(1,2) - Vmx*frac*L + dP*((1-frac)*Vmy + frac*Vey)/mag_v;

Py_calc = Pmid(1,2) -Vmx*frac*L + (Px - Pmid(1,1) - Vmy*frac*L)/((1-frac)*Vmx + frac*Vex) * ((1-frac)*Vmy + frac*Vey);

fprintf(1,'Comparison of Py values:');
fprintf(1,'Py from original point:\n');
disp(Py);
fprintf(1,'Py from intermediate calculations: \n');
disp(Py_calc);

% Solve for s using quadratic formula
a = Vmx*L*(Vex-Vmx) +Vmy*L*(Vey-Vmy);
b = (Py-Pmy)*(Vex-Vmx) - (Px-Pmx)*(Vey-Vmy) + (Vmx^2+Vmy^2)*L;
c = Vmx*(Py-Pmy)-Vmy*(Px-Pmx);
f_calcs = [(-b + (b^2-4*a*c).^0.5)/(2*a); (-b - (b^2-4*a*c).^0.5)/(2*a)];
good_one = (f_calcs>=0).*(f_calcs<=1);
f_calc = f_calcs(good_one~=0);


fprintf(1,'Comparison of frac values:');
fprintf(1,'frac from original point:\n');
disp(frac);
fprintf(1,'frac from calculations: \n');
disp(f_calc);

Vx_calc_1 = (1-f_calc)*Vmx + f_calc*Vex;
Vy_calc_1 = (1-f_calc)*Vmy + f_calc*Vey;

mag = (Vy_calc_1.^2 + Vx_calc_1.^2).^0.5;

Vx_calc = ((1-f_calc)*Vmx + f_calc*Vex)/mag;

d_calc = ((Px-Pmx) - f_calc*Vmy*L)/Vx_calc;

fprintf(1,'Comparison of d values:');
fprintf(1,'d from original point:\n');
disp(dP);
fprintf(1,'d from calculations: \n');
disp(d_calc);

% Plot the results
figure(232);
clf;
hold on;
grid on;
axis equal
%axis([-2 7 -2 7]);

points = [Pend; Pmid];
plot(points(:,1),points(:,2),'r.-','LineWidth',3,'MarkerSize',20);
text(Pmid(:,1),Pmid(:,2),'Pmid');
text(Pend(:,1),Pend(:,2),'Pend');
quiver(Pend(:,1),Pend(:,2),v_end(:,1),v_end(:,2),0,'r','LineWidth',2);
quiver(Pmid(:,1),Pmid(:,2),v_mid(:,1),v_mid(:,2),0,'r','LineWidth',2);
plot(P(:,1),P(:,2),'b.','LineWidth',3,'MarkerSize',20);
text(P(:,1),P(:,2),'P');
% Show the vectors
quiver(Pbase(:,1),Pbase(:,2),unit_v(:,1),unit_v(:,2),0,'b','LineWidth',2);
quiver(Pmid(:,1),Pmid(:,2),unit_v_mid(:,1),unit_v_mid(:,2),0,'b','LineWidth',2);
quiver(Pend(:,1),Pend(:,2),unit_v_end(:,1),unit_v_end(:,2),0,'b','LineWidth',2);

%% Try again


% % Check Px formula (it works)
% Px_calc = frac*L + dP*unit_v(1,1);
% Py_calc = dP*unit_v(1,2);
% 
% % Check vx formula
% vx_calc = frac*Vxe;
% vy_calc = (1-frac)+frac*Vye;
% 



%% Now, test the reverse process 

s = Px/(1+Py/H);

%% Try again
% % ANOTHER BAD METHOD
% Vxe = v_end(1,1);
% Vye = v_end(1,2);
% Vxm = v_mid(1,1);
% Vym = v_mid(1,2);
% 
% % Use cross-product of v_mid and v_end to find theta
% cross_product = crossProduct(unit_v_end,unit_v_mid);
% theta = asin(cross_product);
% H = mag_v_Pmid_to_Pend/tan(theta);
% 
% % Solve for s using quadratic formula
% a = (Vye-Vym);
% b = (Vym+(Vxm-Vxe)* H );
% c = (-H*Vxm);
% s = (-b + (b^2-4*a*c).^0.5)/(2*a);

% % (BAD METHOD)
% % Step 1 - find dp, then project from Pmid to point A
% v_Pmid_to_P = P-Pmid;
% d_p = unit_v_mid*v_Pmid_to_P';
% A = unit_v_mid*d_p + Pmid;
% plot(A(:,1),A(:,2),'g.','LineWidth',3,'MarkerSize',20);
% text(A(:,1),A(:,2),'A');
% 
% % Step 2 - use cross-product of v_mid and v_end to find theta
% cross_product = crossProduct(unit_v_end,unit_v_mid);
% theta = asin(cross_product);
% 
% % Step 3 - use theta ad dp to find d_end, then point B
% d_end = d_p/cos(theta);
% B = unit_v_end*d_end + Pend;
% plot(B(:,1),B(:,2),'c.','LineWidth',3,'MarkerSize',20);
% text(B(:,1),B(:,2),'B');
% 
% % Step 4 - use the dot product of vector from A to P against vector A to B
% % to find d
% v_A_to_B = B - A;
% quiver(A(:,1),A(:,2),v_A_to_B(:,1),v_A_to_B(:,2),0,'c','LineWidth',2);
% v_A_to_P = P - A;
% quiver(A(:,1),A(:,2),v_A_to_P(:,1),v_A_to_P(:,2),0,'b','LineWidth',1);
% mag_v_A_to_B = sum(v_A_to_B.^2,2).^0.5;
% unit_v_A_to_B = v_A_to_B/mag_v_A_to_B;
% d = v_A_to_P*unit_v_A_to_B';
% 
% % Step 5 - use D to solve for the ratio, r, of s/S
% D = sum(v_A_to_B.^2,2).^0.5;
% r = d/D;
% fprintf(1,'Comparison of ratios: \n');
% disp([frac; r]);
% 
% 
% % Step 6 - use the ratio to solve for Pbase and dp
% Pbase_calc = r*v_Pmid_to_Pend + Pmid;
% fprintf(1,'Comparison of Pbase: \n');
% disp([Pbase; Pbase_calc])
% 
% dp = sum((Pbase_calc-P).^2,2).^0.5;
% fprintf(1,'Comparison of distaces: \n');
% disp([dP; dp])


%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§
%% Calculate cross products
function result = crossProduct(v,w)
result = v(:,1).*w(:,2)-v(:,2).*w(:,1);
end
