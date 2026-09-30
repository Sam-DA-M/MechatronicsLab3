clc
clear

%% Mass of each link (in kg)

M_cover = .583;
SamMass = 80.74;
M_small_lever = 0.0214;
M_BigLeverAndSpring3 = 0.10581;
M_rack = 0.009;
M_Lpiece = .001;
M_PaperWheel = .01014;
M_pinion = 0.007;

%% Distances and Lengths

L1 = 0.132;
L2 = .105;
L3 = 0.230;
L4 = .009;
L5 = 0.017;
L6 = L1*.25; %Approximation based on visual observation (Of length 6 relative to length 1)
L7 = L3*.15; %Of L3 relative to L7 (See above)
Lrack = .136;
R_pinion = 0.003;
R_PaperWheel = .063;

%% Dimensions of Springs

d_Spring1 = 0.001;
d_Spring3 = .0025;
d_Spring4 = .0003;

D_Spring1 = .008;
D_Spring3 = .018;
D_Spring4 = .006;

L_Spring1 = .01;
L_Spring3 = .016;
L_Spring4 = .0166;

%% MMoI
%Parallel axis theorem 
%Assumes that L-Piece is a perfect rectangle

MMoI_Biglever = (1/12)*(M_BigLeverAndSpring3*.5)*(L3^2) + (M_BigLeverAndSpring3*.5*(L7)^2);
MMoI_small_lever = (1/12)*M_small_lever*(L1^2) + (M_small_lever*(L6)^2);
MMoI_PaperWheel = .5*M_PaperWheel*(R_PaperWheel)^2;
MMoI_pinion = .5*M_pinion*(R_pinion)^2;
MMoI_L_piece = (1/12)*(M_Lpiece)*(L5^2+L4^2)+(M_Lpiece)*((sqrt(L4^2+L5^2)/2)^2); %Derived by treating L Piece as a rectangular plate


%% No. Coils of Springs
N1 = 9;
N3 = 3.75;
N4 = 40;
%% Calculating Spring Stiffnesses (N/m)

%G = E/(2(1+v)) where G is shear modulus, E is young's modulus, and v is
%Poisson's ratio of AISI 1075
%
E1 = 200*10^9;%Pa or N/m^2
E3 = E1;
E4 = E3;

v1 = .285;
v3 = v1;
v4 = v3;

G1 = E1/(2*(1+v1));
G3 = E1/(2*(1+v1));
G4 = E1/(2*(1+v1));

k_spring1 = ((d_Spring1)^4*G1)/(8*(D_Spring1)^3*N1); 
k_spring3 = ((d_Spring3)^4*G3)/(8*(D_Spring3)^3*N3);
k_spring4 = ((d_Spring4)^4*G4)/(8*(D_Spring4)^3*N4);

%% M(theta)''+ D(theta)' K(theta) = 0 



%% Meq Calculations

Meq1 = MMoI_small_lever/L1^2;
Jeq1 = MMoI_Biglever + Meq1*L2^2;
Meq2 = 2*(Jeq1/L3^2);
Jeq2 = Meq2*L4^2;
Meq3 = (Jeq2)/L5^2+M_rack;
Jeq3 = Meq3*R_pinion^2 + MMoI_PaperWheel + MMoI_pinion;

Jeq = Jeq3;

%% Keq calculations

Keq1rot = .5*k_spring1*L6^2;
Keq2trans = Keq1rot/L1^2;
Keq3rot = Keq2trans*L2^2 + (k_spring1*.5)*L7^2;
Keq4trans = 2*Keq3rot/(L3)^2 + k_spring3;
Keq5rot = Keq4trans*L4^2;
Keq6trans = Keq5rot/L5^2 - k_spring4;
Keq7rot = Keq6trans*(R_pinion)^2;

Keq = Keq7rot; %at pinion

%% Solving for Feq

Fin = .25*SamMass*9.81 + .25*M_cover*9.81; 
Teq1 = Fin*L6;
Feq1 = Teq1/L1;
Teq2 = Feq1*(L2) + Fin*L7;
Feq2 = Teq2/L3; 
Teq3 = 2*Feq2*L4; 
Feq3 = Teq3/L5;
Teq4 = Feq3*R_pinion;

Teq = Teq4; % Lumping at pinion


%% Deq calculations
% Assumptions:
% Coefficient of friction steel/steel =0.6
% Friction force is greatest and thus most significant at initial contact
% Normal force in each link is 10% of the initial normal force divided by 4
% D1 = coef_friction * Force_normal * 4 * 1/pi
% assuming input frequency and amplitude of oscillation to be 1


B1 = 4*(1/pi)*0.6*(SamMass*9.81/4)*L6^2; 
B2 = 4*(1/pi)*0.6*(SamMass*9.81/4)*L6^2; 
Deq1 = B1/L1^2; 
Beq1 = Deq1*(L2)^2+B2;  
Deq2 = (Beq1)/L3^2 + 4*(1/pi)*0.6*(Feq3-Feq2);
Beq2 = 2*Deq2*L4^2; %Account for both sides
Deq3 = Beq2/L5^2 + 4*(1/pi)*0.6*(Feq3);

Beq = Deq3*R_pinion^2;  %at pinion

%% Plot

% x = theta of disk
% x'=diff(x,t)
%Values of different parameters
J = Jeq %kg*m^2
K = Keq %N*m
B = Beq %N.s/m
T = Teq %N*m
oscilationfreq = sqrt(K/J)

syms x(t) %x(t) is is the x in Mx'' + Bx' + Kx = F

eqn1= J*diff(x,t,2) + B*diff(x,t) + K*x == T;
%Note: in the above equation, x'' is written as diff(x,t,2)
%x'is diff(x,t)

Dx=diff(x,t);

%specifyiung initial conditions
%If my weight is 178 lb, and (2pi rad)/300 lb represents a full rotation,
%Then Theta = 3.728 rad.
%therefore, look for angular position 'x' = 3.728 rad
initialCon=[x(0) == 0, Dx(0)==0];

%solving fox X(t)
solutionX=dsolve(eqn1,initialCon);

%Solving for V(t) - Velocity as a function of time
solutionV=diff(solutionX);

%Solving for A(t) - Acceleration as a function of time
solutionA=diff(solutionV);

%Plotting solution between times 0 and 10s
%SUCCESS!!! The Final Angular Position is about just a tenth of a radian off from that which we predicted. 
figure;
fplot(solutionX,[0 1]);
title 'Angular Position theta vs time';
xlabel 'time (s)';
ylabel 'Angular Position (rad)';

figure;
fplot(solutionV,[0 1]);
title 'Angular Velocity vs time';
xlabel 'time (s)';
ylabel 'Angular Velocity (rad/s)';

figure;
fplot(solutionA,[0 1]);
title 'Angular Acceleration vs time';
xlabel 'time (s)';
ylabel 'Angular Acceleration (rad/s^2)';


%laplace transforms