clear;
clc;
close all;
global A thetas A0 c1 c2 c3 d1 d2 d3 Gamma gamma kp a w e1 e2 e3 k;

PRINT = true;
%PRINT = false;

%Simulation time
tfinal = 5;

% Unit vectors
e1 = [1 0 0]';
e2 = [0 1 0]';
e3 = [0 0 1]';

% System matrix
A = [0 1 0;0 0 1;0 0 0];

%Reference
a = [1 1];
w = [1 3];

%% First parameters

kp_1 = 5;
Z_1 = [1];
P_1 = [1 3 3 1];
thetas_1 = [kp_1 P_1(2) P_1(3) P_1(4)]';

k_1 = [1 1 1]';

%Initial conditions
X0_1  = [0 0 0]';
theta0_1 = [0 0 0 0]';
lambda0_1 = [0 0 0]';
eta0_1 = [0 0 0]';
rho0_1 = 1;

%Adaptation gain
Gamma_1 = 1;
gamma_1 = 1;
c1_1 = 1;
c2_1 = 1;
c3_1 = 1;
d1_1 = 1;
d2_1 = 1;
d3_1 = 1;
%% Second parameters

kp_2 = 5;
Z_2 = [1];
P_2 = [1 3 3 1];
thetas_2 = [kp_2 P_2(2) P_2(3) P_2(4)]';

k_2 = [1 1 1]';

%Initial conditions
X0_2  = [0 0 0]';
theta0_2 = [1 1 1 1]';
lambda0_2 = [0 0 0]';
eta0_2 = [0 0 0]';
rho0_2 = 1;

%Adaptation gain
Gamma_2 = 1;
gamma_2 = 1;
c1_2 = 1;
c2_2 = 1;
c3_2 = 1;
d1_2 = 1;
d2_2 = 1;
d3_2 = 1;