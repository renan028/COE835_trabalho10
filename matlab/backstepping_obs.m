%----------------------------------------------------------------------
%
%  COE-835  Controle adaptativo
%
%  Script para simular o trabalho 10
%
%  Backstepping  :  n  = 3     Second and third order plant
%                   n* = 3     Relative degree
%                   np = 4     Adaptive parameters
% Caso com observador completo
%----------------------------------------------------------------------

function dx = backstepping_obs(t,x)

global A thetas A0 Gamma gamma kp a w e1 e3 z1 alpha_bar tau3 u;

X           = x(1:3); y = e1'*X;
theta       = x(4:7);
lambda      = x(8:10);
eta         = x(11:13);
rho         = x(14);

%% Input
yr=0; dyr=0; ddyr=0;dddyr=0;
for i=1:length(a)
    yr = yr + a(i)*sin(w(i)*t);
    dyr = dyr + w(i)*a(i)*cos(w(i)*t);
    ddyr = ddyr - w(i)^2*a(i)*sin(w(i)*t);
    dddyr = dddyr - w(i)^3*a(i)*cos(w(i)*t);
end

Phi = [-y 0 0;0 -y 0;0 0 -y];

eta1 = eta(1);
eta2 = eta(2);
eta3 = eta(3);

theta1 = theta(1);
theta2 = theta(2);
theta3 = theta(3);
theta4 = theta(4);

l1 = lambda(1);
l2 = lambda(2);
l3 = lambda(3);


%% Z
z1_ = eval(subs(z1));
alpha_bar_ = eval(subs(alpha_bar));

%% Filtro eta
deta = A0*eta + e3*y;

%% Atualização 1
drho = - gamma * z1_ * sign(kp) * (dyr + alpha_bar_);

%% Variables 3
tau3_ = eval(subs(tau3));
dtheta = Gamma * tau3_;
u_ = eval(subs(u));


%% Filtros
dlambda = A0*lambda + e3*u_;

%% Planta
F = [e3*u_ Phi];
dX = A*X + F*thetas;

%% Translation
dx = [dX' dtheta' dlambda' deta' drho]';
