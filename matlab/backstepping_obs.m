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

global A thetas A0 c1 c2 c3 d1 d2 d3 Gamma gamma kp a w e1 e2 e3 k;

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

%% Variables 1
xi = -A0^3 * eta;
Xi = -[A0^2*eta A0*eta eta];
v0_1 = lambda(1);
v0_2 = lambda(2);
v0_3 = lambda(3);
omega_bar = [0, (Xi(2,:) - y*e1')]';
omega = [v0_2, (Xi(2,:) - y*e1')]';

%% Z
z1 = y - yr;
alpha_bar = -c1*z1 - d1*z1 - xi(2) - omega_bar'*theta;
alpha_1 = rho * alpha_bar;
z2 = v0_2 - rho*dyr - alpha_1;

%% Filtro eta
deta = A0*eta + e3*y;

%% dalpha/dt
dady = rho * (- c1 - d1 + [0,e1']*theta);
dadeta_deta = rho * (e2' * A0^3 * deta + ...
    [0,e2'*A0^2*deta, e2'*A0*deta, e2'*eye(3)*deta]*theta);
dadyr = rho*(c1 + d1);
dadtheta = - rho * omega_bar';
dadrho = -(c1 + d1)*z1 - e2'*xi - omega_bar'*theta;

%% dz2/dt
dz2dy = - dady;
dz2deta_deta = -dadeta_deta;
dz2dyr = - dadyr;
dz2dtheta = - dadtheta;
dz2drho = - dyr - dadrho;
dz2dlambda1 = 0;
dz2dlambda2 = 1;
dz2ddyr = - rho;

%% dz1/dt
dz1dy = 1;
dz1deta_deta = 0;
dz1dyr = -1;
dz1dtheta = 0;
dz1drho = 0;
dz1dlambda1 = 0;
dz1dlambda2 = 0;
dz1ddyr = 0;

%% ddrho/dt
ddrhody = - gamma*sign(kp)*((dyr+alpha_bar) + z1*dady/rho); 
ddrhodeta = - gamma*sign(kp)*z1*dadeta_deta/rho;
ddrhodyr = - gamma*sign(kp)*(-1*(dyr+alpha_bar) + dadyr/rho);
ddrhodtheta = - gamma*sign(kp)*z1*dadtheta/rho;
ddrhodrho = 0;
ddrhodlambda1 = 0;
ddrhodlambda2 = 0;
ddrhoddyr = - gamma*sign(kp)*z1;

%% Variables 2
tau_1 = (omega - rho*(dyr + alpha_bar)*[e1',0]')*z1;
tau2 = tau_1 - z2 * (dady * omega);

%% Atualização 1
drho = - gamma * z1 * sign(kp) * (dyr + alpha_bar);
beta_2 = k(2)*v0_1 + dady * (xi(2) + omega'*theta) + ...
    dadeta_deta + dadyr * dyr + (dyr + dadrho) * drho;
alpha2 = -c2*z2 + beta_2 + dadtheta*Gamma*tau2 - d2*z2*(dady)^2 - ...
    z1*theta(1);

%% dbeta2/dt
dadeta_dy = rho * (e2' * A0^3 * e3 + ... 
    [0,e2'*A0^2*e3, e2'*A0*e3, e2'*eye(3)*e3]*theta);
dbeta2dy = dady*([0,e1']*theta) + dadeta_dy + ...
    drho*(-(c1+d1)-[0,e1']*theta) + (dyr+dadrho)*ddrhody;
dbeta2deta_deta = dady*dadeta_deta/rho + ...
    [0,e2'*A0^3*deta, e2'*A0^2*deta, e2'*A0*deta]*theta + ...
    drho*(-dadeta_deta/rho) + (dyr + dadrho)*ddrhodeta;
dbeta2dyr = drho*(c1+d1) + (dyr + dadrho)*ddrhodyr;
dbeta2dlambda1 = k(2);
dbeta2dlambda2 = dady*[1, 0, 0, 0]*theta;
dbeta2drho = (dady/rho)*(xi(2)+omega'*theta) + dadeta_deta/rho + ...
    (dadyr/rho)*dyr;
dbeta2ddyr = dadyr + drho + (dyr+dadrho)*ddrhoddyr;
dbeta2dtheta = rho*[0,e1']*(xi(2)+omega'*theta) + dady*omega' + ...
    rho*[0,e2'*A0^2*deta, e2'*A0*deta, e2'*eye(3)*deta] - ...
    drho*omega_bar' + (dyr+dadrho)*ddrhodtheta;

%% dtau2/dt
dtau2dy = [0,e1']'*z1 + omega - dady*[e1',0]'*z1 - ...
    rho*(dyr+alpha_bar)*[e1',0]' - z2*dady*[0,e1']' + (dady)^2*omega;
dtau2deta_deta = [0,e2'*A0^2*deta, e2'*A0*deta, e2'*eye(3)*deta]'*z1 - ...
    dadeta_deta*[e1',0]'*z1 + dadeta_deta*dady*omega - ...
    z2*dady*[0,e2'*A0^2*deta, e2'*A0*deta, e2'*eye(3)*deta]';
dtau2dyr = -omega - dadyr*[e1',0]'*z1 + rho*(dyr+alpha_bar)*[e1',0]' + ...
    dadyr*dady*omega;
dtau2dlambda1 = [1,0,0,0]'*z1 - z2*dady*[1,0,0,0]';
dtau2dlambda2 = -dady*omega;
dtau2drho = -(dyr+alpha_bar)*[e1',0]'*z1 - (-dyr-dadrho)*dady*omega - ...
    z2*(dady)/rho*omega;
dtau2ddyr = -rho*[e1',0]'*z1 + rho*dady*omega;
dtau2dtheta = 0;

%% dadtheta/dt
dadthetady = [0,0,0,0];
dadthetadeta_deta = -rho*[0,e2'*A0^2*deta, e2'*A0*deta, e2'*eye(3)*deta];
dadthetadyr = [0,0,0,0];
dadthetadlambda1 = [0,0,0,0];
dadthetadlambda2 = [0,0,0,0];
dadthetadrho = -[0,e2'*A0^2*eta, e2'*A0*eta, e2'*eye(3)*eta];
dadthetaddyr = [0,0,0,0];
dadthetadtheta = [0,0,0,0];

%% (dady)^2/dt
dady2dy = 0;
dady2deta_deta = 0;
dady2dyr = 0;
dady2dlambda1 = 0;
dady2dlambda2 = 0;
dady2drho = 2*rho*(- c1 - d1 + [0,e1']*theta)^2;
dady2ddyr = 0;
dady2dtheta = 2*(rho^2)*(- c1 - d1 + [0,e1']*theta)*[0,e1'];

%% da2/dt
da2dy = -c2*dz2dy - theta(1)*dz1dy + dbeta2dy + dadthetady*Gamma*tau2 + ...
    dadtheta*Gamma*dtau2dy - d2*dady2dy*z2 - d2*dady^2*dz2dy;
da2deta_deta = -c2*dz2deta_deta - theta(1)*dz1deta_deta + ...
    dbeta2deta_deta + dadthetadeta_deta*Gamma*tau2 + ...
    dadtheta*Gamma*dtau2deta_deta - d2*dady2deta_deta*z2 - ...
    d2*dady^2*dz2deta_deta;
da2dlambda1 = -c2*dz2dlambda1 - theta(1)*dz1dlambda1 + dbeta2dlambda1 + ...
    dadthetadlambda1*Gamma*tau2 + dadtheta*Gamma*dtau2dlambda1 - ...
    d2*dady2dlambda1*z2 - d2*dady^2*dz2dlambda1;
da2dlambda2 = -c2*dz2dlambda2 - theta(1)*dz1dlambda2 + dbeta2dlambda2 + ...
    dadthetadlambda2*Gamma*tau2 + dadtheta*Gamma*dtau2dlambda2 - ...
    d2*dady2dlambda2*z2 - d2*dady^2*dz2dlambda2;
da2dyr = -c2*dz2dyr - theta(1)*dz1dyr + dbeta2dyr + ...
    dadthetadyr*Gamma*tau2 + dadtheta*Gamma*dtau2dyr - ...
    d2*dady2dyr*z2 - d2*dady^2*dz2dyr;
da2ddyr = -c2*dz2ddyr - theta(1)*dz1ddyr + dbeta2ddyr + ...
    dadthetaddyr*Gamma*tau2 + dadtheta*Gamma*dtau2ddyr - ...
    d2*dady2ddyr*z2 - d2*dady^2*dz2ddyr;
da2drho = -c2*dz2drho - theta(1)*dz1drho + dbeta2drho + ...
    dadthetadrho*Gamma*tau2 + dadtheta*Gamma*dtau2drho - ...
    d2*dady2drho*z2 - d2*dady^2*dz2drho;
da2dtheta = -c2*dz2dtheta - theta(1)*dz1dtheta + dbeta2dtheta + ...
    dadthetadtheta*Gamma*tau2 + dadtheta*Gamma*dtau2dtheta - ...
    d2*dady2dtheta*z2 - d2*dady^2*dz2dtheta;

%% Variables 3
z3 = lambda(3) - rho*ddyr - alpha2;
tau3 = tau2 - da2dy*omega*z3;
dtheta = Gamma * tau3;
beta3 = k(3)*v0_1 + da2dy * (xi(2) + omega'*theta) + ...
    da2deta_deta + da2dyr * dyr + (ddyr + da2drho) * drho + ...
    da2dlambda1 * (-k(1)*lambda(1)+lambda(2)) + ...
    da2dlambda2 * (-k(2)*lambda(2)+lambda(3)) + da2ddyr*ddyr;
u = -c3*z3 + beta3 + rho*dddyr + da2dtheta*Gamma*tau3 - ...
    d3*(da2dy)^2*z3 - z2*dadtheta*Gamma*da2dy*omega;


%% Filtros
dlambda = A0*lambda + e3*u;

%% Planta
F = [e3*u Phi];
dX = A*X + F*thetas;

%% Translation
dx = [dX' dtheta' dlambda' deta' drho]';
