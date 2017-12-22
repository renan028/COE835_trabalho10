global z1 z2 alpha_bar tau3 u
syms y yr dyr rho l1 l2 l3 ddyr dddyr real
eta = sym('eta', [3 1], 'real');
theta = sym('theta',[4 1], 'real');

xi = -A0^3*eta;
Xi = -[A0^2*eta A0*eta eta];
omega_bar = [0, (Xi(2,:) - y*e1')]';
omega = [l2, (Xi(2,:) - y*e1')]';

deta = A0*eta + e3*y;

z1 = y - yr;
alpha_bar = -c1*z1 - d1*z1 - xi(2) - omega_bar'*theta;
alpha_1 = rho * alpha_bar;
z2 = l2 - rho*dyr - alpha_1;

tau_1 = (omega - rho*(dyr + alpha_bar)*[e1',0]')*z1;
tau2 = tau_1 - z2 * (diff(alpha_1,y) * omega);

drho = - gamma * z1 * sign(kp) * (dyr + alpha_bar);

dadeta = [diff(alpha_1,eta(1)),diff(alpha_1,eta(2)),diff(alpha_1,eta(3))];
beta_2 = k(2)*l1 + diff(alpha_1,y) * (xi(2) + omega'*theta) + ...
    dadeta*deta + diff(alpha_1,yr) * dyr + ...
    (dyr + diff(alpha_1,rho)) * drho;

dadtheta = [diff(alpha_1,theta(1)),diff(alpha_1,theta(2)), ...
    diff(alpha_1,theta(3)), diff(alpha_1,theta(4))];
alpha2 = -c2*z2 + beta_2 + dadtheta*Gamma*tau2 - ...
    d2*z2*(diff(alpha_1,y))^2 - z1*theta(1);

z3 = l3 - rho*ddyr - alpha2;
tau3 = tau2 - diff(alpha2,y)*omega*z3;

%% beta3
da2deta = [diff(alpha2,eta(1)),diff(alpha2,eta(2)),diff(alpha2,eta(3))];
beta3 = k(3)*l1 + diff(alpha2,y) * (xi(2) + omega'*theta) + ...
    da2deta*deta + diff(alpha2,yr) * dyr + ...
    (ddyr + diff(alpha2,rho)) * drho + ...
    diff(alpha2,l1) * (-k(1)*l1+l2) + ...
    diff(alpha2,l2) * (-k(2)*l2+l3) + ...
    diff(alpha2,dyr) * ddyr;

%% u
da2dtheta = [diff(alpha2,theta(1)),diff(alpha2,theta(2)), ...
    diff(alpha2,theta(3)), diff(alpha2,theta(4))];
u = -c3*z3 + beta3 + rho*dddyr + da2dtheta*Gamma*tau3 - ...
    d3*(diff(alpha2,y))^2*z3 - ...
    z2*dadtheta*Gamma*diff(alpha2,y)*omega;


clear deta drho