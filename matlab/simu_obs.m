%----------------------------------------------------------------------
%
%  COE-835  Controle adaptativo
%
%  Script para simular exemplo 
%
%  Backstepping  :  n  = 2     Second and third order plant
%                   n* = 2     Relative degree
%                   np = 4, 6     Adaptive parameters
% Com observador completo
%----------------------------------------------------------------------

global A thetas A0 c1 c2 c3 d1 d2 d3 Gamma gamma kp a w e1 e2 e3 k;

sim_str = strcat('sim0','_');

%% ------------------------------------------------- Simulation 1 (default)
kp = kp_1;
Z = Z_1;
P = P_1;

thetas = thetas_1;

k = k_1;
A0 = A - k*e1';

c1 = c1_1;
c2 = c2_1;
c3 = c3_1;
d1 = d1_1;
d2 = d2_1;
d3 = d3_1;
Gamma = Gamma_1*eye(4);
gamma = gamma_1;

% Initialization
X0  = X0_1;
theta0 = theta0_1;
lambda0 = lambda0_1;
eta0 = eta0_1;
rho0 = rho0_1;
init = [X0' theta0' lambda0' eta0' rho0]';

[T_1,X_1] = ode23s('backstepping_obs',tfinal,init,'');
y_1      = X_1(:,1);
theta_1 =  X_1(:,4:7);
tiltheta_1 = thetas' - theta_1;
modtt_1 = sqrt(sum(theta_1.^2,2));
r_1 = 0;
for i=1:length(a)
    r_1 = r_1 + a(i)*sin(w(i).*T_1);
end
e0_1 =  y_1 - r_1;


%% --------------------------------------------------- Simulation 2 (gamma)
changed = 1;

kp = kp_1;
Z = Z_1;
P = P_1;

thetas = thetas_1;

k = k_1;
A0 = A - k*e1';

c1 = c1_1;
c2 = c2_1;
c3 = c3_1;
d1 = d1_1;
d2 = d2_1;
d2 = c3_1;
Gamma = Gamma_2*eye(4);
gamma = gamma_2;

% Initialization
X0  = X0_1;
theta0 = theta0_1;
lambda0 = lambda0_1;
eta0 = eta0_1;
rho0 = rho0_1;
init = [X0' theta0' lambda0' eta0' rho0]';

[T_2,X_2] = ode23s('backstepping_obs',tfinal,init,'');
y_2      = X_2(:,1);
theta_2 =  X_2(:,4:7);
tiltheta_2 = thetas' - theta_2;
modtt_2 = sqrt(sum(theta_2.^2,2));
r_2 = 0;
for i=1:length(a)
    r_2 = r_2 + a(i)*sin(w(i).*T_2);
end
e0_2 =  y_2 - r_2;

%Plot
run plot_bkst.m;

% %% -------------------------------------------------- Simulation 3 (planta)
% changed = 2;
% 
% kp = kp_2;
% Z = Z_2;
% P = P_2;
% 
% thetas = thetas_2;
% 
% k = k_1;
% A0 = A - k*e1';
% 
% c1 = c1_1;
% c2 = c2_1;
% d1 = d1_1;
% d2 = d2_1;
% Gamma = Gamma_1*eye(3);
% gamma = gamma_1;
% 
% % Initialization
% X0  = X0_1;
% theta0 = theta0_1;
% lambda0 = lambda0_1;
% eta0 = eta0_1;
% rho0 = rho0_1;
% init = [X0' theta0' lambda0' eta0' rho0]';
% 
% [T_2,X_2] = ode23s('backstepping_obs',tfinal,init,'');
% y_2      = X_2(:,1);
% theta_2 =  X_2(:,3:5);
% tiltheta_2 = thetas' - theta_2;
% modtt_2 = sqrt(sum(theta_2.^2,2));
% r_2 = 0;
% for i=1:length(a)
%     r_2 = r_2 + a(i)*sin(w(i).*T_2);
% end
% e0_2 =  y_2 - r_2;
% 
% %Plot
% run plot_bkst.m;
% 
% %% ------------------------------------------------------ Simulation 4 (y0)
% changed = 3;
% 
% kp = kp_1;
% Z = Z_1;
% P = P_1;
% 
% thetas = thetas_1;
% 
% k = k_1;
% A0 = A - k*e1';
% 
% c1 = c1_1;
% c2 = c2_1;
% d1 = d1_1;
% d2 = d2_1;
% Gamma = Gamma_1*eye(3);
% gamma = gamma_1;
% 
% % Initialization
% X0  = X0_2;
% theta0 = theta0_2;
% lambda0 = lambda0_2;
% eta0 = eta0_2;
% rho0 = rho0_2;
% init = [X0' theta0' lambda0' eta0' rho0]';
% 
% [T_2,X_2] = ode23s('backstepping_obs',tfinal,init,'');
% y_2      = X_2(:,1);
% theta_2 =  X_2(:,3:5);
% tiltheta_2 = thetas' - theta_2;
% modtt_2 = sqrt(sum(theta_2.^2,2));
% 
% r_2 = 0;
% for i=1:length(a)
%     r_2 = r_2 + a(i)*sin(w(i).*T_2);
% end
% e0_2 =  y_2 - r_2;
% 
% %Plot
% run plot_bkst.m;

% %% ------------------------------------------------- Simulation 5 (model)
% changed = 4;
% 
% kp = kp_1;
% Z = Z_1;
% P = P_1;
% 
% thetas = thetas_2;
% 
% k = k_1;
% A0 = A - k*e1';
% 
% c1 = c1_1;
% c2 = c2_1;
% d1 = d1_1;
% d2 = d2_1;
% Gamma = Gamma_1*eye(3);
% gamma = gamma_1;
% 
% % Initialization
% X0  = X0_1;
% theta0 = theta0_1;
% lambda0 = lambda0_1;
% eta0 = eta0_1;
% rho0 = rho0_1;
% init = [X0' theta0' lambda0' eta0' rho0]';
% 
% [T_2,X_2] = ode23s('backstepping_obs',tfinal,init,'');
% y_2      = X_2(:,1);
% theta_2 =  X_2(:,3:5);
% tiltheta_2 = thetas' - theta_2;
% modtt_2 = sqrt(sum(theta_2.^2,2));
% r_2 = 0;
% for i=1:length(a)
%     r_2 = r_2 + a(i)*sin(w(i).*T_2);
% end
% e0_2 =  y_2 - r_2;
% 
% %Plot
% run plot_bkst.m;
