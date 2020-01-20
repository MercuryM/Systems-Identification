a = 1/2;
lamda2 = 9;
miu2 = 1;
var = miu2 + lamda2/(1 - a^2)
gamma0 = var
gamma1 = a * lamda2 /(1 - a^2)
gamma2 = a^2 * lamda2 /(1 - a^2)
theta1 = gamma1/gamma0
var1 = (1 + theta1^2) * gamma0 - 2 * theta1 * gamma1
theta2 = [(gamma1 * gamma0 - gamma2 * gamma1)/(gamma0^2 - gamma1^2); (gamma2 * gamma0 - gamma1^2)/(gamma0^2 - gamma1^2)]
var2 = (1 + theta2(1)^2 + theta2(2)^2) * gamma0 + (2 * theta2(1) * theta2(2) - 2 * theta2(1)) * gamma1 - 2 * theta2(2) * gamma2

% figure (1)
% x = 1:1:2
% yx = [2,2]
% plot (x,yx,'LineWidth',1.5)
% xlabel('\itN')
% ylabel('\theta_0')
% title('Simulated \theta_1^o versus \itN')

pv1 = 10.3308;
pv2 = 10.3123;

FPE1 = pv1 * (2000 + 1)/ (2000 - 1)
FPE2 = pv2 * (2000 + 2)/ (2000 - 2)
AIC1 = 2 * 1 / 2000 + log(pv1)
AIC2 = 2 * 2 / 2000 + log(pv2)
MDL1 = log(2000) / 2000 + log(pv1)
MDL2 = 2* log (2000) / 2000 + log(pv2)
