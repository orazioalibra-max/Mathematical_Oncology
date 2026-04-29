% [3] T. Simmons and D. Levy. Modeling the development of cellular exhaustion
% and tumor-immune stalemate. Bulletin of Mathematical Biology, 85(11):106, 2023.
%
% Data reproduction
clc
clear all
close all

%% Initial data

% N=x(1) Tumor Volume                                    mm^3
% K=x(2) Residual Healthy CTL Population (CD8+)          const
% B=x(3) Progenitor Exhausted T cells                    per k per mm^3
% E=x(4) Terminally Exhausted T cells                    per k per mm^2
% P=x(5) Interleukin-2 Cytokine Concentration            per k per mm^3

alpha=0.3; % Tumor growth rate                           per day
Tc=4472.5; % Tumor carrying capacity                     mm^3
lambda_K=0.111; % Tumor-CTL lysis parameter              mm^3 per k per day
lambda_E=0.0888; % Tumor-Terminally Exhausted lys. param mm^3 per k per day
h=10; % Magnifier for physical interactions              const
delta_K=0.4; % Death rate of CTLs                        per k per day
delta_B=0.41; % Death rate of Progenitor Exhausted cells per k per day   
delta_E=0.72; % Death rate of Terminally Exhausted cells per k per day
delta_P=5.5; % Decay rate of IL-2                        per k per day
m=7; % Progenitor cell divisions                         const
sB=0.0012; % Naive cell to Progenitor cell supply rate   per k per day
rE=0.0542; % Differentiation Rate                        per k per day
rho_K=1; % CTL cytokine secretion rate                   per k per day
rho_B=0.36; % Progenitor cytokine secretion              per k per day

tspan=linspace(0,80);

%% Paper model
odefun= @(t,x) [alpha*x(1)*(1-(x(1)/Tc))-h*((lambda_K*x(1)*x(2)+lambda_E*x(1)*x(4))/(x(1)+x(2)+x(4)+1));...
    h*x(5)*x(2)-delta_K*x(2);...
    2^m*sB*x(1)+h*x(5)*x(3)-rE*x(3)-delta_B*x(3);...
    rE*x(3)-delta_E*x(4);...
    rho_K*x(2)+rho_B*x(3)-h*x(5)*(x(2)+x(3))-delta_P*x(5)];

%x0triv=[0,0,0,0,0];
x0=[3,2*1e-4,0,0,0];

[~,x]=ode23(odefun,tspan, x0);

figure
subplot(5,1,1); plot(tspan, x(:,1)); title('N'); xlabel('Days'); ylabel('Tumor');
subplot(5,1,2); plot(tspan, x(:,2)); title('K'); xlabel('Days'); ylabel('Preexisting');
subplot(5,1,3); plot(tspan, x(:,3)); title('B'); xlabel('Days'); ylabel('Progenitor');
subplot(5,1,4); plot(tspan, x(:,4)); title('E'); xlabel('Days'); ylabel('Terminal');
subplot(5,1,5); plot(tspan, x(:,5)); title('P'); xlabel('Days'); ylabel('IL-2');
drawnow;


%% Extension of the model - homogeneous assumptions

% N=x(1) Tumor Volume                                    mm^3
% K=x(2) Residual Healthy CTL Population (CD8+)          const
% B=x(3) Progenitor Exhausted T cells                    per k per mm^3
% E_0=x(4) Initially Exhausted T cells                   per k per mm^2
% P=x(5) Interleukin-2 Cytokine Concentration            per k per mm^3
% E_half=x(6) Middle-Exausted T cells                    per k per mm^2
% E_1=x(7) Deeply-Exausted T cells                       per k per mm^2


lambda_y=[lambda_E,lambda_E,lambda_E];
delta_y=[delta_E,delta_E,delta_E];
r_t=5*rE;

odefun2= @(t,x) [alpha*x(1)*(1-(x(1)/Tc))-h*((lambda_K*x(1)*x(2)+lambda_y*x(1)*[x(4);x(6);x(7)])/(x(1)+x(2)+x(4)+x(6)+x(7)+1));...
    h*x(5)*x(2)-delta_K*x(2);...
    2^m*sB*x(1)+h*x(5)*x(3)-rE*x(3)-delta_B*x(3);...
    rE*x(3)-r_t*x(4)-delta_y(1)*x(4);...
    rho_K*x(2)+rho_B*x(3)-h*x(5)*(x(2)+x(3))-delta_P*x(5);...
    r_t*x(4)-r_t*x(6)-delta_y(2)*x(6);...
    r_t*x(6)-delta_y(3)*x(7)];

x0_E=[3,2*1e-4,0,0,0,0,0];

[~,xHom]=ode23(odefun2,tspan, x0_E);

figure
subplot(4,2,2); plot(tspan, xHom(:,1)); title({'Extended model','N'}); xlabel('Days'); ylabel('Tumor');
subplot(4,2,4); plot(tspan, xHom(:,2)); title('K'); xlabel('Days'); ylabel('Preexisting');
subplot(4,2,6); plot(tspan, xHom(:,3)); title('B'); xlabel('Days'); ylabel('Progenitor');
subplot(4,2,8); plot(tspan, xHom(:,5)); title('P'); xlabel('Days'); ylabel('IL-2');
subplot(4,2,1); plot(tspan, x(:,1)); title({'Paper model','N'}); xlabel('Days'); ylabel('Tumor');
subplot(4,2,3); plot(tspan, x(:,2)); title('K'); xlabel('Days'); ylabel('Preexisting');
subplot(4,2,5); plot(tspan, x(:,3)); title('B'); xlabel('Days'); ylabel('Progenitor');
subplot(4,2,7); plot(tspan, x(:,5)); title('P'); xlabel('Days'); ylabel('IL-2');

figure
subplot(3,2,1); plot(tspan, x(:,4)); title('E'); xlabel('Days'); ylabel('Terminal');
subplot(3,2,2); plot(tspan, xHom(:,4)); title('E0'); xlabel('Days'); ylabel('Initially exhausted');
subplot(3,2,3); plot(tspan, xHom(:,4)+xHom(:,6)+xHom(:,7)); title('E = E0 + E1/2 + E1'); xlabel('Days'); ylabel('Sum');
subplot(3,2,4); plot(tspan, xHom(:,6)); title('E1/2'); xlabel('Days'); ylabel('Middle exhausted');
subplot(3,2,6); plot(tspan, xHom(:,7)); title('E1'); xlabel('Days'); ylabel('Terminally exhausted');
drawnow;

%% Extension of the model - heterogeneous assumptions

% N=x(1) Tumor Volume                                    mm^3
% K=x(2) Residual Healthy CTL Population (CD8+)          const
% B=x(3) Progenitor Exhausted T cells                    per k per mm^3
% E_0=x(4) Initially Exhausted T cells                   per k per mm^2
% P=x(5) Interleukin-2 Cytokine Concentration            per k per mm^3
% E_half=x(6) Middle-Exausted T cells                    per k per mm^2
% E_1=x(7) Deeply-Exausted T cells                       per k per mm^2


lambda_y=[lambda_K,lambda_K+0.5*(lambda_E-lambda_K),lambda_E];
delta_y=[delta_K,delta_K+0.5*(delta_E-delta_K),delta_E];
r_t=5*rE;

odefun2= @(t,x) [alpha*x(1)*(1-(x(1)/Tc))-h*((lambda_K*x(1)*x(2)+lambda_y*x(1)*[x(4);x(6);x(7)])/(x(1)+x(2)+x(4)+x(6)+x(7)+1));...
    h*x(5)*x(2)-delta_K*x(2);...
    2^m*sB*x(1)+h*x(5)*x(3)-rE*x(3)-delta_B*x(3);...
    rE*x(3)-r_t*x(4)-delta_y(1)*x(4);...
    rho_K*x(2)+rho_B*x(3)-h*x(5)*(x(2)+x(3))-delta_P*x(5);...
    r_t*x(4)-r_t*x(6)-delta_y(2)*x(6);...
    r_t*x(6)-delta_y(3)*x(7)];

x0triv_E=[0,0,0,0,0,0,0];
x0_E=[3,2*1e-4,0,0,0,0,0];

[~,xHom]=ode23(odefun2,tspan, x0_E);

figure
subplot(4,2,2); plot(tspan, xHom(:,1)); title({'Extended model','N'}); xlabel('Days'); ylabel('Tumor');
subplot(4,2,4); plot(tspan, xHom(:,2)); title('K'); xlabel('Days'); ylabel('Preexisting');
subplot(4,2,6); plot(tspan, xHom(:,3)); title('B'); xlabel('Days'); ylabel('Progenitor');
subplot(4,2,8); plot(tspan, xHom(:,5)); title('P'); xlabel('Days'); ylabel('IL-2');
subplot(4,2,1); plot(tspan, x(:,1)); title({'Paper model','N'}); xlabel('Days'); ylabel('Tumor');
subplot(4,2,3); plot(tspan, x(:,2)); title('K'); xlabel('Days'); ylabel('Preexisting');
subplot(4,2,5); plot(tspan, x(:,3)); title('B'); xlabel('Days'); ylabel('Progenitor');
subplot(4,2,7); plot(tspan, x(:,5)); title('P'); xlabel('Days'); ylabel('IL-2');

figure
subplot(3,2,1); plot(tspan, x(:,4)); title('E'); xlabel('Days'); ylabel('Terminal');
subplot(3,2,2); plot(tspan, xHom(:,4)); title('E0'); xlabel('Days'); ylabel('Initially exhausted');
subplot(3,2,3); plot(tspan, xHom(:,4)+xHom(:,6)+xHom(:,7)); title('E0 + E1/2 + E1'); xlabel('Days'); ylabel('sum');
subplot(3,2,4); plot(tspan, xHom(:,6)); title('E1/2'); xlabel('Days'); ylabel('Middle exhausted');
subplot(3,2,6); plot(tspan, xHom(:,7)); title('E1'); xlabel('Days'); ylabel('Terminally exhausted');
drawnow;