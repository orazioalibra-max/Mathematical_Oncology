%Extension of the model - continuous population

function PDEmodel()
    clc
    clear all
    close all
    % N=x(1) Tumor Volume                                    mm^3
    % K=x(2) Residual Healthy CTL Population (CD8+)          per k per mm^3
    % B=x(3) Progenitor Exhausted T cells                    per k per mm^3
    % P=x(4) Interleukin-2 Cytokine Concentration            per k per mm^3
    % E=x(5) Terminally Exhausted T cells                    per k per mm^3
    
    alpha=0.3; % Tumor growth rate                           per day
    Tc=4472.5; % Tumor carrying capacity                     mm^3
    lambda_K=0.111; % Tumor-CTL lysis parameter              mm^3 per k per day
    lambda_E=0.0888; % Tumor-Terminally Exhausted lys. param mm^3 per k per day
    hphy=10; % Magnifier for physical interactions           const
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
    r=1001;
    sx=0;
    dx=1;
    y=linspace(sx,dx,r)';
    h=(dx-sx)/(r-1);
    x0=[3,2*1e-4,0,0,zeros(1,r)];
    
    v_y=5*rE;%rE+4*y*rE;
    delta_y=delta_K+y.*(delta_E-delta_K);
    lambda_y=lambda_K+y.*(lambda_E-lambda_K);
    
    %First Block: PDE
    Apde=spdiags([v_y,-v_y]/h,-1:0,r,r)-spdiags(delta_y,0,r,r);
    Apde=[sparse(r,4),Apde];
    Apde(1,3)=rE/h;
    
    %Second Block: ODEs
    w=h*[0.5, ones(1, r-2), 0.5]; %weights

    opts = odeset(...
        'Jacobian',    @(t,x) Jacobian(t,x), ...
        'AbsTol',      1e-8, ...
        'RelTol',      1e-6, ...
        'NonNegative', 1:(4+r));
    
    [t_sol, X]=ode15s(@(t,x) ODE(t,x), tspan, x0, opts);

    %% Plot
    N=X(:,1);
    K=X(:,2);
    B=X(:,3);
    P=X(:,4);
    E=X(:,5:end);
    
    % Time dynamic ODE
    figure    
    subplot(2,2,1); plot(t_sol, N, 'k', 'LineWidth', 2); title('Tumor Volume (N)'); xlabel('Time (days)'); ylabel('Volume (mm^3)');
    subplot(2,2,2); plot(t_sol, K, 'b', 'LineWidth', 2); title('Healthy CTLs (K)'); xlabel('Time (days)'); ylabel('Population');
    subplot(2,2,3); plot(t_sol, B, 'g', 'LineWidth', 2); title('Progenitor Exhausted Cells (B)'); xlabel('Time (days)'); ylabel('Population');
    subplot(2,2,4); plot(t_sol, P, 'm', 'LineWidth', 2); title('IL-2 Cytokine (P)'); xlabel('Time (days)'); ylabel('Concentration');

    % Space-time evolution of E(t,y) (Heatmap)
    figure
    [Y_grid, T_grid]=meshgrid(y, t_sol);
    surf(T_grid, Y_grid, E, 'EdgeColor', 'none');
    view(2); %from above
    colormap parula; colorbar;
    title('Distribution Terminally Exhausted T cells E(y,t)');
    xlabel('Time (days)');
    ylabel('Exhaustion Stage (y)');
    zlabel('Density E');
    
    %% ODE
    function dxdt = ODE(~, x)
        Etot=w*x(5:end); %trapezoidal rule
        integral=(w.*lambda_y')*x(5:end); %trapezoidal rule
        
        N_num=lambda_K*x(1)*x(2)+x(1)*integral;
        N_den=x(1)+x(2)+Etot+1;
        kill=hphy*N_num/N_den;
        
        dN=alpha*x(1)*(1-(x(1)/Tc))-kill;
        dK=hphy*x(4)*x(2)-delta_K*x(2);
        dB=2^m*sB*x(1)+hphy*x(4)*x(3)-rE*x(3)-delta_B*x(3);
        dP=rho_K*x(2)+rho_B*x(3)-hphy*x(4)*(x(2)+x(3))-delta_P*x(4);
        dE=Apde*x;
        
        dxdt=[dN; dK; dB; dP; dE];
    end

    %% Jacobian
    function J = Jacobian(~,x)

        Etot=w* x(5:end); %trapezoidal rule
        integral=(w.*lambda_y')*x(5:end); %trapezoidal rule
        
        N_num=lambda_K*x(1)*x(2)+x(1)*integral;
        N_den=x(1)+x(2)+Etot+1;

        dkill_dN= hphy*((lambda_K*x(2)+integral)*N_den-N_num)/N_den^2;
        dkill_dK=hphy*(lambda_K*x(1)*N_den -N_num)/N_den^2;
        df_dE=x(1)*(lambda_y'.*w); 
        dg_dE=w;
        dkill_dE=hphy*(df_dE*N_den-N_num*dg_dE)/N_den^2;
        
        DN=[alpha*(1-2*x(1)/Tc)-dkill_dN,-dkill_dK,0,0,-dkill_dE];  
        DK=[0,hphy*x(4)-delta_K,0,hphy*x(2),zeros(1,r)];
        DB=[2^m*sB,0,hphy*x(4)-rE-delta_B,hphy*x(3),zeros(1,r)];
        DP=[0,rho_K-hphy*x(4),rho_B-hphy*x(4),-hphy*(x(2)+x(3))-delta_P,zeros(1,r)];
        
        Aode=[DN;DK;DB;DP];
        J = sparse([Aode; Apde]);
    end

end