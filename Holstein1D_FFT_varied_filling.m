%% Holstein1D_FFT
%% Set LaTeX as default text interpreter
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
maxNumCompThreads(16);
%% Set System Parameters
N_k = 16; % number of sites
N_w_init = 50; % (minimum number of fermion matsubara frequencies - 1)
w_E = 0.5; % phonon frequency
t = 1.0; % hopping amplitude
lambda_0 = 1.5; % lambda_0 is the interaction parameter
tol_G = 1e-8; % tolerance for self-consistent calculation of G
tol_n = 1e-8; % tolerance for self-consistent calculation of n
U = 0; % Hubbard U term
%% Execute non-renormalized calculation for varying n
%{
T_list = linspace(0.01,1,21);
T_list = flip(T_list);
n_list = 0.1.*(2:10); % Occupancy per site
N_w_list = 1000*round(N_w_init./(1000*T_list));

% Preallocate zero vectors for SP susceptibilities
Chi_SP_array = zeros(length(T_list),length(n_list));

% Preallocate zero vectors for CDW susceptibilities
Chi_CDW_array = Chi_SP_array;

% Preallocate zero vector for q_max
q_max = zeros(length(T_list),length(n_list));

w = waitbar(0,'My Progress Bar'); % create a new waitbar, w with 0% progress
total_time = 0;
for jj = 1:length(n_list)
    disp(['-------------Begin n Iteration ',num2str(jj),'/',num2str(length(n_list)),'--------------'])
    n = n_list(jj);
    % Use Sigma0 for initial guess
    Sigma_trial_init = [];
    % initial 2nd guess for mu
    mu_2_guess = -0.1;
    for j = 1:length(T_list)
        tic
        disp(['-------------Begin T Iteration ',num2str(j),'/',num2str(length(T_list)),'--------------'])
        T = T_list(j);
        N_w = N_w_list(j);
        % Non-Renormalized Calculation
        calculate_Pi = []; % Exclude Phonon Self-Energy
        % Initialize without Phonon Self-Energy
        H = Holstein1D_FFT(N_k,N_w,w_E,lambda_0,T,t,n,tol_G,tol_n,Sigma_trial_init,calculate_Pi,U);
        % Converge mu, G
        [G, Sigma, D_nq, mu] = H.converge_mu(mu_2_guess);
        % Assign new initial guess for Sigma on next T iteration
        %Sigma_trial_init = Sigma;
        % Assign new initial 2nd guess for mu on next T iteration
        mu_2_guess = mu;
        % Converge Lambda
        F = abs(G).^2;
        Lambda = H.converge_Lambda(F,D_nq);
        % Calculate Chi_SP
        Chi_SP_array(j,jj) = H.Chi_SP(F,Lambda);
        % Calculate Chi_CDW
        [Chi_out, q] = H.Chi_CDW_maxq();
        Chi_CDW_array(j,jj) = mean(Chi_out);
        q_max(j,jj) = mean(abs(q));
        
        disp(['-------------End T Iteration ',num2str(j),'/',num2str(length(T_list)),'--------------'])
        total_time = total_time + toc;
        str = ['n iteration: ',num2str(jj),'/',num2str(length(n_list)), ' Complete'];
        str = [str newline 'T Iteration: ',num2str(j),'/',num2str(length(T_list)), ' Complete'];
        str = [str newline 'Time Elapsed = ',num2str(total_time)];
        w = waitbar((j + (jj-1)*(length(T_list)))/(length(n_list)*length(T_list)),w,str); % update the wait bar each iteration
    end
    disp(['-------------End n Iteration ',num2str(jj),'/',num2str(length(n_list)),'--------------'])
    
end
disp(['----------------Total Time Elapsed = ',num2str(total_time),' s------------------'])
close(w)
Chi_CDW_array(Chi_CDW_array < 0) = NaN;
%% Plotting Results
%{
% Plot SP Susceptibility 2D Mesh Plot
fig1 = figure(1);
mesh(n_list,T_list,Chi_SP_array);
xlabel '$n$'
ylabel '$T$'
zlabel '$\chi^{SP}$'
title '1D Lattice, SP Susceptibility, Non-Renormalized'
zlim([0 1]);
[lambda0_num, lambda0_den] = rat(lambda_0);
dim = [.15 .025 .3 .3];
str = [num2str(N_k),' Sites' newline num2str(N_w+1),' Matsubara Frequencies'];
str = [str newline '$\lambda_0$ = ',num2str(lambda0_num),'/',num2str(lambda0_den)...
    ', $\omega_E$ = ',num2str(w_E), ', $U$ = ', num2str(U)];
annotation('textbox',dim,'String',str,'FitBoxToText','on',...
                'interpreter','latex','FontSize',14,'linewidth',2,'backgroundcolor','w');
set(gca,'linewidth',2,'FontSize',18);
savefig(fig1, ['Chi_SP_vs_T_n_mesh_plot_Nk=',num2str(N_k),'_Nw=',num2str(N_w+1),...
    '_lambda0=',num2str(lambda0_num),'_',num2str(lambda0_den),'_Tmin=',num2str(min(T_list)),...
    '_Tmax=',num2str(max(T_list)),'_n_min=',num2str(min(n_list)),'_n_max=',...
    num2str(max(n_list)),'_tol=',num2str(tol_G),'_wE=',num2str(w_E),'_U=',num2str(U),'.fig'])

% Plot CDW Susceptibility 2D Mesh Plot
fig2 = figure(2);
mesh(n_list,T_list,Chi_CDW_array);
xlabel '$n$'
ylabel '$T$'
zlabel '$\chi^{CDW}$'
title '1D Lattice, CDW Susceptibility, Non-Renormalized'
zlim([0 10]);
[lambda0_num, lambda0_den] = rat(lambda_0);
dim = [.15 .025 .3 .3];
str = [num2str(N_k),' Sites' newline num2str(N_w+1),' Matsubara Frequencies'];
str = [str newline '$\lambda_0$ = ',num2str(lambda0_num),'/',num2str(lambda0_den)...
    ', $\omega_E$ = ',num2str(w_E), ', $U$ = ', num2str(U)];
annotation('textbox',dim,'String',str,'FitBoxToText','on',...
                'interpreter','latex','FontSize',14,'linewidth',2,'backgroundcolor','w');
set(gca,'linewidth',2,'FontSize',18);
savefig(fig2, ['Chi_CDW_vs_T_n_mesh_plot_Nk=',num2str(N_k),'_Nw=',num2str(N_w+1),...
    '_lambda0=',num2str(lambda0_num),'_',num2str(lambda0_den),'_Tmin=',num2str(min(T_list)),...
    '_Tmax=',num2str(max(T_list)),'_n_min=',num2str(min(n_list)),'_n_max=',...
    num2str(max(n_list)),'_tol=',num2str(tol_G),'_wE=',num2str(w_E),'_U=',num2str(U),'.fig'])
%}
%
% Deal with nonsense values, Plot and Save
H.plotnsave_Chi_SP_varying_n(N_w_init,T_list,n_list,Chi_SP_array,1)
H.plotnsave_Chi_CDW_varying_n(N_w_init,T_list,n_list,Chi_CDW_array,10,q_max)
%}
%% Execute renormalized calculation for varying n
T_list = linspace(0.01,1,101);
T_list = flip(T_list);
n_list = 1; % Occupancy per site
N_w_list = 10*round(N_w_init./(10*T_list));

% Preallocate zero vectors for SP susceptibilities
Chi_SP_renorm_array = zeros(length(T_list),length(n_list));

% Preallocate zero vectors for CDW susceptibilities
Chi_CDW_renorm_array = Chi_SP_renorm_array;

% Preallocate zero vector for q_max
q_max_Pi = zeros(length(T_list),length(n_list));
mu_Pi = zeros(1,length(T_list)+1);
w = waitbar(0,'My Progress Bar'); % create a new waitbar, w with 0% progress
total_time = 0;
for jj = 1:length(n_list)
    disp(['-------------Begin n Iteration ',num2str(jj),'/',num2str(length(n_list)),'--------------'])
    n = n_list(jj);
    % Use Sigma0 for initial guess
    Sigma_trial_init = [];
    % initial 2nd guess for mu
    mu_1_guess = 0;
    mu_2_guess = -1;
    Chi_out = 0;
    for j = 1:length(T_list)
        tic
        disp(['-------------Begin T Iteration ',num2str(j),'/',num2str(length(T_list)),'--------------'])
        str = ['n iteration: ',num2str(jj),'/',num2str(length(n_list)), ' Executing'];
        str = [str newline 'T Iteration: ',num2str(j),'/',num2str(length(T_list)), ' Executing'];
        str = [str newline 'Time Elapsed = ',num2str(total_time)];
        w = waitbar((j + (jj-1)*(length(T_list)))/(length(n_list)*length(T_list)),w,str); % update the wait bar each iteration
        T = T_list(j);
        N_w = N_w_list(j);
        % Renormalized Calculation
        calculate_Pi = 1; % Include Phonon Self-Energy
        % Initialize with Phonon Self-Energy
        H = Holstein1D_FFTv2(N_k,N_w,w_E,lambda_0,T,t,n,tol_G,tol_n,Sigma_trial_init,calculate_Pi,U);
        % Converge mu, G
        [G_Pi, Sigma_Pi, D_nq_Pi, Pi, mu_Pi(j+1)] = H.converge_mu(mu_1_guess, mu_2_guess);
        if isnan(G_Pi)
            Chi_SP_renorm_array(j,jj) = NaN;
            Chi_CDW_renorm_array(j,jj) = NaN;
            continue
        end
        % Assign new initial guess for Sigma on next T iteration
        Sigma_trial_init = Sigma_Pi;
        % Assign new initial guesses for mu on next T iteration
        mu_1_guess = mu_Pi(j);
        mu_2_guess = mu_Pi(j+1);
        % Converge Lambda
        F_Pi = abs(G_Pi).^2;
        Lambda_Pi = H.converge_Lambda(F_Pi,D_nq_Pi);
        % Calculate Chi_SP
        Chi_SP_renorm_array(j,jj) = H.Chi_SP(F_Pi,Lambda_Pi);
        % Calculate Chi_CDW
        [Chi_out, q] = H.Chi_CDW_maxq(Pi);
        Chi_out = mean(Chi_out);
        Chi_CDW_renorm_array(j,jj) =  Chi_out;
        
        disp(['-------------End T Iteration ',num2str(j),'/',num2str(length(T_list)),'--------------'])
        total_time = total_time + toc;

    end
    disp(['-------------End n Iteration ',num2str(jj),'/',num2str(length(n_list)),'--------------'])
    H.plot_Chi_CDW_vs_Chi_SP_fixed_n(Chi_SP_renorm_array(:,jj),Chi_CDW_renorm_array(:,jj),T_list,N_w_init)
    hold off
end
disp(['----------------Total Time Elapsed = ',num2str(total_time),' s------------------'])
close(w)
%% Deal with nonsense values, Plot and Save
%% Plotting Results
% Plot Susceptibilities
%{
fig3 = figure(3);
mesh(n_list,T_list,Chi_SP_renorm_array);
xlabel '$n$'
ylabel '$T$'
zlabel '$\chi^{SP}$'
title '1D Lattice, SP Susceptibility, Renormalized'
zlim([0 1]);
[lambda0_num, lambda0_den] = rat(lambda_0);
dim = [.15 .025 .3 .3];
str = [num2str(N_k),' Sites' newline num2str(N_w+1),' Matsubara Frequencies'];
str = [str newline '$\lambda_0$ = ',num2str(lambda0_num),'/',num2str(lambda0_den)...
    ', $\omega_E$ = ',num2str(w_E), ', $U$ = ', num2str(U)];
annotation('textbox',dim,'String',str,'FitBoxToText','on',...
                'interpreter','latex','FontSize',14,'linewidth',2,'backgroundcolor','w');
set(gca,'linewidth',2,'FontSize',18);
savefig(fig3, ['Renorm_Chi_SP_vs_T_n_mesh_plot_Nk=',num2str(N_k),'_Nw=',num2str(N_w+1),...
    '_lambda0=',num2str(lambda0_num),'_',num2str(lambda0_den),'_Tmin=',num2str(min(T_list)),...
    '_Tmax=',num2str(max(T_list)),'_n_min=',num2str(min(n_list)),'_n_max=',...
    num2str(max(n_list)),'_tol=',num2str(tol_G),'_wE=',num2str(w_E),'_U=',num2str(U),'.fig'])

% Plot CDW Susceptibility 2D Mesh Plot
fig4 = figure(4);
mesh(n_list,T_list,Chi_CDW_renorm_array);
xlabel '$n$'
ylabel '$T$'
zlabel '$\chi^{CDW}$'
title '1D Lattice, CDW Susceptibility, Renormalized'
zlim([0 10]);
[lambda0_num, lambda0_den] = rat(lambda_0);
dim = [.15 .025 .3 .3];
str = [num2str(N_k),' Sites' newline num2str(N_w+1),' Matsubara Frequencies'];
str = [str newline '$\lambda_0$ = ',num2str(lambda0_num),'/',num2str(lambda0_den)...
    ', $\omega_E$ = ',num2str(w_E), ', $U$ = ', num2str(U)];
annotation('textbox',dim,'String',str,'FitBoxToText','on',...
                'interpreter','latex','FontSize',14,'linewidth',2,'backgroundcolor','w');
set(gca,'linewidth',2,'FontSize',18);
savefig(fig4, ['Renorm_Chi_CDW_vs_T_n_mesh_plot_Nk=',num2str(N_k),'_Nw=',num2str(N_w+1),...
    '_lambda0=',num2str(lambda0_num),'_',num2str(lambda0_den),'_Tmin=',num2str(min(T_list)),...
    '_Tmax=',num2str(max(T_list)),'_n_min=',num2str(min(n_list)),'_n_max=',...
    num2str(max(n_list)),'_tol=',num2str(tol_G),'_wE=',num2str(w_E),'_U=',num2str(U),'.fig'])
%}
%
H.plotnsave_Chi_SP_varying_n(N_w_init,T_list,n_list,Chi_SP_renorm_array,1)
H.plotnsave_Chi_CDW_varying_n(N_w_init,T_list,n_list,Chi_CDW_renorm_array,10,q_max_Pi)
