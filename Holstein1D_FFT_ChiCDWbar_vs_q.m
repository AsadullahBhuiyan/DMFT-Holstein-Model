%% Holstein1D_FFT
%% Set LaTeX as default text interpreter
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
maxNumCompThreads(16);
%% Set System Parameters
N_k = 256; % number of sites
N_w_init = 50; % (minimum number of fermion matsubara frequencies - 1)
w_E = 0.5; % phonon frequency
t = 1.0; % hopping amplitude
lambda_0 = 1.5; % lambda_0 is the interaction parameter
tol_G = 1e-8; % tolerance for self-consistent calculation of G
tol_n = 1e-8; % tolerance for self-consistent calculation of n
U = 0; % Hubbard U term
%% Execute Non-interacting calculation for Chi_CDW_vs_q_vs_T_vs_n
%{
T_list = [0.0002 0.02 0.01];
T_list = flip(T_list);
n_list = 0.2:0.1:0.9; % Occupancy per site
N_w_list = 100*round(N_w_init./(100*T_list));
Sigma_trial_init = [];
lambda_0 = 0;
w = waitbar(0,'My Progress Bar'); % create a new waitbar, w with 0% progress
total_time = 0;
for j = 1:length(T_list)
    disp(['-------------Begin T Iteration ',num2str(j),'/',num2str(length(T_list)),'--------------'])
    T = T_list(j);
    N_w = N_w_list(j);
    fig = figure;
    % Preallocate zero vectors for CDW susceptibilities
    Chi_CDW_renorm_array = zeros(length(n_list),N_k);
    % initial guesses for mu
    mu_1_guess = 0;
    mu_2_guess = 1;
    for jj = 1:length(n_list)
        tic
        disp(['-------------Begin n Iteration ',num2str(jj),'/',num2str(length(n_list)),'--------------'])
        n = n_list(jj);
        % Use Sigma0 for initial guess
        % Non-Renormalized Calculation
        calculate_Pi = []; % Exclude Phonon Self-Energy
        % Initialize without Phonon Self-Energy
        H = Holstein1D_FFT(N_k,N_w,w_E,lambda_0,T,t,n,tol_G,tol_n,Sigma_trial_init,calculate_Pi,U);
        % Converge mu, G
        [G, ~, ~, mu] = H.converge_mu(mu_1_guess,mu_2_guess);
        % Assign new initial guess for Sigma on next T iteration
        k_ = H.k;
        for jjj = 1:length(k_)
            Chi_CDW_0_q(jj,jjj) = H.Chi_CDW_0(-mu,k_(jjj));
        end
        % Plot for fixed n
        plot(k_/pi,Chi_CDW_0_q(jj,:),'linewidth',2,'DisplayName',['$n = $ ',num2str(n)])
        hold on
        
        disp(['-------------End n Iteration ',num2str(jj),'/',num2str(length(n_list)),'--------------'])
        total_time = total_time + toc;
        str = ['T iteration: ',num2str(j),'/',num2str(length(T_list)), ' Complete'];
        str = [str newline 'n Iteration: ',num2str(jj),'/',num2str(length(n_list)), ' Complete'];
        str = [str newline 'Time Elapsed = ',num2str(total_time)];
        w = waitbar((jj + (j-1)*(length(n_list)))/(length(n_list)*length(T_list)),w,str); % update the wait bar each iteration
    end
    xlabel '$q/\pi$'
    ylabel '$\chi_0^{CDW}(q)$'
    xlim([-0.1 1])
    grid on
    legend show
    dim = [.15 .025 .3 .3];
    str = ['1D Lattice, $N$ = ',num2str(N_k) newline num2str(N_w_init),' Matsubara Frequencies Min'];
    str = [str newline '$\lambda_0$ = ',num2str(lambda_0),', $\omega_E$ = ',num2str(w_E),', $U$ = ',num2str(U)...
        , ', T = ',num2str(T)];
    str = [str newline 'Non-Interacting'];
    annotation('textbox',dim,'String',str,'FitBoxToText','on',...
        'interpreter','latex','FontSize',14,'linewidth',2,'backgroundcolor','w');
    set(gca,'linewidth',2,'FontSize',18);
    savefig(fig, ['Non_int_Chi_CDW_vs_q_n_plot_Nk=',num2str(N_k),'_minNw=',num2str(N_w_init+1),...
        '_lambda0=',num2str(lambda_0),'_T=',num2str(T),'_n_min=',num2str(min(n_list)),'_n_max=',...
        num2str(max(n_list)),'_tol=',num2str(tol_G),'_wE=',num2str(w_E),'_U=',num2str(U),'.fig'])
    hold off
    disp(['-------------End T Iteration ',num2str(j),'/',num2str(length(T_list)),'--------------'])
end
disp(['----------------Total Time Elapsed = ',num2str(total_time),' s------------------'])
close(w)
%}
%% Execute non-renormalized calculation for Chi_CDW_vs_q_vs_T_vs_n
%{
T_list = [0.001 0.005 0.01 0.02 0.05 0.1];
T_list = flip(T_list);
n_list = 0.1.*(2:10); % Occupancy per site
N_w_list = N_w_init.*ones(1,length(T_list));

w = waitbar(0,'My Progress Bar'); % create a new waitbar, w with 0% progress
total_time = 0;
for j = 1:length(T_list)
    disp(['-------------Begin T Iteration ',num2str(j),'/',num2str(length(T_list)),'--------------'])
    T = T_list(j);
    N_w = N_w_list(j);
    fig = figure;
    % Preallocate zero vectors for CDW susceptibilities
    Chi_CDW_array = zeros(length(n_list),N_k);
    mu_2_guess = -0.1;
    for jj = 1:length(n_list)
        tic
        disp(['-------------Begin n Iteration ',num2str(jj),'/',num2str(length(n_list)),'--------------'])
        n = n_list(jj);
        % Use Sigma0 for initial guess
        Sigma_trial_init = [];
        % Non-Renormalized Calculation
        calculate_Pi = []; % Exclude Phonon Self-Energy
        % Initialize without Phonon Self-Energy
        H = Holstein1D_FFT(N_k,N_w,w_E,lambda_0,T,t,n,tol_G,tol_n,Sigma_trial_init,calculate_Pi,U);
        % Converge mu, G
        [G, Sigma, ~, mu] = H.converge_mu(mu_2_guess);
        % Assign new initial guess for Sigma on next T iteration
        Sigma_trial_init = Sigma;
        % Assign new initial 2nd guess for mu on next T iteration
        mu_2_guess = mu;
        % Calculate Pi
        Pi_nq = H.Pi(G);
        % Calculate Chi_CDW_bar at Pi(v_n = 0,:)
        Chi_CDW_bar = (-1/H.coupling_constant^2).*real(Pi_nq(ceil(end/2),:));
        Chi_CDW_array(jj,:) = Chi_CDW_bar./(1 - (lambda_0 - U).*Chi_CDW_bar);
        % Plot for fixed n
        plot(H.k(end/2:end)./pi,Chi_CDW_array(jj,end/2:end),'linewidth',2,'DisplayName',['$n = $ ',num2str(n)])
        hold on
        
        disp(['-------------End n Iteration ',num2str(jj),'/',num2str(length(n_list)),'--------------'])
        total_time = total_time + toc;
        str = ['T iteration: ',num2str(j),'/',num2str(length(T_list)), ' Complete'];
        str = [str newline 'n Iteration: ',num2str(jj),'/',num2str(length(n_list)), ' Complete'];
        str = [str newline 'Time Elapsed = ',num2str(total_time)];
        w = waitbar((jj + (j-1)*(length(n_list)))/(length(n_list)*length(T_list)),w,str); % update the wait bar each iteration
    end
    xlabel '$q/\pi$'
    ylabel '$\chi^{CDW}(q)$'
    xlim([0 1])
    grid on
    legend show
    dim = [.15 .025 .3 .3];
    str = ['1D Lattice, ', num2str(N_w_init),' Matsubara Frequencies Min'];
    str = [str newline '$\lambda_0$ = ',num2str(lambda_0),', $\omega_E$ = ',num2str(w_E),', $U$ = ',num2str(U),...
         ', T = ',num2str(T)];
    annotation('textbox',dim,'String',str,'FitBoxToText','on',...
        'interpreter','latex','FontSize',14,'linewidth',2,'backgroundcolor','w');
    set(gca,'linewidth',2,'FontSize',18);
    savefig(fig, ['Chi_CDW_vs_q_n_plot_Nk=',num2str(N_k),'_minNw=',num2str(N_w_init+1),...
        '_lambda0=',num2str(lambda_0),'_T=',num2str(T),'_n_min=',num2str(min(n_list)),'_n_max=',...
        num2str(max(n_list)),'_tol=',num2str(tol_G),'_wE=',num2str(w_E),'_U=',num2str(U),'.fig'])
    hold off
    disp(['-------------End T Iteration ',num2str(j),'/',num2str(length(T_list)),'--------------'])
end
disp(['----------------Total Time Elapsed = ',num2str(total_time),' s------------------'])
close(w)
%}
%% Execute renormalized calculation for Chi_CDW_vs_q_vs_T_vs_n
T_list = [0.02 0.05 0.1];
T_list = flip(T_list);
n_list = 0.2:0.1:1.0; % Occupancy per site
N_w_list = 100*round(N_w_init./(100*T_list));
mu_Pi = zeros(1,length(T_list)+1);
Sigma_trial_init = [];
w = waitbar(0,'My Progress Bar'); % create a new waitbar, w with 0% progress
total_time = 0;
for j = 1:length(T_list)
    disp(['-------------Begin T Iteration ',num2str(j),'/',num2str(length(T_list)),'--------------'])
    T = T_list(j);
    N_w = N_w_list(j);
    fig = figure;
    % Preallocate zero vectors for CDW susceptibilities
    Chi_CDW_renorm_array = zeros(length(n_list),N_k);
    % initial guesses for mu
    mu_1_guess = 0;
    mu_2_guess = -1;
    for jj = 1:length(n_list)
        tic
        disp(['-------------Begin n Iteration ',num2str(jj),'/',num2str(length(n_list)),'--------------'])
        n = n_list(jj);
        % Use Sigma0 for initial guess
        % Non-Renormalized Calculation
        calculate_Pi = 1; % Include Phonon Self-Energy
        % Initialize with Phonon Self-Energy
        H = Holstein1D_FFT(N_k,N_w,w_E,lambda_0,T,t,n,tol_G,tol_n,Sigma_trial_init,calculate_Pi,U);
        % Converge mu, G
        [G_Pi, Sigma_Pi, ~, mu_Pi(j+1)] = H.converge_mu(mu_1_guess,mu_2_guess);
        if isnan(G_Pi)
            continue
        end
        % Assign new initial guess for Sigma on next T iteration
        %Sigma_trial_init = Sigma_Pi;
        % Assign new initial guesses for mu on next T iteration
        mu_1_guess = mu_Pi(j);
        mu_2_guess = mu_Pi(j+1);
        % Calculate Pi
        Pi_nq = H.Pi(G_Pi);
        % Calculate Chi_CDW_bar at Pi(v_n = 0,:)
        Chi_CDW_bar_renorm = (-1/H.coupling_constant^2).*real(Pi_nq(ceil(end/2),:));
        Chi_CDW_renorm_array(jj,:) = Chi_CDW_bar_renorm./(1 - (lambda_0 - U).*Chi_CDW_bar_renorm);
        % Plot for fixed n
        plot(H.k./pi,Chi_CDW_renorm_array(jj,:),'linewidth',2,'DisplayName',['$n = $ ',num2str(n)])
        hold on
        
        disp(['-------------End n Iteration ',num2str(jj),'/',num2str(length(n_list)),'--------------'])
        total_time = total_time + toc;
        str = ['T iteration: ',num2str(j),'/',num2str(length(T_list)), ' Complete'];
        str = [str newline 'n Iteration: ',num2str(jj),'/',num2str(length(n_list)), ' Complete'];
        str = [str newline 'Time Elapsed = ',num2str(total_time)];
        w = waitbar((jj + (j-1)*(length(n_list)))/(length(n_list)*length(T_list)),w,str); % update the wait bar each iteration
    end
    xlabel '$q/\pi$'
    ylabel '$\chi^{CDW}(q)$'
    xlim([-0.1 1])
    grid on
    legend show
    dim = [.15 .025 .3 .3];
    str = ['1D Lattice, $N$ = ',num2str(N_k) newline num2str(N_w_init),' Matsubara Frequencies Min'];
    str = [str newline '$\lambda_0$ = ',num2str(lambda_0),', $\omega_E$ = ',num2str(w_E),', $U$ = ',num2str(U)...
        , ', T = ',num2str(T)];
    str = [str newline 'Renormalized'];
    annotation('textbox',dim,'String',str,'FitBoxToText','on',...
        'interpreter','latex','FontSize',14,'linewidth',2,'backgroundcolor','w');
    set(gca,'linewidth',2,'FontSize',18);
    savefig(fig, ['Renorm_Chi_CDW_vs_q_n_plot_Nk=',num2str(N_k),'_minNw=',num2str(N_w_init+1),...
        '_lambda0=',num2str(lambda_0),'_T=',num2str(T),'_n_min=',num2str(min(n_list)),'_n_max=',...
        num2str(max(n_list)),'_tol=',num2str(tol_G),'_wE=',num2str(w_E),'_U=',num2str(U),'.fig'])
    hold off
    disp(['-------------End T Iteration ',num2str(j),'/',num2str(length(T_list)),'--------------'])
end
disp(['----------------Total Time Elapsed = ',num2str(total_time),' s------------------'])
close(w)