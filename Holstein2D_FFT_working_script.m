%% Holstein1D_FFT
%% Set LaTeX as default text interpreter
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
maxNumCompThreads(16);
H_Pi = Holstein2D_FFT();
H = Holstein2D_FFT();
%% Set System Parameters
N_kx = 64; % number of sites
N_ky = N_kx;
W = 8;
w_cutoff = 10*W;
N_w_init = 10*ceil(w_cutoff/(10*pi)); % (minimum number of fermion matsubara frequencies - 1)
w_E = 0.5; % phonon frequency
t = 1.0; % hopping amplitude
lambda_0 = 1.5; % lambda_0 is the interaction parameter
tol_G = 1e-8; % tolerance for self-consistent calculation of G
tol_n = 1e-8; % tolerance for self-consistent calculation of n
U = 0; % Hubbard U term
n = 1;
H_Pi.init_T_indep_objects(N_kx,N_ky,w_E,lambda_0,t,n,tol_G,tol_n,U);
H.init_T_indep_objects(N_kx,N_ky,w_E,lambda_0,t,n,tol_G,tol_n,U);
%%
T_list = flip(linspace(0.01,1,51));
N_w_list = 10.*round(N_w_init./(T_list*10));
Chi_SP_renorm_list = zeros(1,length(T_list));
Chi_CDW_renorm_list = zeros(1,length(T_list));

Chi_SP_list = zeros(1,length(T_list));
Chi_CDW_list = zeros(1,length(T_list));

w = waitbar(0,'My Progress Bar'); % create a new waitbar, w with 0% progress
total_time = 0;
for j = 1:length(T_list)
    tic
    disp(['-------------Begin T Iteration ',num2str(j),'/',num2str(length(T_list)),'--------------'])
    w = waitbar(j/length(T_list),w,['Iteration: ',num2str(j),'/',num2str(length(T_list)),' Executing',...
        newline 'T = ',num2str(T), newline 'Time Elapsed = ',num2str(total_time)]); % update the wait bar each iteration
    T = T_list(j);
    N_w = N_w_list(j);
    H_Pi.init_T_dep_objects(N_kx,N_ky,N_w,T,t);
    %% Phonon Self-Energy Renormalized Calculation
    H_Pi.renorm_calc = 1; % Include Phonon Self-Energy
    % Initialize with Phonon Self-Energy
    % Converge mu, G
    [G_Pi, Sigma_Pi, D_nq_Pi, mu_Pi] = H_Pi.converge_mu();
    if ~isnan(G_Pi)
        % Assign new initial guess for Sigma on next T iteration
        H_Pi.Sigma_trial_init = Sigma_Pi;
        % Converge Lambda
        F_Pi = abs(G_Pi).^2;
        Lambda_Pi = H_Pi.converge_Lambda(F_Pi,D_nq_Pi);
        % Calculate (Renorm) Chi_SP
        Chi_SP_renorm_list(j) = H_Pi.Chi_SP(F_Pi,Lambda_Pi);
        % Calculate (Renorm) Chi_CDW
        [Chi_out, qx_max(j), qy_max(j)] = H_Pi.Chi_CDW_maxq();
        Chi_CDW_renorm_list(j) = Chi_out;
    else
        Chi_SP_renorm_list(T_list <= T) = NaN;
        Chi_CDW_renorm_list(T_list <= T) = NaN;
        Chi_SP_list(T_list <= T) = NaN;
        Chi_CDW_list(T_list <= T) = NaN;
        disp('-------------- CDW has Diverged! Terminating Loop ---------------------')
        break
    end

    disp('-------------------- Renormalized Calculation Complete ----------------------')
    %
    %% Unrenormalized Calculation
    H.init_T_dep_objects(N_kx,N_ky,N_w,T,t);
    H.calculate_Pi = []; % Exclude Phonon Self-Energy
    Sigma_trial_init = [];
    % Converge mu, G
    [G, Sigma, D_nq, mu] = H.converge_mu();
    if ~isnan(G)
        % Assign new initial guess for Sigma on next T iteration
        H.Sigma_trial_init = Sigma;
        % Converge Lambda
        F = abs(G).^2;
        Lambda = H.converge_Lambda(F,D_nq);
        % Calculate Chi_SP
        Chi_SP_list(j) = H.Chi_SP(F,Lambda);
        % Calculate Chi_CDW
        [Chi_out, qx_max(j), qy_max(j)] = H.Chi_CDW_maxq();
        Chi_CDW_list(j) = Chi_out;
    else
        Chi_SP_renorm_list(T_list <= T) = NaN;
        Chi_CDW_renorm_list(T_list <= T) = NaN;
        Chi_SP_list(T_list <= T) = NaN;
        Chi_CDW_list(T_list <= T) = NaN;
        disp('-------------- CDW has Diverged! Terminating Loop ---------------------')
        break
    end
    %}
    disp(['-------------End T Iteration ',num2str(j),'/',num2str(length(T_list)),'--------------'])
    total_time = total_time + toc;
end
disp(['----------------Total Time Elapsed = ',num2str(total_time),' s------------------'])
close(w)
Chi_CDW_list(Chi_CDW_list < 0) = NaN;
%% Plotting
figure
plot(T_list,Chi_SP_renorm_list,'-b','linewidth',2,'DisplayName','Renormalized $\chi^{SP}$')
hold on
plot(T_list,Chi_SP_list,'--b','linewidth',2,'DisplayName','Unrenormalized $\chi^{SP}$')
str = [num2str(N_kx),' $\times$ ',num2str(N_ky),' Sites' newline '$\omega_{\mathrm{cutoff}}$ = $', num2str(w_cutoff),'t$'];
str = [str newline '$\lambda_0$ = ',num2str(lambda_0),', $\omega_E$ = ',num2str(w_E), ', $n$ = ', num2str(n)];
dim = [.15 .025 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on',...
    'interpreter','latex','FontSize',14,'linewidth',2,'backgroundcolor','w');
legend show
xlabel '$T$'
ylim([0,1/2])
xlim([0,1])
grid on
set(gca,'linewidth',2,'FontSize',18)
hold off

figure
plot(T_list,Chi_CDW_renorm_list,'-r','linewidth',2,'DisplayName','Renormalized $\chi^{CDW}$')
hold on
plot(T_list,Chi_CDW_list,'--r','linewidth',2,'DisplayName','Unrenormalized $\chi^{CDW}$')
str = [num2str(N_kx),' $\times$ ',num2str(N_ky),' Sites' newline '$\omega_{\mathrm{cutoff}}$ = $', num2str(w_cutoff),'t$'];
str = [str newline '$\lambda_0$ = ',num2str(lambda_0),', $\omega_E$ = ',num2str(w_E), ', $n$ = ', num2str(n)];
dim = [.15 .025 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on',...
    'interpreter','latex','FontSize',14,'linewidth',2,'backgroundcolor','w');
legend show
xlabel '$T$'
ylim([0,20])
xlim([0,1])
grid on
set(gca,'linewidth',2,'FontSize',18)
hold off
%%
%{
% Plot and deal with nonsense values
Chi_CDW_list(Chi_CDW_list <= 0) = NaN;
Chi_SP_renorm_list(Chi_CDW_renorm_list < 0) = NaN;
Chi_CDW_renorm_list(Chi_CDW_renorm_list < 0) = NaN;
y_lim = 20;
H.plotnsave_Chi_SP(T_list,Chi_SP_renorm_list,Chi_SP_list,Chi_SP_0_list,1)
H.plotnsave_Chi_CDW(T_list,Chi_CDW_renorm_list,Chi_CDW_list,Chi_CDW_0_list,10)
%}