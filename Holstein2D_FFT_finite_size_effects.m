%% Holstein2D_FFT
%% Set LaTeX as default text interpreter
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
maxNumCompThreads(16);
H = Holstein2D_FFT();
%% Set System Parameters
N_kx_list = 2.^(4:6); % number of sites
N_ky_list = N_kx_list;
W = 8;
w_cutoff = 100*W;
N_w_init = 10*ceil(w_cutoff/(10*pi)); % (minimum number of fermion matsubara frequencies - 1)
H.w_E = 1.0; % phonon frequency
H.t = 1.0; % hopping amplitude
H.lambda_0 = 2.4; % lambda_0 is the interaction parameter
H.coupling_constant = sqrt(H.lambda_0*H.w_E/2); %gw_E
H.tol_G = 1e-8; % tolerance for self-consistent calculation of G
H.tol_n = 1e-8; % tolerance for self-consistent calculation of n
H.tol_Lambda = 1e-3; % tolerance for self-consistent vertex calculation of Lambda
H.U = 0; % Hubbard U term
%% Execute renormalized calculation for varying n to first divergence of SP
H.n = 1;
CDW_cell = cell(1,length(N_kx_list));
SP_cell = CDW_cell;
T_cell = SP_cell;
H.renorm_calc = 1; % Execute phonon self energy renormalized calculation
% Preallocate zero vector for q_max
w = waitbar(0,'My Progress Bar'); % create a new waitbar, w with 0% progress
total_time = 0;
for jj = 1:length(N_kx_list)
    N_kx = N_kx_list(jj);
    N_ky = N_kx;
    H.init_wavevectors(N_kx,N_ky);
    T_list = linspace(1e-2,0.4,21);
    T_list = flip(T_list);
    % Preallocate zero vectors for SP susceptibilities
    Chi_SP_array = zeros(length(T_list),1);
    % Preallocate zero vectors for CDW susceptibilities
    Chi_CDW_array = Chi_SP_array;
    N_w_list = 10*round(N_w_init./(10*T_list));
    q_x = zeros(length(T_list),1);
    q_y = q_x;
    disp(['-------------Begin N_k Iteration ',num2str(jj),'/',num2str(length(N_kx_list)),'--------------'])
    % Use Sigma0 for initial guess
    H.Sigma_trial_init = [];
    % initial 2nd guess for mu
    mu_1_guess = 0;
    mu_2_guess = -1;
    mu_Pi = zeros(1,length(T_list)+1);
    for j = 1:length(T_list)
        tic
        T = T_list(j);
        N_w = N_w_list(j);
        disp(['-------------Begin T Iteration ',num2str(j),'/',num2str(length(T_list)),'--------------'])
        str = ['$N_{kx}$ iteration: ',num2str(jj),'/',num2str(length(N_kx_list)), ' Executing'];
        str = [str newline 'T = ',num2str(T),', $N_w$ = ',num2str(N_w),', Iteration: ',num2str(j),'/',num2str(length(T_list)), ' Executing'];
        str = [str newline 'Time Elapsed = ',num2str(total_time)];
        w = waitbar((j + (jj-1)*(length(T_list)))/(length(N_kx_list)*length(T_list)),w,str); % update the wait bar each iteration
        H.init_T_dep_objects(N_kx,N_ky,N_w,T);
        % Converge mu, G
        [G_Pi, Sigma_Pi, D_nq_Pi, mu_Pi(j+1)] = H.converge_mu(mu_1_guess, mu_2_guess);
        if isnan(G_Pi)
            T_list = T_list(T_list > T);
            Chi_SP_array = Chi_SP_array(T_list > T);
            Chi_CDW_array = Chi_CDW_array(T_list > T);
            disp('-------------- CDW has Diverged! Terminating Loop ---------------------')
            break
        end
        % Assign new initial guess for Sigma on next T iteration
        H.Sigma_trial_init = Sigma_Pi;
        % Assign new initial guesses for mu on next T iteration
        mu_1_guess = mu_Pi(j);
        mu_2_guess = mu_Pi(j+1);
        % Converge Lambda
        F_Pi = abs(G_Pi).^2;
        Lambda_Pi = H.converge_Lambda(F_Pi,D_nq_Pi);
        % Calculate Chi_SP
        Chi_SP_out = H.Chi_SP(F_Pi,Lambda_Pi);
        Chi_SP_array(j) = Chi_SP_out;
        
        % Calculate Chi_CDW
        [Chi_CDW_out, ~, ~] = H.Chi_CDW_maxq();
        Chi_CDW_array(j) = Chi_CDW_out;
        
        if Chi_SP_out == inf
            % exit inner loop if the SP susceptibility diverges
            T_list = T_list(T_list >= T);
            Chi_SP_array = Chi_SP_array(T_list >= T);
            Chi_CDW_array = Chi_CDW_array(T_list >= T);
            disp('-------------- SP has Diverged! Terminating Loop ---------------------')
            break
        elseif 1/Chi_CDW_out < 1e-14 || (j>1 && Chi_CDW_array(j) < Chi_CDW_array(j-1))
            % exit inner loop if the CDW susceptibility diverges
            T_list = T_list(T_list > T);          
            Chi_SP_array = Chi_SP_array(T_list > T);
            Chi_CDW_array = Chi_CDW_array(T_list > T);
            disp('-------------- CDW has Diverged! Terminating Loop ---------------------')
            break
        end
        
        disp(['-------------End T Iteration ',num2str(j),'/',num2str(length(T_list)),'--------------'])
        total_time = total_time + toc;
        
    end
    SP_cell{jj} = Chi_SP_array;
    CDW_cell{jj} = Chi_CDW_array;
    T_cell{jj} = T_list;
    disp(['-------------End N_kx Iteration ',num2str(jj),'/',num2str(length(N_kx_list)),'--------------'])
end
disp(['----------------Total Time Elapsed = ',num2str(total_time),' s------------------'])
close(w)
%% Plot CDW Susceptibilities
fig = figure;
for j = 1:length(N_kx_list)
    plot(T_cell{j},1./CDW_cell{j},'-o','linewidth',2,'DisplayName',['N = $',num2str(N_kx_list(j)),'^2$']);
    hold on
end
xlabel '$T$'
ylabel '1/$\chi^{CDW}(\pi,\pi)$'
xlim([0 max(T_list)])
ylim([0 0.1])
grid on
legend show
dim = [.15 .025 .3 .3];
str = ['2D Lattice, ', num2str(N_w_init),', $\omega_c = ',num2str(w_cutoff/W),'W$'];
str = [str newline '$\lambda_0$ = ',num2str(H.lambda_0),', $\omega_E$ = ',num2str(H.w_E),', $U$ = ',num2str(H.U), ', $n$ = 1'];
annotation('textbox',dim,'String',str,'FitBoxToText','on',...
    'interpreter','latex','FontSize',18,'linewidth',2,'backgroundcolor','w');
set(gca,'linewidth',2,'FontSize',24);
savefig(fig, ['Finite_Size_Effects/Finite_Size_Renorm_Chi_CDW_nonrecip_vs_T_n=1','_minNkx=',num2str(min(N_kx_list)),'_maxNk=',num2str(max(N_kx_list)),...
    '_Nw=',num2str(N_w_init), '_lambda0=',num2str(H.lambda_0),'_tol=',num2str(H.tol_G),'_wE=',num2str(H.w_E),'_U=',num2str(H.U),'.fig'])
%% Plot SP Susceptibilites
fig = figure;
for j = 1:length(N_kx_list)
    plot(T_cell{j},1./SP_cell{j},'-o','linewidth',2,'DisplayName',['N = $',num2str(N_kx_list(j)),'^2$']);
    hold on
end
xlabel '$T$'
ylabel '$1/\chi^{SP}$'
xlim([0 max(T_list)])
ylim([0 5])
grid on
legend show
dim = [.15 .025 .3 .3];
str = ['2D Lattice, ', num2str(N_w_init),', $\omega_c = ',num2str(w_cutoff/W),'W$'];
str = [str newline '$\lambda_0$ = ',num2str(H.lambda_0),', $\omega_E$ = ',num2str(H.w_E),', $U$ = ',num2str(H.U), ', $n$ = 1'];
annotation('textbox',dim,'String',str,'FitBoxToText','on',...
    'interpreter','latex','FontSize',18,'linewidth',2,'backgroundcolor','w');
set(gca,'linewidth',2,'FontSize',24);
savefig(fig, ['Finite_Size_Effects/Finite_Size_Renorm_Chi_SP_vs_T_n=1','_minNkx=',num2str(min(N_kx_list)),'_maxNk=',num2str(max(N_kx_list)),...
    '_Nw=',num2str(N_w_init), '_lambda0=',num2str(H.lambda_0),'_tol=',num2str(H.tol_G),'_wE=',num2str(H.w_E),'_U=',num2str(H.U),'.fig'])

