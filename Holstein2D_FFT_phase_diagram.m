%% Holstein2D_FFT
%% Set LaTeX as default text interpreter
clear all
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
maxNumCompThreads(16);
%% Set System Parameters
H = Holstein2D_FFT();
W = 8;
w_cutoff = 10*W;
N_w_init = 10*ceil(w_cutoff/(10*pi)); % (minimum number of fermion matsubara frequencies - 1)
N_kx = 16; % number of sites
N_ky = N_kx;
H.init_wavevectors(N_kx,N_ky);
H.w_E = 0.5; % phonon frequency
H.t = 1.0; % hopping amplitude
H.lambda_0 = 2.5; % lambda_0 is the interaction parameter
H.coupling_constant = sqrt(H.lambda_0*H.w_E/2); %gw_E
H.tol_G = 1e-8; % tolerance for self-consistent calculation of G
H.tol_n = 1e-8; % tolerance for self-consistent calculation of n
H.tol_Lambda = 1e-3; % tolerance for self-consistent vertex calculation of Lambda
H.U = 0; % Hubbard U term
H.renorm_calc = 1; % Execute phonon self energy renormalized calculation
%% Execute renormalized calculation for varying n to first divergence of SP
n_list = [0.02 0.05:0.05:0.7 0.72:0.02:0.88 0.9 0.95 1.0];
q_max_array = cell(1,length(n_list));
T_SP = zeros(1,length(n_list));
T_CDW = T_SP;
% Preallocate zero vector for q_max
w = waitbar(0,'My Progress Bar'); % create a new waitbar, w with 0% progress
total_time = 0;
for jj = 1:length(n_list)
    T_list = linspace(0.001,0.25,101);
    T_list = flip(T_list);
    % Preallocate zero vectors for SP susceptibilities
    Chi_SP_array = zeros(length(T_list),1);
    % Preallocate zero vectors for CDW susceptibilities
    Chi_CDW_array = Chi_SP_array;
    N_w_list = 10*ceil(N_w_init./(10*T_list));
    disp(['-------------Begin n Iteration ',num2str(jj),'/',num2str(length(n_list)),'--------------'])
    H.n = n_list(jj);
    % Use Sigma0 for initial guess
    H.Sigma_trial_init = [];
    % initial 2nd guess for mu
    mu_1_guess = 0;
    mu_2_guess = -2;
    mu = zeros(1,length(T_list)+1);
    q_x = zeros(1,length(T_list));
    q_y = q_x;
    for j = 1:length(T_list)
        tic
        T = T_list(j);
        N_w = N_w_list(j);
        H.init_T_dep_objects(N_kx,N_ky,N_w,T)
        disp(['-------------Begin T Iteration ',num2str(j),'/',num2str(length(T_list)),'--------------'])
        str = ['$n$ iteration: ',num2str(jj),'/',num2str(length(n_list)), ' Executing'];
        str = [str newline 'T = ',num2str(T),', $N_w$ = ',num2str(N_w),', Iteration: ',num2str(j),'/',num2str(length(T_list)), ' Executing'];
        str = [str newline 'Time Elapsed = ',num2str(total_time)];
        w = waitbar((j + (jj-1)*(length(T_list)))/(length(n_list)*length(T_list)),w,str); % update the wait bar each iteration
        % Converge mu, G
        [G, Sigma, D_nq, mu(j+1)] = H.converge_mu(mu_1_guess, mu_2_guess);
        if isnan(G)
            T_list = T_list(T_list > T);
            q_x = q_x(T_list > T);
            q_y = q_y(T_list > T);            
            Chi_SP_array = Chi_SP_array(T_list > T);
            Chi_CDW_array = Chi_CDW_array(T_list > T);
            disp('-------------- CDW has Diverged! Terminating Loop ---------------------')
            break
        end
        % Assign new initial guess for Sigma on next T iteration
        H.Sigma_trial_init = Sigma;
        % Assign new initial guesses for mu on next T iteration
        mu_1_guess = mu(j);
        mu_2_guess = mu(j+1);
        % Converge Lambda
        F = abs(G).^2;
        Lambda = H.converge_Lambda(F,D_nq);
        % Calculate Chi_SP
        Chi_SP_out = H.Chi_SP(F,Lambda);
        Chi_SP_array(j) = Chi_SP_out;
        
        % Calculate Chi_CDW
        [Chi_CDW_out, q_x(j), q_y(j)] = H.Chi_CDW_maxq();
        Chi_CDW_array(j) = Chi_CDW_out;
        
        if Chi_SP_out == inf
            % exit inner loop if the SP susceptibility diverges
            T_list = T_list(T_list >= T);
            q_x = q_x(T_list >= T);
            q_y = q_y(T_list >= T);
            Chi_SP_array = Chi_SP_array(T_list >= T);
            Chi_CDW_array = Chi_CDW_array(T_list >= T);
            disp('-------------- SP has Diverged! Terminating Loop ---------------------')
            break
        elseif 1/Chi_CDW_out < 0 %|| (j>1 && Chi_CDW_array(j) < Chi_CDW_array(j-1))
            % exit inner loop if the CDW susceptibility diverges
            T_list = T_list(T_list > T);
            q_x = q_x(T_list > T);
            q_y = q_y(T_list > T);            
            Chi_SP_array = Chi_SP_array(T_list > T);
            Chi_CDW_array = Chi_CDW_array(T_list > T);
            disp('-------------- CDW has Diverged! Terminating Loop ---------------------')
            break
        end
        
        disp(['-------------End T Iteration ',num2str(j),'/',num2str(length(T_list)),'--------------'])
        total_time = total_time + toc;
        
    end
    disp(['-------------End n Iteration ',num2str(jj),'/',num2str(length(n_list)),'--------------'])
    q_max_array{jj} = [q_x(end); q_y(end)]; % Store max CDW wavevector at T_CDW_final
    [T_SP(jj), T_CDW(jj)] = H.linear_linear_fit_plot_Tc(w_cutoff,N_w_init,T_list,Chi_SP_array,Chi_CDW_array,q_max_array{jj});
    pause(1.0);

end
disp(['----------------Total Time Elapsed = ',num2str(total_time),' s------------------'])
close(w)
%% Phase Diagram
%T_SP(T_CDW > T_SP) = NaN;
%T_CDW(T_SP > T_CDW) = NaN;
H.phase_diagram(w_cutoff,N_w_init,n_list,q_max_array,T_SP,T_CDW);