classdef Holstein1D_FFTv2 < handle
    % Contains methods pertinent to solving the 1D Holstein Model
    properties
        k
        w_m
        v_n
        tau
        t
        w_E
        E_kxmesh
        lambda_0
        Beta
        n
        coupling_constant
        tol
        G0_hf
        D0_non_int
        tol_G
        tol_n
        Sigma_trial_init
        calculate_Pi
        U
    end
    
    methods
        function self = Holstein1D_FFTv2(N_k,N_w,w_E,lambda_0,T,t,n_,tol_G_,tol_n_,Sigma_trial_init_,calculate_Pi,U)
            self.t = t;
            % hopping amplitude
            
            self.k = (2*pi/N_k).*linspace(-N_k/2,N_k/2,N_k+1);
            self.k = self.k(2:end);
            % k is the electron wavevector in units of a
            % k in (-N_k/2, N_k/2], where N_k is the number of sites
            
            m = linspace(-N_w/2,N_w/2,N_w+1);
            % m in [-N_w, N_w] (N_w + 1 values of m)
            
            self.w_m = self.t*pi*T.*(2.*m-1);
            % w_m is the fermionic matsubara frequency
            
            self.v_n = self.t*pi*T.*(2.*m);
            % v_n is the bosonic matsubare frequency
            
            self.w_E = w_E*self.t;
            % w_E is the phonon frequency
            
            self.tau = (0:N_w)./(T*(N_w+1));
            % tau is the imaginary time vector
            
            self.lambda_0 = lambda_0;
            % lambda is the interaction parameter
            
            self.Beta = 1/T*t;
            % inverse temperature
            
            self.n = n_;
            % n is the occupation number
            
            self.coupling_constant = sqrt(self.lambda_0*self.w_E/2); % gw_E
            
            self.tol_G = tol_G_;
            % tolerance for self-consistent calculations of G, Lambda
            
            self.tol_n = tol_n_;
            % tolerance for self-consistent calculations of n
            
            self.Sigma_trial_init = Sigma_trial_init_;
            % Initial guess for Sigma
            
            [k_xmesh,w_m_ymesh] = meshgrid(self.k,self.w_m);
            self.E_kxmesh = -2*self.t.*cos(k_xmesh);
            self.G0_hf = 1./(1i*w_m_ymesh - self.E_kxmesh);
            % Construct G0 for mu = 0 (half-filling)
            
            self.D0_non_int = repmat(-2*self.w_E.*(1./(self.w_E^2 + (self.v_n').^2)), [1 N_k]);
            % Non-interacting phonon propagator
            
            self.calculate_Pi = calculate_Pi;
            % Property that determines whether Pi is to be included in
            % self-consistent calculations
            % leave as empty variable if Pi is to be excluded
            
            self.U = U;
            % Hubbard U term for on-site repulsion
            
            self.tau = repmat((0:1:N_w).'.*self.Beta/(N_w + 1), [1 N_k]);
            % Matsubara Fourier transform of G0
        end
        
        function [T_SP, T_CDW] = find_phase_Tc(self,T_init_1,T_init_2,N_w_init,tol)
            counter = 0;
            mu_2_guess = -0.1;
            
            % First Guess
            T_1 = T_init_1;
            T_1_SP = T_1;
            T_1_CDW = T_1;
            N_w = 100*round(N_w_init/(100*T_1));
            % init T dep objects
            self.init_T_dep_objects(T_1,N_w)
            % Phonon Self-Energy Renormalized Calculation
            % Converge mu, G
            [G_Pi, Sigma_Pi, D_nq_Pi, mu_Pi] = self.converge_mu(mu_2_guess);
            mu_2_guess = mu_Pi;
            % Assign new initial guess for Sigma on next T iteration
            self.Sigma_trial_init_Pi = Sigma_Pi;
            % Converge Lambda
            F_Pi = abs(G_Pi).^2;
            Lambda_Pi = H.converge_Lambda(F_Pi,D_nq_Pi);
            % Calculate (Renorm) Chi_SP
            Chi_SP_1 = H.Chi_SP(F_Pi,Lambda_Pi);
            % Calculate (Renorm) Chi_CDW
            Chi_CDW_1 = H.Chi_CDW_qmax();
            
            counter = counter + 1;
            disp([num2str(counter),' Iterations Completed to Find Tc'])
            disp(['|1/Chi_SP|=',num2str(abs(1/Chi_SP_1))])
            disp(['|1/Chi_CDW|=',num2str(abs(1/Chi_CDW_1))])
            
            % Second Guess
            T_2 = T_init_2;
            T_2_SP = T_2;
            T_2_CDW = T_2;
            N_w = 100*round(N_w_init/(100*T_2));
            % init T dep objects
            self.init_T_dep_objects(T_2,N_w)
            % Phonon Self-Energy Renormalized Calculation
            % Converge mu, G
            [G_Pi, Sigma_Pi, D_nq_Pi, ~] = self.converge_mu(mu_2_guess);
            mu_2_guess_CDW = mu_Pi;
            mu_2_guess_SP = mu_Pi;
            % Assign new initial guess for Sigma on next T iteration
            self.Sigma_trial_init_Pi = Sigma_Pi;
            % Converge Lambda
            F_Pi = abs(G_Pi).^2;
            Lambda_Pi = H.converge_Lambda(F_Pi,D_nq_Pi);
            % Calculate (Renorm) Chi_SP
            Chi_SP_2 = 1./H.Chi_SP(F_Pi,Lambda_Pi);
            % Calculate (Renorm) Chi_CDW
            Chi_CDW_2 = 1./H.Chi_CDW_qmax();
            
            counter = counter + 1;
            disp([num2str(counter),' Iterations Completed to Find Tc'])
            disp(['|1/Chi_SP|=',num2str(abs(1/Chi_SP_2))])
            disp(['|1/Chi_CDW|=',num2str(abs(1/Chi_CDW_2))])
            
            while (abs(1/Chi_SP_2) > tol || abs(1/Chi_CDW_2) > tol)
                % Calculate improved guess T for CDW
                T_3_CDW = T_2_CDW - (1/Chi_CDW_2)*(T_2_CDW - T_1_CDW)...
                    /(1/Chi_CDW_2 - 1/Chi_CDW_1);
                
                % Calculate improved guess T for SP
                T_3_SP = T_2_SP - (1/Chi_SP_2)*(T_2_SP - T_1_SP)...
                    /(1/Chi_SP_2 - 1/Chi_SP_1);
                
                % Calculate improved guess for 1/Chi_SP
                N_w = 100*round(N_w_init/(100*T_3_SP));
                % init T dep objects
                self.init_T_dep_objects(T_3_SP,N_w)
                % Phonon Self-Energy Renormalized Calculation
                % Converge mu, G
                [G_Pi, ~, D_nq_Pi, ~] = self.converge_mu(mu_2_guess_SP);
                %mu_2_guess_SP = mu_Pi;
                % Assign new initial guess for Sigma on next T iteration
                self.Sigma_trial_init_Pi = [];
                % Converge Lambda
                F_Pi = abs(G_Pi).^2;
                Lambda_Pi = H.converge_Lambda(F_Pi,D_nq_Pi);
                % Calculate (Renorm) Chi_SP
                Chi_SP_3 = H.Chi_SP(F_Pi,Lambda_Pi);
                
                % Recursively Assign Variables
                T_1_SP = T_2_SP;
                T_2_SP = T_3_SP;
                Chi_SP_1 = Chi_SP_2;
                Chi_SP_2 = Chi_SP_3;
                
                % Calculate improved guess for 1/Chi_SP
                N_w = 100*round(N_w_init/(100*T_3_CDW));
                % init T dep objects
                self.init_T_dep_objects(T_3_CDW,N_w)
                % Phonon Self-Energy Renormalized Calculation
                % Converge mu, G
                [~, ~, ~, ~] = self.converge_mu(mu_2_guess_CDW);
                %mu_2_guess_CDW = mu_Pi;
                % Assign new initial guess for Sigma on next T iteration
                self.Sigma_trial_init_Pi = [];
                % Calculate (Renorm) Chi_CDW
                Chi_CDW_3 = H.Chi_CDW_qmax();
                
                % Recursively Assign Variables
                T_1_CDW = T_2_CDW;
                T_2_CDW = T_3_CDW;
                Chi_CDW_1 = Chi_CDW_2;
                Chi_CDW_2 = Chi_CDW_3;
                
                counter = counter + 1;
                disp([num2str(counter),' Iterations Completed to Find Tc'])
                disp(['|1/Chi_SP|=',num2str(abs(1/Chi_SP_2))])
                disp(['|1/Chi_CDW|=',num2str(abs(1/Chi_CDW_2))])
            end
            
            if abs(1/Chi_SP_2) < tol
                T_SP = T_3_SP;
                T_CDW = NaN;
                
            elseif abs(1/Chi_CDW_2) < tol
                T_SP = NaN;
                T_CDW = T_3_CDW;
            else
                T_SP = NaN;
                T_CDW = NaN;
            end
            
        end
        
        function init_T_dep_objects(self,T,N_w)
            m = linspace(-N_w/2,N_w/2,N_w+1);
            % m in [-N_w, N_w] (N_w + 1 values of m)
            
            self.w_m = self.t*pi*T.*(2.*m-1);
            % w_m is the fermionic matsubara frequency
            
            self.v_n = self.t*pi*T.*(2.*m);
            % v_n is the bosonic matsubare frequency
            self.Beta = 1/(T*self.t);
            % inverse temperature
            
            [k_xmesh,w_m_ymesh] = meshgrid(self.k,self.w_m);
            self.E_kxmesh = -2*self.t.*cos(k_xmesh);
            self.G0_hf = 1./(1i*w_m_ymesh - self.E_kxmesh);
            % Construct G0 for mu = 0 (half-filling)
            
            self.D0 = repmat(-2*self.w_E.*(1./(self.w_E^2 + (self.v_n').^2)), [1 N_k]);
            % Non-interacting phonon propagator
            
        end
        
        function out = acquire_figure_data(~,fig)
            dataObjs = findobj(fig,'-property','YData');
            out = dataObjs.YData;
        end
        
        function plot_Chi_CDW_vs_Chi_SP_fixed_n(self,Chi_SP,Chi_CDW,T_list,N_w_init)
            N_k = length(self.k);
            tol_ = self.tol_G;
            lambda0_ = self.lambda_0;
            w_E_ = self.w_E;
            U_ = self.U;
            n_ = self.n;
            
            fig = figure;
            set(fig,'color','w','Units', 'Normalized', 'OuterPosition', [0 0 1/2 3/4])
            plot(T_list,Chi_SP,'-','linewidth',2,'DisplayName','SP')
            hold on
            plot(T_list,Chi_CDW,'-','linewidth',2,'DisplayName','CDW')
            dim = [.15 .025 .3 .3];
            str = [num2str(N_k),' Sites' newline num2str(N_w_init),' Matsubara Frequencies Min'];
            str = [str newline '$\lambda_0$ = ',num2str(lambda0_),', $n$ = ',num2str(n_)...
                ', $\omega_E$ = ',num2str(w_E_),', $U$ = ',num2str(U_)];
            
            annotation('textbox',dim,'String',str,'FitBoxToText','on',...
                'interpreter','latex','FontSize',14,'linewidth',2,'backgroundcolor','w');
            legend show
            xlabel '$T$'
            ylim([0, 10])
            xlim([0, 1])
            grid on
            set(gca,'linewidth',2,'FontSize',18)
            
            if isempty(self.calculate_Pi)
                title(str)
                savefig(fig, ['SP_vs_CDW_vs_T_n=',num2str(n_),'_Nk=',num2str(N_k),'_Nw=',num2str(N_w_init),...
                    '_lambda0=',num2str(lambda0_),'_Tmin=',num2str(min(T_list)),...
                    '_Tmax=',num2str(max(T_list)),'_tol=',num2str(tol_),'_U=',num2str(U_),'.fig'])
                
            else
                title([str,' Renormalized'])
                savefig(fig, ['Renorm_SP_vs_CDW_vs_T_n=',num2str(n_),'_Nk=',num2str(N_k),'_Nw=',num2str(N_w_init),...
                    '_lambda0=',num2str(lambda0_),'_Tmin=',num2str(min(T_list)),...
                    '_Tmax=',num2str(max(T_list)),'_tol=',num2str(tol_),'_U=',num2str(U_),'.fig'])
                
            end
        end
        
        function plotnsave_Chi_SP_varying_n(self,N_w_init,T_list,n_list,Chi_array,y_lim)
            [lambda0_num, lambda0_den] = rat(self.lambda_0);
            N_k = length(self.k);
            tol_ = self.tol_G;
            w_E_ = self.w_E;
            U_ = self.U;
            
            fig = figure;
            set(fig,'color','w','Units', 'Normalized', 'OuterPosition', [0 0 1/2 3/4])
            for j = 1:length(n_list)
                label = ['$n = $',num2str(n_list(j))];
                plot(T_list,Chi_array(:,j),'-','linewidth',2,'DisplayName',label)
                hold on
            end
            dim = [.15 .025 .3 .3];
            str = [num2str(N_k),' Sites' newline num2str(N_w_init),' Matsubara Frequencies Min'];
            str = [str newline '$\lambda_0$ = ',num2str(lambda0_num),'/',num2str(lambda0_den)...
                ', $\omega_E$ = ',num2str(w_E_),', $U$ = ',num2str(U_)];
            
            annotation('textbox',dim,'String',str,'FitBoxToText','on',...
                'interpreter','latex','FontSize',14,'linewidth',2,'backgroundcolor','w');
            legend show
            xlabel '$T$'
            ylim([0, y_lim])
            xlim([0,1])
            str = '1D Lattice, SP Suscepitibility for Varying n';
            grid on
            set(gca,'linewidth',2,'FontSize',18)
            
            if isempty(self.calculate_Pi)
                title(str)
                savefig(fig, ['Chi_SP_vs_T_varying_n_min_n=',num2str(min(n_list)),...
                    '_max_n=',num2str(max(n_list)),'_Nk=',num2str(N_k),'_Nw=',num2str(N_w_init),...
                    '_lambda0=',num2str(lambda0_num),'_',num2str(lambda0_den),'_Tmin=',num2str(min(T_list)),...
                    '_Tmax=',num2str(max(T_list)),'_tol=',num2str(tol_),'_U=',num2str(U_),'.fig'])
                
            else
                title([str,' Renormalized'])
                savefig(fig, ['Renorm_Chi_SP_vs_T_varying_n_min_n=',num2str(min(n_list)),...
                    '_max_n=',num2str(max(n_list)),'_Nk=',num2str(N_k),'_Nw=',num2str(N_w_init),...
                    '_lambda0=',num2str(lambda0_num),'_',num2str(lambda0_den),'_Tmin=',num2str(min(T_list)),...
                    '_Tmax=',num2str(max(T_list)),'_tol=',num2str(tol_),'_U=',num2str(U_),'.fig'])
                
            end
            
        end
        
        function plotnsave_Chi_CDW_varying_n(self,N_w_init,T_list,n_list,Chi_array,y_lim,q_max)
            [lambda0_num, lambda0_den] = rat(self.lambda_0);
            N_k = length(self.k);
            tol_ = self.tol_G;
            w_E_ = self.w_E;
            U_ = self.U;
            
            fig = figure;
            set(fig,'color','w','Units', 'Normalized', 'OuterPosition', [0 0 1/2 3/4])
            for j = 1:length(n_list)
                label = ['$n = $',num2str(n_list(j))];
                plot(T_list,Chi_array(:,j),'-','linewidth',2,'DisplayName',label)
                hold on
            end
            dim = [.15 .025 .3 .3];
            str = [num2str(N_k),' Sites' newline num2str(N_w_init),' Matsubara Frequencies Min'];
            str = [str newline '$\lambda_0$ = ',num2str(lambda0_num),'/',num2str(lambda0_den)...
                ', $\omega_E$ = ',num2str(w_E_),', $U$ = ',num2str(U_)];
            
            annotation('textbox',dim,'String',str,'FitBoxToText','on',...
                'interpreter','latex','FontSize',14,'linewidth',2,'backgroundcolor','w');
            legend show
            xlabel '$T$'
            ylim([0,y_lim])
            xlim([0,1])
            str = '1D Lattice, CDW Suscepitibility for Varying n';
            grid on
            set(gca,'linewidth',2,'FontSize',18)
            
            if isempty(self.calculate_Pi)
                title(str)
                savefig(fig, ['Chi_CDW_vs_T_varying_n_min_n=',num2str(min(n_list)),...
                    '_max_n=',num2str(max(n_list)),'_Nk=',num2str(N_k),'_Nw=',num2str(N_w_init),...
                    '_lambda0=',num2str(lambda0_num),'_',num2str(lambda0_den),'_Tmin=',num2str(min(T_list)),...
                    '_Tmax=',num2str(max(T_list)),'_tol=',num2str(tol_),'_U=',num2str(U_),'.fig'])
                
            else
                title([str,' Renormalized'])
                savefig(fig, ['Renorm_Chi_CDW_vs_T_varying_n_min_n=',num2str(min(n_list)),...
                    '_max_n=',num2str(max(n_list)),'_Nk=',num2str(N_k),'_Nw=',num2str(N_w_init),...
                    '_lambda0=',num2str(lambda0_num),'_',num2str(lambda0_den),'_Tmin=',num2str(min(T_list)),...
                    '_Tmax=',num2str(max(T_list)),'_tol=',num2str(tol_),'_U=',num2str(U_),'.fig'])
                
                
            end
            
        end
        
        
        function plotnsave_Chi_CDW(self,N_w_init,T_list,Chi_CDW_renorm,Chi_CDW,Chi_CDW_0,y_lim)
            [lambda0_num, lambda0_den] = rat(self.lambda_0);
            N_w = length(self.w_m);
            N_k = length(self.k);
            tol_ = self.tol_G;
            n_ = self.n;
            w_E_ = self.w_E;
            
            fig = figure;
            set(fig,'color','w','Units', 'Normalized', 'OuterPosition', [0 0 1/2 3/4])
            plot(T_list,Chi_CDW_renorm,'-o','linewidth',2,'DisplayName','$\chi^{CDW}$, Renormalized')
            hold on
            plot(T_list,Chi_CDW,'-s','linewidth',2,'DisplayName','$\chi^{CDW}$')
            plot(T_list,Chi_CDW_0,'-d','linewidth',2,'DisplayName','$\chi_0^{CDW}$')
            dim = [.15 .025 .3 .3];
            str = [num2str(N_k),' Sites' newline num2str(N_w_init + 1),' Matsubara Frequencies Min'];
            str = [str newline '$\lambda_0$ = ',num2str(lambda0_num),'/',num2str(lambda0_den)...
                ', $\omega_E$ = ',num2str(w_E_), ', $n$ = ', num2str(n_)];
            
            annotation('textbox',dim,'String',str,'FitBoxToText','on',...
                'interpreter','latex','FontSize',14,'linewidth',2,'backgroundcolor','w');
            legend show
            xlabel '$T$'
            ylim([0, y_lim])
            xlim([0,1])
            title '1D Lattice, CDW Suscepitibility'
            grid on
            set(gca,'linewidth',2,'FontSize',18)
            
            savefig(fig, ['Chi_CDW_vs_T_n=',num2str(self.n),'_Nk=',num2str(N_k),'_minNw=',num2str(N_w_init + 1),...
                '_lambda0=',num2str(lambda0_num),'_',num2str(lambda0_den),'_Tmin=',num2str(min(T_list)),...
                '_Tmax=',num2str(max(T_list)),'_tol=',num2str(tol_),'.fig'])
            
            ylim('auto')
            set(gca,'XScale','log','YScale','log')
            savefig(fig, ['Log_Scale_Chi_CDW_vs_T_n=',num2str(n_),'_Nk=',num2str(N_k),'_minNw=',num2str(N_w_init + 1),...
                '_lambda0=',num2str(lambda0_num),'_',num2str(lambda0_den),'_Tmin=',num2str(min(T_list)),...
                '_Tmax=',num2str(max(T_list)),'_tol=',num2str(tol_),'.fig'])
            
            
            Chi_CDW_Table = table(T_list',Chi_CDW_0',Chi_CDW',Chi_CDW_renorm','VariableNames',{'T','Chi_CDW_0','Chi_CDW','Chi_CDW_Renorm'});
            save(['Chi_CDW_vs_T__n=',num2str(n_),'_Nk=',num2str(N_k),'_minNw=',num2str(N_w_init + 1),...
                '_lambda0=',num2str(lambda0_num),'_',num2str(lambda0_den),'_Tmin=',num2str(min(T_list)),...
                '_Tmax=',num2str(max(T_list)),'_tol=',num2str(tol_),'.mat'],'Chi_CDW_Table')
        end
        
        function plotnsave_Chi_SP(self,N_w_init,T_list,Chi_SP_renorm,Chi_SP,Chi_SP_0,y_lim)
            [lambda0_num, lambda0_den] = rat(self.lambda_0);
            N_w = length(self.w_m);
            N_k = length(self.k);
            tol_ = self.tol_G;
            n_ = self.n;
            w_E_ = self.w_E;
            
            fig = figure;
            set(fig,'color','w','Units', 'Normalized', 'OuterPosition', [0 0 1/2 3/4])
            plot(T_list,Chi_SP_renorm,'-o','linewidth',2,'DisplayName','$\chi^{SP}$, Renormalized')
            hold on
            plot(T_list,Chi_SP,'-s','linewidth',2,'DisplayName','$\chi^{SP}$')
            plot(T_list,Chi_SP_0,'-d','linewidth',2,'DisplayName','$\chi_0^{SP}$')
            dim = [.15 .025 .3 .3];
            str = [num2str(N_k),' Sites' newline num2str(N_w_init + 1),' Matsubara Frequencies Min'];
            str = [str newline '$\lambda_0$ = ',num2str(lambda0_num),'/',num2str(lambda0_den)...
                ', $\omega_E$ = ',num2str(w_E_), ', $n$ = ', num2str(n_)];
            %{
            if ~isempty(T_c_index)
                str = [str newline '$T_c$ = ',num2str(T_c)];
            end
            
                %}
                annotation('textbox',dim,'String',str,'FitBoxToText','on',...
                    'interpreter','latex','FontSize',14,'linewidth',2,'backgroundcolor','w');
                legend show
                xlabel '$T$'
                %ylabel '$\chi^{SP}$'
                ylim([0, y_lim])
                xlim([0,1])
                title '1D Lattice, SP Susceptibility'
                grid on
                set(gca,'linewidth',2,'FontSize',18)
                
                savefig(fig, ['Chi_SP_vs_T_n=',num2str(n_),'_Nk=',num2str(N_k),'_minNw=',num2str(N_w_init + 1),...
                    '_lambda0=',num2str(lambda0_num),'_',num2str(lambda0_den),'_Tmin=',num2str(min(T_list)),...
                    '_Tmax=',num2str(max(T_list)),'_tol=',num2str(tol_),'.fig'])
                
                ylim('auto')
                set(gca,'XScale','log','YScale','log')
                savefig(fig, ['Log_Scale_Chi_SP_vs_T_n=',num2str(n_),'_Nk=',num2str(N_k),'_minNw=',num2str(N_w_init + 1),...
                    '_lambda0=',num2str(lambda0_num),'_',num2str(lambda0_den),'_Tmin=',num2str(min(T_list)),...
                    '_Tmax=',num2str(max(T_list)),'_tol=',num2str(tol_),'.fig'])
                
                
                Chi_SP_Table = table(T_list',Chi_SP_0',Chi_SP',Chi_SP_renorm','VariableNames',{'T','Chi_SP_0','Chi_SP','Chi_SP_Renorm'});
                save(['Chi_SP_vs_T__n=',num2str(n_),'_Nk=',num2str(N_k),'_minNw=',num2str(N_w_init + 1),...
                    '_lambda0=',num2str(lambda0_num),'_',num2str(lambda0_den),'_Tmin=',num2str(min(T_list)),...
                    '_Tmax=',num2str(max(T_list)),'_tol=',num2str(tol_),'.mat'],'Chi_SP_Table')
        end
        
        function out = Chi_CDW_0(self,mu,q)
            k_ = self.k;
            N_k = length(k_);
            E_k = -2*self.t.*cos(k_);
            C_k = ones(N_k,1);
            beta = self.Beta;
            f_E_k = 1./(exp(beta.*(E_k - mu)) + 1);
            k_index = repmat(1:N_k,[1,3]);
            q_index = round(q/(2*pi/N_k));
            kq_index = k_index((N_k + 1 + q_index):(2*N_k + q_index));
            Pi0_numer = f_E_k(1,:) - f_E_k(1,kq_index);
            Pi0_denom = E_k - E_k(1,kq_index);
            Pi0_summand = Pi0_numer./Pi0_denom;
            lhopital = f_E_k(1,abs(Pi0_denom) < 1e-14);
            Pi0_summand(1,abs(Pi0_denom) < 1e-14) = beta.*lhopital.*(-1 + lhopital);
            out = -real((2/N_k)*Pi0_summand*C_k);
        end
        
        function [Chi_out, q_max] = Chi_CDW_maxq(self,Pi_nq)
            q = self.k;
            Pi_vn0 = Pi_nq(ceil(end/2),:);
            Chi_CDW_bar = -1/(self.coupling_constant^2).*real(Pi_vn0);
            Chi_CDW_ = Chi_CDW_bar./(1-(self.lambda_0 - self.U).*Chi_CDW_bar);
            q_max_cond = (abs(Pi_vn0) == max(abs(Pi_vn0)));
            Chi_out = Chi_CDW_(q_max_cond);
            q_max = q(q_max_cond);
        end
        
        function out = Chi_CDW(self,q,Pi_nq)
            N_k = length(self.k);
            Pi_vn0 = Pi_nq(ceil(end/2),:);
            q_index = round(q/(2*pi/N_k) + N_k/2);
            Chi_CDW_bar = -1/(self.coupling_constant^2).*real(Pi_vn0(q_index));
            out = Chi_CDW_bar./(1-(self.lambda_0 - self.U).*Chi_CDW_bar);
        end
        
        function out = Chi_SP_0(self,mu)
            beta = self.Beta;
            k_ = self.k;
            N_k = length(self.k);
            E_k = -2.*self.t.*cos(k_);
            summand_numer = 1 - 2./(exp(beta.*(E_k - mu)) + 1);
            summand_denom = (2*(E_k-mu));
            zero_index = find(abs(summand_denom) < 1e-14);
            sum_vector = summand_numer./summand_denom;
            
            if ~isempty(zero_index)
                sum_vector(zero_index) = beta/4;
            end
            
            out = (1/N_k).*sum(sum_vector);
            
        end
        
        function out = Chi_SP(self,F,Lambda)
            beta = self.Beta;
            N_k = length(self.k);
            N_w = length(self.w_m);
            C_mk = ones(N_k*N_w,1);
            out = real(1/(N_k*beta)*reshape(F.*Lambda,[],1).'*C_mk);
        end
        
        function out = converge_Lambda(self,F,D_nq)
            disp('--------------------- Begin Lambda Convergence ----------------------')
            % Initiate objects for anderson mixing
            M = 3;
            res_cell = cell(1,M);
            trial_cell = cell(1,M);
            
            N_w = length(self.w_m);
            N_k = length(self.k);
            beta = self.Beta;
            U_ = self.U;
            Lambda_factor = self.coupling_constant^2/(beta*N_k);
            Lambda_trial = ones(size(F));
            Lambda_calc = [];
            Lambda_tol = self.tol_G;
            F_Lambda_trial = F.*Lambda_trial;
            C_mk = ones(N_k*N_w,1);
            trial_diff_1 = 1;
            counter = 0;
            divergence_counter = 0;
            divergence_lim = 100;
            iter_lim = 100000;
            r_phase_factor = repmat(exp(-1i*pi*(2/N_k - 1).*(0:N_k-1)),[N_w 1]);
            j_phase_factor = repmat(exp(-1i*pi*(N_w-1)/N_w.*(0:N_w-1).'),[1 N_k]);
            while abs(trial_diff_1) > Lambda_tol && divergence_counter < divergence_lim && counter < iter_lim
                counter = counter + 1;
                if counter <= M
                    Lambda_calc = 1 - (Lambda_factor).*(ifft2(r_phase_factor.*conj(j_phase_factor)...
                        .*fft2(D_nq).*fft2(F_Lambda_trial)))...
                        - (U_/(N_k*beta))*reshape(F_Lambda_trial,[],1).'*C_mk;
                    res_2 = Lambda_calc - Lambda_trial;
                    trial_diff_2 = max(res_2,[],'all');
                    if abs(trial_diff_2) > abs(trial_diff_1)
                        divergence_counter = divergence_counter + 1;
                    else
                        divergence_counter = 0;
                    end
                    res_cell{counter} = res_2;
                    trial_cell{counter} = Lambda_trial;
                    trial_diff_1 = trial_diff_2;
                    Lambda_trial = Lambda_calc;
                    F_Lambda_trial = F.*Lambda_trial;
                elseif counter == M + 1
                    Lambda_calc = 1 - (Lambda_factor).*(ifft2(r_phase_factor.*conj(j_phase_factor)...
                        .*fft2(D_nq).*fft2(F_Lambda_trial)))...
                        - (U_/(N_k*beta))*reshape(F_Lambda_trial,[],1).'*C_mk;
                    res_2 = (Lambda_calc - Lambda_trial);
                    trial_diff_2 = max(res_2,[],'all');
                    if abs(trial_diff_2) > abs(trial_diff_1)
                        divergence_counter = divergence_counter + 1;
                    else
                        divergence_counter = 0;
                    end
                    trial_diff_1 = trial_diff_2;
                    res_1 = res_2;
                    b = 1;
                    Lambda_trial = self.anderson_mixing(b,res_cell,res_1,trial_cell,Lambda_calc);
                    F_Lambda_trial = F.*Lambda_trial;
                else
                    Lambda_calc = 1 - (Lambda_factor).*(ifft2(r_phase_factor.*conj(j_phase_factor)...
                        .*fft2(D_nq).*fft2(F_Lambda_trial)))...
                        - (U_/(N_k*beta))*reshape(F_Lambda_trial,[],1).'*C_mk;
                    trial_diff_2 = max(abs(Lambda_calc - Lambda_trial),[],'all');
                    
                    res_cell = res_cell(2:end);
                    res_cell{M} = (res_1);
                    trial_cell = trial_cell(2:end);
                    trial_cell{M} = (Lambda_trial);
                    
                    if abs(trial_diff_2) < abs(trial_diff_1)
                        if counter <= 10
                            Lambda_trial = self.anderson_mixing(b,res_cell,res_1,trial_cell,Lambda_calc);
                        else
                            Lambda_trial = Lambda_calc;
                        end
                        divergence_counter = 0;
                    elseif abs(trial_diff_2) > abs(trial_diff_1)
                        Lambda_trial = Lambda_calc;
                        divergence_counter = divergence_counter + 1;
                    end
                    
                    trial_diff_1 = trial_diff_2;
                    res_1 = res_2;
                    F_Lambda_trial = F.*Lambda_trial;
                end
                disp([num2str(counter),' Iterations Completed to Converge Lambda'])
                disp(['|trial_diff|=',num2str(abs(trial_diff_1))])
            end
            
            if divergence_counter >= divergence_lim || counter >= iter_lim
                Lambda_calc(:) = 1e10;
            end
            
            out = Lambda_calc;
            disp('--------------------- End Lambda Convergence ----------------------')
            
        end
        
        function [G_out, Sigma_out, D_out, Pi_out, mu_out] = converge_mu(self,mu_1_guess,mu_2_guess)
            
            disp('--------------------- Begin mu, G Convergence ----------------------')
            beta = self.Beta;
            tol_n_ = self.tol_n;
            n_ = self.n;
            N_k = length(self.k);
            N_w = length(self.w_m);
            C_k = ones(N_k,1);
            C_mk = ones(N_k*N_w,1);
            k_ = self.k;
            E_k = -2*self.t.*cos(k_);
            G0_hf_ = self.G0_hf;
            
            % Check for half-filling
            if n_ == 1
                mu_out = 0;
                [G_out, Sigma_out, D_out, Pi_out] = self.converge_G(G0_hf_,mu_out);
                return
            end
            
            %First Guess
            mu_1 = mu_1_guess;
            n0_1 = (2/N_k)*(1./(exp(beta.*(E_k - mu_1)) + 1))*C_k;
            G0_1 = 1./(1./G0_hf_ + mu_1);
            [G_1, Sigma_mk, D_nq, Pi_nq] = self.converge_G(G0_1,mu_1);
            if isnan(G_1)
                G_out = NaN;
                Sigma_out = NaN;
                D_out = NaN;
                Pi_out = NaN;
                mu_out = NaN;
                return
            end
            n_1 = n0_1 + (2/(beta*N_k))*reshape(G_1-G0_1,[],1).'*C_mk;
            %n_1 = 2/(beta*N_k).*reshape(G_1,[],1).'*C_mk;
            delta_n_1 = real(n_ - n_1);
            
            counter = 1;
            disp([num2str(counter),' Iterations Completed to Converge mu'])
            disp(['|delta_n|=',num2str(abs(delta_n_1))])
            
            if abs(delta_n_1) < tol_n_
                G_out = G_1;
                Sigma_out = Sigma_mk;
                D_out = D_nq;
                Pi_out = Pi_nq;
                mu_out = mu_1;
                return
            end
            
            % Second Guess
            mu_2 = mu_2_guess;
            n0_2 = (2/N_k)*(1./(exp(beta.*(E_k - mu_2)) + 1))*C_k;
            G0_2 = 1./(1./G0_hf_ + mu_2);
            [G_2, Sigma_mk, D_nq, Pi_nq] = self.converge_G(G0_2,mu_2);
            if isnan(G_2)
                G_out = NaN;
                Sigma_out = NaN;
                Pi_out = NaN;
                D_out = NaN;
                mu_out = NaN;
                return
            end
            n_2 = n0_2 + (2/(beta*N_k))*reshape(G_2-G0_2,[],1).'*C_mk;
            %n_2 = 2/(beta*N_k).*reshape(G_2,[],1).'*C_mk;
            delta_n_2 = real(n_ - n_2);
            
            counter = counter + 1;
            disp([num2str(counter),' Iterations Completed to Converge mu'])
            disp(['|delta_n|=',num2str(abs(delta_n_2))])
            
            if abs(delta_n_2) < tol_n_
                G_out = G_2;
                Sigma_out = Sigma_mk;
                D_out = D_nq;
                Pi_out = Pi_nq;
                mu_out = mu_2;
                return
            end
            
            while abs(delta_n_2) > tol_n_
                if abs(delta_n_2 - delta_n_1) < 1e-14
                    delta_n_2 = delta_n_1 + 1;
                end
                % Apply Discrete Newton-Raphson to improve (3+j)rd guess for mu
                mu_3 = mu_2 - (mu_2-mu_1)*(delta_n_2)/(delta_n_2-delta_n_1);
                
                % (3+j)rd Guess
                n0_3 = (2/N_k)*(1./(exp(beta.*(E_k - mu_3)) + 1))*C_k;
                G0_3 = 1./(1./G0_hf_ + mu_3);
                [G_3, Sigma_mk, D_nq, Pi_nq] = self.converge_G(G0_3,mu_3);
                if isnan(G_3)
                    G_out = NaN;
                    Sigma_out = NaN;
                    D_out = NaN;
                    Pi_out = NaN;
                    mu_out = NaN;
                    return
                end
                n_3 = n0_3 + (2/(beta*N_k))*reshape(G_3-G0_3,[],1).'*C_mk;
                %n_3 = 2/(beta*N_k)*reshape(G_3,[],1).'*C_mk;
                delta_n_3 = real(n_ - n_3);
                
                % Recursively assign values to previous guesses
                mu_1 = mu_2;
                mu_2 = mu_3;
                delta_n_1 = delta_n_2;
                delta_n_2 = delta_n_3;
                
                counter = counter + 1;
                disp([num2str(counter),' Iterations Completed to Converge mu'])
                disp(['|delta_n|=',num2str(abs(delta_n_2))])
            end
            
            G_out = G_3;
            Sigma_out = Sigma_mk;
            D_out = D_nq;
            Pi_out = Pi_nq;
            mu_out = mu_2;
            disp('--------------------- End mu, G Convergence ----------------------')
            
        end
        
        
        
        function [G_out, Sigma_out, D_out, Pi_out] = converge_G(self,G0,mu)
            disp('--------------------- Begin G Convergence ----------------------')
            % Utilizes Anderson mixing
            % Condition for CDW_bar used to identify CDW divergence
            
            counter = 1;
            Pi_counter = 1;
            tol_ = self.tol_G;
            M = 3; % number of previous iterations to be mixed
            res_cell = cell(1,M);
            trial_cell = cell(1,M);
            
            % Construct Fourier transforms of G0
            beta = self.Beta;
            tau_ = self.tau;
            E_k = self.E_kxmesh;
            G0_tau = -exp((beta-tau_).*(E_k-mu))./(exp(beta.*(E_k-mu)) + 1);
            G0_tau_conj = exp(-(beta-tau_).*(E_k-mu))./(exp(-beta.*(E_k-mu)) + 1);
            
            % Construct Sigma0 as guess for first T iteration)
            Sigma_trial_1 = self.Sigma_trial_init;
            if isempty(Sigma_trial_1)
                Sigma_trial_1 = self.Sigma0(mu);
            elseif ~isempty(Sigma_trial_1)
                % pad previous temperature iteration with Sigma0 if
                % number of matsubara frequencies are different
                N_w_prev = length(Sigma_trial_1(:,1));
                N_w = length(self.w_m);
                Sigma0_matrix = self.Sigma0(mu);
                Sigma_trial_1 = [Sigma0_matrix(1:round((N_w - N_w_prev)/2),:); Sigma_trial_1;...
                    Sigma0_matrix((round((N_w + N_w_prev)/2 + 1)):N_w,:)];
            end
            
            D0 = self.D0_non_int;
            near_unity = 1 - 10^(round(log10(1/beta)));
            
            % First Guess
            Sigma_trial = Sigma_trial_1;
            trial_cell{1} = Sigma_trial;
            G_trial = 1./(1./G0 - Sigma_trial);
            Pi_nq = self.Pi(G_trial,G0,G0_tau,G0_tau_conj);
            dressed_CDW_bar = -((self.lambda_0 - self.U)...
                /self.coupling_constant^2).*real(Pi_nq(ceil(end/2),:));
            if sum(dressed_CDW_bar > 1) ~= 0
                Pi_nq(dressed_CDW_bar > 1) = near_unity./(-((self.lambda_0 - self.U)...
                /self.coupling_constant^2));
            end
            if ~isempty(self.calculate_Pi)
                D_nq = 1./(1./D0 - Pi_nq);
            else
                D_nq = D0;
            end
            Sigma_calc = Sigma(self,G_trial,G0,G0_tau,D_nq);
            res_1 = Sigma_calc - Sigma_trial;
            res_cell{1} = res_1;
            trial_diff_1 = max(res_1,[],'all');
            disp([num2str(counter),' Iterations Completed to Converge G'])
            disp(['|trial_diff|=',num2str(abs(trial_diff_1))])
            
            % Initial damping factor
            b = 0.1;
            
            iter_lim_init = 500;
            iter_lim = iter_lim_init;
            % (2 + j)th Guess
            while abs(trial_diff_1) > tol_ && counter <= iter_lim
                counter = counter + 1;
                if counter <= M
                    Sigma_trial = Sigma_trial + b.*res_1;
                    trial_cell{counter} = Sigma_trial;
                    G_trial = 1./(1./G0 - Sigma_trial);
                    Pi_nq = self.Pi(G_trial,G0,G0_tau,G0_tau_conj);
                    dressed_CDW_bar = -((self.lambda_0 - self.U)...
                        /self.coupling_constant^2).*real(Pi_nq(ceil(end/2),:));
                    if sum(dressed_CDW_bar > 1) ~= 0 && Pi_counter <= 10
                        Pi_nq(ceil(end/2), dressed_CDW_bar > 1) = near_unity./(-((self.lambda_0 - self.U)...
                /self.coupling_constant^2));
                    end
                    if ~isempty(self.calculate_Pi)
                        D_nq = 1./(1./D0 - Pi_nq);
                    else
                        D_nq = D0;
                    end
                    Sigma_calc = Sigma(self,G_trial,G0,G0_tau,D_nq);
                    res_2 = Sigma_calc - Sigma_trial;
                    res_cell{counter} = res_2;
                    trial_diff_2 = max(res_2,[],'all');
                    b = b/(1 - trial_diff_2/trial_diff_1);
                    if abs(trial_diff_2) < abs(trial_diff_1)
                        Pi_counter = Pi_counter + 1;
                    else
                        Pi_counter = 0;
                    end
                    trial_diff_1 = trial_diff_2;
                    res_1 = res_2;
                    
                elseif counter == M + 1
                    Sigma_trial = Sigma_trial + b.*res_1;
                    G_trial = 1./(1./G0 - Sigma_trial);
                    Pi_nq = self.Pi(G_trial,G0,G0_tau,G0_tau_conj);
                    dressed_CDW_bar = -((self.lambda_0 - self.U)...
                        /self.coupling_constant^2).*real(Pi_nq(ceil(end/2),:));
                    if sum(dressed_CDW_bar > 1) ~= 0 && Pi_counter <= 10
                        Pi_nq(ceil(end/2), dressed_CDW_bar > 1) = near_unity./(-((self.lambda_0 - self.U)...
                /self.coupling_constant^2));
                    end
                    if ~isempty(self.calculate_Pi)
                        D_nq = 1./(1./D0 - Pi_nq);
                    else
                        D_nq = D0;
                    end
                    Sigma_calc = Sigma(self,G_trial,G0,G0_tau,D_nq);
                    res_2 = Sigma_calc - Sigma_trial;
                    trial_diff_2 = max(res_2,[],'all');
                    if abs(trial_diff_2) < abs(trial_diff_1)
                        Pi_counter = Pi_counter + 1;
                    else
                        Pi_counter = 0;
                    end
                    b = b/(1 - trial_diff_2/trial_diff_1);
                    trial_diff_1 = trial_diff_2;
                    res_1 = res_2;
                    
                else
                    % Initialize anderson mixing
                    b = 1;
                    Sigma_trial_1 = ...
                        self.anderson_mixing(b,res_cell,res_1,trial_cell,Sigma_trial);
                    trial_cell = trial_cell(2:end);
                    trial_cell{M} = Sigma_trial;
                    res_cell = res_cell(2:end);
                    res_cell{M} = res_1;
                    G_trial = 1./(1./G0 - Sigma_trial);
                    Pi_nq = self.Pi(G_trial,G0,G0_tau,G0_tau_conj);
                    dressed_CDW_bar = -((self.lambda_0 - self.U)...
                        /self.coupling_constant^2).*real(Pi_nq(ceil(end/2),:));
                    if sum(dressed_CDW_bar > 1) ~= 0 && Pi_counter <= 10
                        Pi_nq(ceil(end/2), dressed_CDW_bar > 1) = near_unity./(-((self.lambda_0 - self.U)...
                /self.coupling_constant^2));
                        Pi_counter = Pi_counter + 1;
                    elseif sum(dressed_CDW_bar > 1) ~= 0 && Pi_counter > 10
                        disp('---------------- CDW is Divergent! ---------------')
                        disp('---------------- Outputting NaNs -----------------')
                        pause(1.0)
                        G_out = NaN;
                        Sigma_out = NaN;
                        D_out = NaN;
                        Pi_out = NaN;
                        return
                    end
                    if ~isempty(self.calculate_Pi)
                        D_nq = 1./(1./D0 - Pi_nq);
                    else
                        D_nq = D0;
                    end
                    Sigma_calc = Sigma(self,G_trial,G0,G0_tau,D_nq);
                    res_2 = Sigma_calc - Sigma_trial_1;
                    trial_diff_2 = max(res_2,[],'all');
                    if abs(trial_diff_2) < abs(trial_diff_1)
                        Pi_counter = Pi_counter + 1;
                    else
                        Pi_counter = 0;
                    end
                    if abs(trial_diff_1) < 1e-2
                        iter_lim = iter_lim + 1;
                    else
                        iter_lim = iter_lim_init;
                    end
                    trial_diff_1 = trial_diff_2;
                    res_1 = res_2;
                    Sigma_trial = Sigma_trial_1;
                end
                
                disp([num2str(counter),' Iterations Completed to Converge G'])
                disp(['|trial_diff|=',num2str(abs(trial_diff_1))])
                
            end
            
            if counter > iter_lim
                disp('Convergence Failed! Electron Self-Energy Too Large!')
                G_out = NaN.*ones(size(G0));
                Sigma_out = NaN.*ones(size(G0));
                D_out = NaN.*ones(size(G0));
                Pi_out = NaN.*ones(size(G0));
                return
            end
            
            G_out = G_trial;
            Pi_out = self.Pi(G_trial,G0,G0_tau,G0_tau_conj);
            if ~isempty(self.calculate_Pi)
                D_out = 1./(1./D0 - Pi_out);
            else
                D_out = D0;
            end
            Sigma_out = self.Sigma(G_trial,G0,G0_tau,D_out);
            disp('--------------------- End G Convergence ----------------------')
            
        end
        
        function out = anderson_mixing(self, b, prev_res_cell, new_res, prev_trial_cell, new_trial)
            M = length(prev_res_cell);
            N_k = length(self.k);
            
            cell_reshape = @(x)reshape(x,[],1);
            prev_res_cell = cellfun(cell_reshape,prev_res_cell,'UniformOutput',false);
            prev_trial_cell = cellfun(cell_reshape,prev_trial_cell,'UniformOutput',false);
            
            new_res = reshape(new_res,[],1);
            new_trial = reshape(new_trial,[],1);
            
            delta_res = zeros(length(prev_res_cell{1}),M);
            delta_trial = delta_res;
            delta_res_mat = zeros(M);
            delta_res_vec = zeros(M,1);
            
            j = 1;
            while j <= M
                delta_res(:,j) = new_res - prev_res_cell{M + 1 - j};
                delta_trial(:,j) = new_trial - prev_trial_cell{M + 1 - j};
                j = j + 1;
            end
            
            j = 1;
            while j <= M
                delta_res_mat(j,j) = delta_res(:,j)'*delta_res(:,j);
                delta_res_vec(j,1) = delta_res(:,j)'*new_res;
                jj = j + 1;
                while jj <= M
                    delta_res_mat(jj,j) = delta_res(:,jj)'*delta_res(:,j);
                    delta_res_mat(j,jj) = delta_res_mat(jj,j);
                    jj = jj + 1;
                end
                j = j + 1;
            end
            
            if rcond(delta_res_mat) > 1e-20
                delta_res_coeff = delta_res_mat\delta_res_vec;
            else
                delta_res_coeff = pinv(delta_res_mat)*delta_res_vec;
            end
            
            new_res_bar = new_res;
            new_trial_bar = new_trial;
            j = 1;
            while j <= length(delta_res_coeff)
                new_res_bar = new_res_bar - delta_res_coeff(j,1).*delta_res(:,j);
                new_trial_bar = new_trial_bar - delta_res_coeff(j,1).*delta_trial(:,j);
                j = j + 1;
            end
            out = reshape(new_trial_bar + b.*new_res_bar,[],N_k);
            
        end
        
        function Sigma_mk = Sigma(self,G_mk,G0_mk,G0_tau,D_nq)
            % FUNCTION DESCRIPTION ----------------------------------------------------
            % Calculates electron self-energy using fft to evaluate propagator convolutions
            % ARGUMENTS ---------------------------------------------------------------
            % G_mk is the single-electron propagator under trial
            % G0_mk is non-interacting single-electron propagator
            % G0_tau is the analytic Fermionic Fourier transform of G0_mk
            % D_nq is the phonon self-energy renormalized phonon propagator
            % OUTPUTS -----------------------------------------------------------------
            % Sigma_mk is the calculated electron self-energy
            
            N_k = length(self.k);
            N_w = length(self.w_m);
            beta = self.Beta;
            factor_Sigma = self.coupling_constant^2/(N_k*beta);
            
            r_phase_factor = repmat(exp(-1i*pi*(2/N_k - 1).*(0:N_k-1)),[N_w 1]);
            j_phase_factor_1 = repmat(exp(-1i*pi*(N_w-1)/N_w.*(0:N_w-1).'),[1 N_k]);
            j_phase_factor_2 = repmat(exp(1i*pi.*(0:N_w-1)).', [1 N_k]);
            
            Sigma_mk = -(factor_Sigma).*ifft2(conj(j_phase_factor_1).*r_phase_factor.*fft2(D_nq)...
                .*conj(j_phase_factor_2).*fft(j_phase_factor_2.*fft(G_mk-G0_mk,N_w,1) + beta*G0_tau,N_k,2));
        end
        
        function Pi_nq = Pi(self,G_mk,G0_mk,G0_tau,G0_tau_conj)
            % FUNCTION DESCRIPTION ----------------------------------------------------
            % Calculates phonon self-energy using fft to evaluate propagator convolutions
            % ARGUMENTS ---------------------------------------------------------------
            % G_mk is the single-electron propagator under trial
            % G0_mk is non-interacting single-electron propagator
            % G0_tau is the analytic Fermionic Fourier transform of G0_mk
            % G0_tau_conj is the analytic Fermionic Fourier transform of conj(G0_mk)
            % OUTPUTS -----------------------------------------------------------------
            % Pi_nq is the calculated phonon self-energy
            
            N_k = length(self.k);
            N_w = length(self.w_m);
            beta = self.Beta;
            factor_Pi = 2*self.coupling_constant^2/(N_k*beta);
            
            r_phase_factor = repmat(exp(-1i*pi*(2/N_k - 1).*(0:N_k-1)), [N_w 1]);
            j_phase_factor_1 = repmat(exp(-1i*pi*(N_w-1)/N_w.*(0:N_w-1).'), [1 N_k]);
            j_phase_factor_2 = repmat(exp(1i*pi.*(0:N_w-1)).', [1 N_k]);
            
            Pi_nq = (factor_Pi).*ifft2(j_phase_factor_1.*r_phase_factor...
                .*fft(j_phase_factor_2.*fft(G_mk-G0_mk,N_w,1) + beta.*G0_tau,N_k,2)...
                .*fft(j_phase_factor_2.*fft(conj(G_mk-G0_mk),N_w,1) + beta.*G0_tau_conj,N_k,2));
        end
        
        
        function out = Pi0_vn0(self,mu)
            q = self.k;
            k_ = q;
            N_k = length(q);
            N_q = N_k;
            t_ = self.t;
            beta = self.Beta;
            
            E_k = -2*t_.*cos(k_xmesh);
            f_E_k = 1./(exp(beta.*(E_k - mu)) + 1);
            
            Pi0_factor = 2*self.coupling_constant^2/N_k;
            
            C_k = ones(N_k,1);
            ones_vector = ones(1,N_k);
            sum_vector = zeros(1,N_q);
            
            j = 1;
            while j <= N_q
                kq = k_ + q(j).*ones_vector;
                E_kq = -2*t_.*cos(kq);
                f_E_kq = 1./(exp(beta.*(E_kq - mu)) + 1);
                numer = (f_E_k - f_E_kq);
                denom = (E_k - E_kq);
                summand = numer./denom;
                zero_kq_index = find(abs(denom) < 1e-14);
                if isempty(zero_kq_index)
                    sum_vector(1,j) = summand*C_k;
                else
                    lhopital = 1./(exp(beta.*(E_kq(1,zero_kq_index) - mu)) + 1);
                    summand(1,zero_kq_index) = beta.*lhopital.*(-1 + lhopital);
                    sum_vector(1,j) = summand*C_k;
                end
                j = j + 1;
            end
            out = Pi0_factor.*sum_vector;
        end
        
        function out = Pi0(self,mu)
            v_n_ = self.v_n;
            N_v = length(v_n_);
            q = self.k;
            N_k = length(q);
            N_q = N_k;
            t_ = self.t;
            beta = self.Beta;
            
            [k_xmesh,v_n_ymesh] = meshgrid(self.k,v_n_);
            zero_nn_index = find(v_n_ == 0);
            E_k_xmesh = -2*t_.*cos(k_xmesh);
            f_E_k_xmesh = 1./(exp(beta.*(E_k_xmesh - mu)) + 1);
            
            Pi0_factor = 2*self.coupling_constant^2/N_k;
            
            C_k = ones(N_k,1);
            ones_matrix = ones(N_v,N_k);
            sum_matrix = zeros(N_v,N_q);
            
            j = 1;
            while j <= N_q
                kq_xmesh = k_xmesh + q(j).*ones_matrix;
                E_kq_xmesh = -2*t_.*cos(kq_xmesh);
                f_E_kq_xmesh = 1./(exp(beta.*(E_kq_xmesh - mu)) + 1);
                numer = (f_E_k_xmesh - f_E_kq_xmesh);
                denom = (1i.*(v_n_ymesh) + E_k_xmesh - E_kq_xmesh);
                summand = numer./denom;
                zero_kq_index = find(abs(denom(zero_nn_index,:)) < 1e-14);
                if isempty(zero_kq_index)
                    sum_matrix(:,j) = summand*C_k;
                else
                    lhopital = 1./(exp(beta.*(E_kq_xmesh(1,zero_kq_index) - mu)) + 1);
                    summand(zero_nn_index,zero_kq_index) = beta.*lhopital.*(-1 + lhopital);
                    sum_matrix(:,j) = summand*C_k;
                end
                j = j + 1;
            end
            out = Pi0_factor.*sum_matrix;
            
        end
        
        function out = Sigma0(self,mu)
            k_ = self.k;
            N_k = length(k_);
            beta = self.Beta;
            [k_xmesh, w_m_ymesh] = meshgrid(k_,self.w_m);
            E_k_xmesh = -2*self.t.*cos(k_xmesh);
            C_k = ones(N_k,1);
            BE_w_E = 1./(exp(beta.*self.w_E) - 1);
            FD_E_k_mu = 1./(exp(beta.*(E_k_xmesh - mu)) + 1);
            Sigma0_mkk = (1 + BE_w_E - FD_E_k_mu)./(1i.*w_m_ymesh - self.w_E - (E_k_xmesh-mu))...
                + (BE_w_E + FD_E_k_mu)./(1i.*w_m_ymesh + self.w_E - (E_k_xmesh-mu));
            Sigma0_m = Sigma0_mkk*C_k;
            out = (self.coupling_constant^2/N_k).*repmat(Sigma0_m,1,N_k);
        end
    end
end

