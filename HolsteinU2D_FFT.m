classdef HolsteinU2D_FFT < handle
    % Contains methods pertinent to solving the 2D Hubbard-Holstein Model
    % according to Marsiglio's formulation
    properties
        kx
        ky
        k_pair_array
        E_kxkymesh
        w_m
        v_n
        t
        w_E
        lambda_0
        Beta
        n
        gwE
        tol
        Pi_vn0
        G0_hf
        D0
        tol_G
        tol_n
        tol_Lambda
        Sigma_trial_init
        renorm_calc
        U
        rx_phase_factor 
        ry_phase_factor
        jtau_phase_factor
        D_divergence_counter
    end
    
    methods
        function self = HolsteinU2D_FFT()
            disp('---------------------------- Holstein2D_FFT class Initialized ----------------------------')
        end
        
        function init_wavevectors(self,N_kx,N_ky)
            % Initialize T independent class properties
            % Call once per set of T iterations
            
            self.kx = (2*pi/N_kx).*linspace(-N_kx/2,N_kx/2,N_kx+1);
            self.kx = self.kx(2:end);
            % kx is the electron wavevector in the x direction in units of a
            % kx in (-N_kx/2, N_kx/2], where N_kx is the number of sites
            % in the x direction
            
            self.ky = (2*pi/N_ky).*linspace(-N_ky/2,N_ky/2,N_ky+1);
            self.ky = self.ky(2:end);
            % ky is the electron wavevector in the y direction in units of a
            % ky in (-N_ky/2, N_ky/2], where N_ky is the number of sites
            % in the y direction
            
            k_pair_array_ = cell(N_ky,N_kx);
            jy = 1;
            while jy <= N_ky
                jx = 1;
                while jx <= N_kx
                    k_pair_array_{jy,jx} = [self.kx(jx) self.ky(jy)];
                    jx = jx + 1;
                end
                jy = jy + 1;
            end
            self.k_pair_array = k_pair_array_;
            
        end
        
        function init_T_dep_objects(self,N_kx,N_ky,N_w,T)
            % initalize T dependent class properties
            % Call once per T iteration
            
            m = linspace(-N_w/2,N_w/2,N_w+1);
            % m in [-N_w, N_w] (N_w + 1 values of m)
            
            self.w_m = self.t*pi*T.*(2.*m-1);
            % w_m is the fermionic matsubara frequency
            
            self.v_n = self.t*pi*T.*(2.*m);
            % v_n is the bosonic matsubare frequency
            
            self.Beta = 1/T*self.t;
            % inverse temperature
            
            [kx_xmesh,ky_ymesh,w_m_zmesh] = meshgrid(self.kx,self.ky,self.w_m);
            self.E_kxkymesh = -2*self.t.*(cos(kx_xmesh) + cos(ky_ymesh));
            self.G0_hf = 1./(1i*w_m_zmesh - self.E_kxkymesh);
            % Construct G0 for mu = 0 (half-filling)
            
            [~,~,v_n_zmesh] = meshgrid(self.kx,self.ky,self.v_n);
            self.D0 = -2*self.w_E.*(1./(self.w_E^2 + (v_n_zmesh).^2));
            % Non-interacting phonon propagator
            
            
            N_w = length(self.w_m);
            [self.rx_phase_factor, self.ry_phase_factor, self.jtau_phase_factor] = ...
                meshgrid(exp(-1i*pi*(2/N_kx - 1).*(0:N_kx-1)),... % x real space phase factor
                exp(-1i*pi*(2/N_ky - 1).*(0:N_ky-1)),... % y real space phase factor
                exp(-1i*pi*(N_w-1)/N_w.*(0:N_w-1))); % imaginary time phase factor
            % Generate phase factors for convenient numerical
            % ordering with FFTs
            
        end
        
        
        function phase_diagram(self,w_cutoff,N_w_init,n_list,q_max_array,T_SP,T_CDW)
            
            %T_SP(T_CDW > T_SP) = NaN;
            %T_CDW(T_SP > T_CDW) = NaN;
            
            N_kx = length(self.kx);
            N_ky = length(self.ky);
            N_k = N_kx*N_ky;
            tol_ = self.tol_G;
            lambda_0_ = self.lambda_0;
            w_E_ = self.w_E;
            U_ = self.U;
            
            fig = figure;
            plot(n_list,T_SP,'-*b','linewidth',2,'DisplayName','$T_{SP}$')
            hold on
            plot(n_list,T_CDW,'-om','linewidth',2,'DisplayName','$T_{CDW}$')
            hold on
            for j = 1:length(n_list)
                if T_CDW(j) > T_SP(j)
                    if ~isempty(q_max_array{j})
                        q_x_pi = round(100*q_max_array{j}(1)/pi)/100;
                        q_y_pi = round(100*q_max_array{j}(2)/pi)/100;
                        if abs(q_x_pi - 1) < 1e-8 && abs(q_y_pi - 1) < 1e-8
                            string = '$(\pi, \pi)$';
                        elseif abs(q_x_pi - 1) < 1e-8 && abs(q_y_pi - 1) > 1e-8
                            string = ['$(\pi,',num2str(q_y_pi),'\pi)$'];
                        elseif abs(q_x_pi - 1) > 1e-8 && abs(q_y_pi - 1) < 1e-8
                            string = ['$(',num2str(q_x_pi),'\pi,','\pi)$'];
                        else
                            string = ['$(',num2str(q_x_pi),'\pi,',num2str(q_y_pi),'\pi)$'];
                        end
                        text(n_list(j),(T_CDW(j) + 0.5e-2),string,'FontSize',16);
                    end
                end
            end
            legend show
            set(legend,'location','best')
            xlabel '$n$'
            ylabel '$T_c$'
            xlim([0 1]);
            ylim([0 0.2]);
            %str = [num2str(N_kx),' $\times$ ',num2str(N_ky),' Sites, ','$\omega_{\mathrm{cutoff}}$ = $', num2str(w_cutoff/8),'W$'];
            %str = [str newline '$\lambda_0$ = ',num2str(lambda_0_),', $\omega_E$ = ',num2str(w_E_), ', $U$ = ', num2str(U_)];
            str = ['$\omega_E$ = $',num2str(w_E_),'$'];
            str = [str newline '$U$ = $',num2str(U_),'$'];
            str = [str newline '$\lambda_0$ = $',num2str(lambda_0_),'$'];
            str = [str newline '$N$ = $',num2str(N_kx),'^2$'];
            dim = [.15 .025 .3 .3];
            annotation('textbox',dim,'String',str,'FitBoxToText','on',...
                'interpreter','latex','FontSize',24,'linewidth',2,'backgroundcolor','w','EdgeColor','w');            
            set(gca,'XMinorTick','on','YMinorTick','on','linewidth',2,'FontSize',24)
            savefig(fig, ['2D_Phase_diagram_min_n=',num2str(min(n_list)),'_max_n=',num2str(max(n_list)),...
                '_Nk=',num2str(N_k),'_Nw=',num2str(N_w_init),'_lambda0=',num2str(lambda_0_),'_wE=',num2str(w_E_)...
                ,'_tol=',num2str(tol_),'_U=',num2str(U_),'.fig'])
        end
        
        
        function [T_SP, T_CDW] = linear_ising_fit_plot_Tc(self,w_cutoff,N_w_init,T_list,Chi_SP,Chi_CDW,q_max)
            % Linear fit to 1/Chi_SP
            % Ising fit to 1/Chi_CDW
            
            recip_SP = 1./Chi_SP(Chi_SP ~= 0);
            recip_CDW = 1./Chi_CDW(Chi_CDW ~= 0);
            
            % Linear Fit to T_SP near Tc
            T_list_near_Tc_SP = T_list(end-2:end).';
            recip_SP_near_Tc_SP = recip_SP(end-2:end);
            fit_SP = fit(T_list_near_Tc_SP,recip_SP_near_Tc_SP,'poly1',...
                'Upper',[inf 0],'Lower',[0 -inf]);
            coeff_vec = coeffvalues(fit_SP);
            T_SP = -coeff_vec(2)/coeff_vec(1);

            % Ising Fit to T_CDW near Tc
            %idx_T_0p1 = find(abs(T_list - 0.1) == min(abs(T_list - 0.1)));
            T_list_near_Tc_CDW = T_list(end-2:end).';
            recip_CDW_near_Tc_CDW = recip_CDW(end-2:end);
            starting_points = [1 T_list(end)];
            upper_bounds = [inf 1];
            lower_bounds = [0 1e-14];
            g = 7/4; % 2D Ising universality class critical exponent
            Ising_fit = fittype(@(A,Tc,x) A.*abs((x - Tc)).^g);
            fit_CDW = fit(T_list_near_Tc_CDW,recip_CDW_near_Tc_CDW,Ising_fit,...
                'Upper',upper_bounds,'Lower',lower_bounds,'StartPoint',starting_points);
            coeff_vec = coeffvalues(fit_CDW);
            T_CDW = coeff_vec(2);
            %g = 1;
            
            % Plot results
            N_kx = length(self.kx);
            N_ky = length(self.ky);
            N_k = N_kx*N_ky;
            tol_ = self.tol_G;
            lambda_0_ = self.lambda_0;
            w_E_ = self.w_E;
            U_ = self.U;
            n_ = self.n;
            T_list_fit = linspace(0.001,1,201);
            
            fig = figure;
            set(fig,'color','w','Units', 'Normalized', 'OuterPosition', [0 0 1/2 3/4])
            plot(T_list,recip_SP,'-','linewidth',2,'DisplayName','$1/\chi_{SP}$','Color','b')
            hold on
            plot(T_list_fit,fit_SP(T_list_fit),'--','linewidth',2,'DisplayName','$1/\chi_{SP}$ Linear Fit','Color','b')
            plot(T_list,recip_CDW,'-','linewidth',2,'DisplayName','$1/\chi_{CDW}(q_{max})$','Color','m')
            plot(T_list_fit,fit_CDW(T_list_fit),'--','linewidth',2,'DisplayName','$1/\chi_{CDW}(q_{max})$ Ising Fit','Color','m')
            
            % display system parameters in textbox
            q_x_pi = round(100*q_max(1)/pi)/100;
            q_y_pi = round(100*q_max(2)/pi)/100;
            if abs(q_x_pi - 1) < 1e-8 && abs(q_y_pi - 1) < 1e-8
                q_max_str = '$(\pi, \pi)$';
            elseif abs(q_x_pi - 1) < 1e-8 && abs(q_y_pi - 1) > 1e-8
                q_max_str = ['$(\pi,',num2str(q_y_pi),'\pi)$'];
            elseif abs(q_x_pi - 1) > 1e-8 && abs(q_y_pi - 1) < 1e-8
                q_max_str = ['$(',num2str(q_x_pi),'\pi,','\pi)$'];
            else
                q_max_str = ['$(',num2str(q_x_pi),'\pi,',num2str(q_y_pi),'\pi)$'];
            end
            dim = [.15 .025 .3 .3];
            str = [num2str(N_kx),' $\times$ ',num2str(N_ky),' Sites, ','$\omega_{\mathrm{cutoff}}$ = $', num2str(w_cutoff/8),'W$',', $n$ = ',num2str(n_),];
            str = [str newline '$\lambda_0$ = ',num2str(lambda_0_),', $\omega_E$ = ',num2str(w_E_),', $U$ = ',num2str(U_)];
            str = [str newline '$T_{SP}$ = $',num2str(T_SP),'$','; $T_{CDW}$ = $',num2str(T_CDW),'$',...
                '$\vec{q}_{\mathrm{max}}$ at $T^{\mathrm{CDW}}_{\mathrm{final}}$ = ',q_max_str];
            annotation('textbox',dim,'String',str,'FitBoxToText','on',...
                'interpreter','latex','FontSize',14,'linewidth',2,'backgroundcolor','w');            
            legend show
            xlabel '$T$'
            ylim([0, 5])
            xlim([0, 0.2])
            grid on
            set(gca,'XMinorTick','on','YMinorTick','on','linewidth',2,'FontSize',18)
            
            if isempty(self.renorm_calc)
                savefig(fig, ['Linear_Ising_Fits_Tc_recip_SP_vs_CDW_vs_T_n=',num2str(n_),'_Nk=',num2str(N_k),'_Nw=',num2str(N_w_init),...
                    '_lambda0=',num2str(lambda_0_),'_Tmin=',num2str(min(T_list)),'_wE=',num2str(w_E_),...
                    '_Tmax=',num2str(max(T_list)),'_tol=',num2str(tol_),'_U=',num2str(U_),'.fig'])
                
            else
                savefig(fig, ['Renorm_Linear_Ising_Fits_Tc_recip_SP_vs_CDW_vs_T_n=',num2str(n_),'_Nk=',num2str(N_k),'_Nw=',num2str(N_w_init),...
                    '_lambda0=',num2str(lambda_0_),'_Tmin=',num2str(min(T_list)),'_wE=',num2str(w_E_),...
                    '_Tmax=',num2str(max(T_list)),'_tol=',num2str(tol_),'_U=',num2str(U_),'.fig'])
                
            end
            
            %{
            if T_SP > T_CDW
                T_CDW = NaN;
            else
                T_SP = NaN;
            end
            %}
            
        end
        
        function [T_SP, T_CDW] = linear_linear_fit_plot_Tc(self,w_cutoff,N_w_init,T_list,Chi_SP,Chi_CDW,q_max)
            % Linear fit to 1/Chi_SP
            % Ising fit to 1/Chi_CDW
            
            recip_SP = 1./Chi_SP(Chi_SP ~= 0);
            recip_CDW = 1./Chi_CDW(Chi_CDW ~= 0);
            
            % Linear Fit to T_SP near Tc
            T_list_near_Tc_SP = T_list(end-3:end).';
            recip_SP_near_Tc_SP = recip_SP(end-3:end);
            fit_SP = fit(T_list_near_Tc_SP,recip_SP_near_Tc_SP,'poly1',...
                'Upper',[inf 0],'Lower',[0 -inf]);
            coeff_vec = coeffvalues(fit_SP);
            T_SP = -coeff_vec(2)/coeff_vec(1);

            % Ising Fit to T_CDW near Tc
            %idx_T_0p1 = find(abs(T_list - 0.1) == min(abs(T_list - 0.1)));
            T_list_near_Tc_CDW = T_list(end-3:end).';
            recip_CDW_near_Tc_CDW = recip_CDW(end-3:end);
            fit_CDW = fit(T_list_near_Tc_CDW,recip_CDW_near_Tc_CDW,'poly1',...
                'Upper',[inf 0],'Lower',[0 -inf]);
            coeff_vec = coeffvalues(fit_CDW);
            T_CDW = -coeff_vec(2)/coeff_vec(1);
            g = 1;
            
            % Plot results
            N_kx = length(self.kx);
            N_ky = length(self.ky);
            N_k = N_kx*N_ky;
            tol_ = self.tol_G;
            lambda_0_ = self.lambda_0;
            w_E_ = self.w_E;
            U_ = self.U;
            n_ = self.n;
            T_list_fit = linspace(0.001,1,201);
            
            fig = figure;
            set(fig,'color','w','Units', 'Normalized', 'OuterPosition', [0 0 1/2 3/4])
            plot(T_list,recip_SP,'-','linewidth',2,'DisplayName','$1/\chi_{SP}$','Color','b')
            hold on
            plot(T_list_fit,fit_SP(T_list_fit),'--','linewidth',2,'DisplayName','$1/\chi_{SP}$ Linear Fit','Color','b')
            plot(T_list,recip_CDW,'-','linewidth',2,'DisplayName','$1/\chi_{CDW}(\vec{q}_{max})$','Color','m')
            plot(T_list_fit,fit_CDW(T_list_fit),'--','linewidth',2,'DisplayName','$1/\chi_{CDW}(\vec{q}_{max})$ Linear Fit','Color','m')
            
            % display system parameters in textbox
            q_x_pi = round(100*q_max(1)/pi)/100;
            q_y_pi = round(100*q_max(2)/pi)/100;
            if abs(q_x_pi - 1) < 1e-8 && abs(q_y_pi - 1) < 1e-8
                q_max_str = '$(\pi, \pi)$';
            elseif abs(q_x_pi - 1) < 1e-8 && abs(q_y_pi - 1) > 1e-8
                q_max_str = ['$(\pi,',num2str(q_y_pi),'\pi)$'];
            elseif abs(q_x_pi - 1) > 1e-8 && abs(q_y_pi - 1) < 1e-8
                q_max_str = ['$(',num2str(q_x_pi),'\pi,','\pi)$'];
            else
                q_max_str = ['$(',num2str(q_x_pi),'\pi,',num2str(q_y_pi),'\pi)$'];
            end
            dim = [.15 .025 .3 .3];
            str = [num2str(N_kx),' $\times$ ',num2str(N_ky),' Sites, ','$\omega_{\mathrm{cutoff}}$ = $', num2str(w_cutoff/8),'W$',', $n$ = ',num2str(n_),];
            str = [str newline '$\lambda_0$ = ',num2str(lambda_0_),', $\omega_E$ = ',num2str(w_E_),', $U$ = ',num2str(U_)];
            str = [str newline '$T_{SP}$ = $',num2str(T_SP),'$','; $T_{CDW}$ = $',num2str(T_CDW),'$'];
            str = [str newline '$\vec{q}_{\mathrm{max}}$ at $T^{\mathrm{CDW}}_{\mathrm{final}}$ = ',q_max_str];
            annotation('textbox',dim,'String',str,'FitBoxToText','on',...
                'interpreter','latex','FontSize',14,'linewidth',2,'backgroundcolor','w');
            legend show
            xlabel '$T$'
            ylim([0, 1])
            xlim([0, 0.3])
            grid on
            set(gca,'XMinorTick','on','YMinorTick','on','linewidth',2,'FontSize',18)
            
            if isempty(self.renorm_calc)
                savefig(fig, ['Linear_Linear_Fits_Tc_recip_SP_vs_CDW_vs_T_n=',num2str(n_),'_Nk=',num2str(N_k),'_Nw=',num2str(N_w_init),...
                    '_lambda0=',num2str(lambda_0_),'_Tmin=',num2str(min(T_list)),'_wE=',num2str(w_E_),...
                    '_Tmax=',num2str(max(T_list)),'_tol=',num2str(tol_),'_U=',num2str(U_),'.fig'])
                
            else
                savefig(fig, ['Renorm_Linear_Linear_Ising_Fits_Tc_recip_SP_vs_CDW_vs_T_n=',num2str(n_),'_Nk=',num2str(N_k),'_Nw=',num2str(N_w_init),...
                    '_lambda0=',num2str(lambda_0_),'_Tmin=',num2str(min(T_list)),'_wE=',num2str(w_E_),...
                    '_Tmax=',num2str(max(T_list)),'_tol=',num2str(tol_),'_U=',num2str(U_),'.fig'])
                
            end
            
            %{
            if T_SP > T_CDW
                T_CDW = NaN;
            else
                T_SP = NaN;
            end
            %}
            
        end
        
        
        function [T_SP, T_CDW] = linear_fit_plot_Tc(self,N_w_init,T_list,Chi_SP,Chi_CDW_bar)
            % Linear fit to 1/Chi_SP and 1 - (lambda_0 - U)Chi_CDW
            
            recip_SP = 1./Chi_SP(~isnan(Chi_SP));
            denom_CDW = 1 - (self.lambda_0 - self.U).*Chi_CDW_bar(~isnan(Chi_CDW_bar));
            
            % Linear Fit to T_SP near Tc
            idx_SP = find(recip_SP == min(recip_SP));
            T_list_near_Tc_SP = T_list((ceil(idx_SP-3):idx_SP)).';
            recip_SP_near_Tc_SP = recip_SP((ceil(idx_SP-3):idx_SP));
            fit_SP = fit(T_list_near_Tc_SP,recip_SP_near_Tc_SP,'poly1',...
                'Upper',[inf 0],'Lower',[0 -inf]);           
            coeff_vec = coeffvalues(fit_SP);
            T_SP = -coeff_vec(2)/coeff_vec(1);

            % Linear Fit to T_CDW near Tc
            idx_CDW = find(denom_CDW == min(denom_CDW));
            T_list_near_Tc_CDW = T_list((ceil(idx_CDW-3):idx_CDW)).';
            denom_CDW_near_Tc_CDW = denom_CDW((ceil(idx_CDW-3):idx_CDW));
            fit_CDW = fit(T_list_near_Tc_CDW,denom_CDW_near_Tc_CDW,'poly1',...
                'Upper',[inf 0],'Lower',[0 -inf]);          
            coeff_vec = coeffvalues(fit_CDW);
            T_CDW = -coeff_vec(2)/coeff_vec(1);

            % Plot results
            N_k = length(self.k);
            tol_ = self.tol_G;
            lambda0_ = self.lambda_0;
            w_E_ = self.w_E;
            U_ = self.U;
            n_ = self.n;
            T_list_fit = linspace(0.001,1,201);

            fig = figure;
            set(fig,'color','w','Units', 'Normalized', 'OuterPosition', [0 0 1/2 3/4])
            plot(T_list,1./Chi_SP,'-','linewidth',2,'DisplayName','$1/\chi_{SP}$','Color','b')
            hold on
            plot(T_list_fit,fit_SP(T_list_fit),'--','linewidth',2,'DisplayName','$1/\chi_{SP}$ Fit','Color','b')
            plot(T_list,1 - (self.lambda_0 - self.U).*Chi_CDW_bar,'-','linewidth',2,'DisplayName','$1-(\lambda_0-U)\chi_{CDW}(q_{max})$','Color','m')
            plot(T_list_fit,fit_CDW(T_list_fit),'--','linewidth',2,'DisplayName','$1-(\lambda_0-U)\chi_{CDW}(q_{max})$ Fit','Color','m')
            dim = [.15 .025 .3 .3];
            str = [num2str(N_k),' Sites' newline num2str(N_w_init),' Matsubara Frequencies Min'];
            str = [str newline '$\lambda_0$ = ',num2str(lambda0_),', $n$ = ',num2str(n_)...
                ', $\omega_E$ = ',num2str(w_E_),', $U$ = ',num2str(U_)];
            str = [str newline '$T_{SP}$ = $',num2str(T_SP),'$','; $T_{CDW}$ = $',num2str(T_CDW),'$'];
            annotation('textbox',dim,'String',str,'FitBoxToText','on',...
                'interpreter','latex','FontSize',14,'linewidth',2,'backgroundcolor','w');
            legend show
            xlabel '$T$'
            ylim([0, 1])
            xlim([0, 0.1])
            grid on
            set(gca,'XMinorTick','on','YMinorTick','on','linewidth',2,'FontSize',18)
            
            if isempty(self.renorm_calc)
                savefig(fig, ['Linear_Fits_Tc_recip_SP_vs_CDW_vs_T_n=',num2str(n_),'_Nk=',num2str(N_k),'_Nw=',num2str(N_w_init),...
                    '_lambda0=',num2str(lambda0_),'_Tmin=',num2str(min(T_list)),'_wE=',num2str(w_E_),...
                    '_Tmax=',num2str(max(T_list)),'_tol=',num2str(tol_),'_U=',num2str(U_),'.fig'])
                
            else
                savefig(fig, ['Renorm_Linear_Fits_Tc_recip_SP_vs_CDW_vs_T_n=',num2str(n_),'_Nk=',num2str(N_k),'_Nw=',num2str(N_w_init),...
                    '_lambda0=',num2str(lambda0_),'_Tmin=',num2str(min(T_list)),'_wE=',num2str(w_E_),...
                    '_Tmax=',num2str(max(T_list)),'_tol=',num2str(tol_),'_U=',num2str(U_),'.fig'])
                
            end
            
            if T_SP > T_CDW
                T_CDW = NaN;
            else
                T_SP = NaN;
            end
            
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
            ylim([0, 2])
            xlim([0, 1])
            grid on
            set(gca,'XMinorTick','on','YMinorTick','on','linewidth',2,'FontSize',18)
            
            if isempty(self.renorm_calc)
                savefig(fig, ['SP_vs_CDW_vs_T_n=',num2str(n_),'_Nk=',num2str(N_k),'_Nw=',num2str(N_w_init),...
                    '_lambda0=',num2str(lambda0_),'_Tmin=',num2str(min(T_list)),'_wE=',num2str(w_E_),...
                    '_Tmax=',num2str(max(T_list)),'_tol=',num2str(tol_),'_U=',num2str(U_),'.fig'])
                
            else
                savefig(fig, ['Renorm_SP_vs_CDW_vs_T_n=',num2str(n_),'_Nk=',num2str(N_k),'_Nw=',num2str(N_w_init),...
                    '_lambda0=',num2str(lambda0_),'_Tmin=',num2str(min(T_list)),'_wE=',num2str(w_E_),...
                    '_Tmax=',num2str(max(T_list)),'_tol=',num2str(tol_),'_U=',num2str(U_),'.fig'])
                
            end
        end
        
        function [ydata_out, xdata_out] = acquire_figure_data(~,fig)
            dataObjs = findobj(fig,'-property','YData');
            ydata = cell(1,length(dataObjs));
            for j = 1:length(dataObjs)
                ydata{j} = dataObjs(j).YData;
            end
            ydata_out = ydata;
            xdata_out = dataObjs.XData;
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
            set(gca,'XMinorTick','on','YMinorTick','on','linewidth',2,'FontSize',18)
            
            if isempty(self.renorm_calc)
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
            set(gca,'XMinorTick','on','YMinorTick','on','linewidth',2,'FontSize',18)
            
            if isempty(self.renorm_calc)
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
        
        function [Chi_out, qx_max, qy_max] = Chi_CDW_bar_maxq(self)
            Pi_vn0_ = self.Pi_vn0;
            Chi_CDW_bar = -1/(self.gwE^2).*real(Pi_vn0_);
            q_max_cond = (Chi_CDW_bar == max(Chi_CDW_bar,[],'all'));
            Chi_out = mean(Chi_CDW_bar(q_max_cond));
            q_max_pair_list = cell2mat(self.k_pair_array(q_max_cond));
            qx_max = mean(abs(q_max_pair_list(:,1)));
            qy_max = mean(abs(q_max_pair_list(:,2)));
        end
        
        function [Chi_out, qx_max, qy_max] = Chi_CDW_maxq(self)
            Pi_vn0_ = self.Pi_vn0;
            Chi_CDW_bar = -1/(self.gwE^2).*real(Pi_vn0_);
            q_max_cond = (Chi_CDW_bar == max(Chi_CDW_bar,[],'all'));
            Chi_CDW_ = Chi_CDW_bar./(1-(self.lambda_0 - self.U).*Chi_CDW_bar);
            Chi_out = mean(Chi_CDW_(q_max_cond));
            q_max_pair_list = cell2mat(self.k_pair_array(q_max_cond));
            qx_max = mean(abs(q_max_pair_list(:,1)));
            qy_max = mean(abs(q_max_pair_list(:,2)));
        end
        
        
        function out = Chi_SP(self,F,Lambda)
            if max(Lambda,[],'all') == inf
                out = inf;
            else
                beta = self.Beta;
                N_k = length(self.kx)*length(self.ky);
                out = real(1/(N_k*beta)*sum(F.*Lambda, [1 2 3]));
            end
        end
        
        function out = converge_Lambda(self,F,D)
            disp('--------------------- Begin Lambda Convergence ----------------------')
            N_kx = length(self.kx);
            N_ky = length(self.ky);
            N_k = N_kx*N_ky;
            beta = self.Beta;
            gw_E = self.gwE;
            Lambda_factor = gw_E^2/(beta*N_k);
            U_ = self.U;
            Lambda_trial = ones(size(F));
            Lambda_calc = [];
            F_Lambda_trial = F.*Lambda_trial;
            trial_diff_1 = 1;
            counter = 0;
            divergence_counter = 0;
            divergence_lim = 100;
            iter_lim = 100000;  
            rx_phase = self.rx_phase_factor;
            ry_phase = self.ry_phase_factor;
            jtau_phase = self.jtau_phase_factor;
            
            while abs(trial_diff_1) > self.tol_Lambda && divergence_counter < divergence_lim && counter < iter_lim
                counter = counter + 1;
                Lambda_calc = 1 - (Lambda_factor).*(ifftn(rx_phase.*ry_phase.*conj(jtau_phase)...
                    .*fftn(D + U_/(gw_E^2)).*fftn(F_Lambda_trial)));
                trial_diff_2 = max(abs(Lambda_calc-Lambda_trial),[],'all');
                if abs(trial_diff_2) > abs(trial_diff_1)
                    divergence_counter = divergence_counter + 1;
                else
                    divergence_counter = 0;
                end
                trial_diff_1 = trial_diff_2;
                Lambda_trial = Lambda_calc;
                F_Lambda_trial = F.*Lambda_trial;
                disp([num2str(counter),' Iterations Completed to Converge Lambda'])
                disp(['|trial_diff|=',num2str(abs(trial_diff_1))])
            end

        
        if divergence_counter >= divergence_lim || counter >= iter_lim
            Lambda_calc(:) = inf;
            disp('--------------------- Lambda Has Diverged ! ----------------------')
        end
        
        out = Lambda_calc;
        disp('--------------------- End Lambda Convergence ----------------------')
        
    end
    
        function [G_out, Sigma_out, D_out, mu_out] = converge_mu(self,mu_1_guess,mu_2_guess)
            
            disp('--------------------- Begin mu, G Convergence ----------------------')
            beta = self.Beta;
            tol_n_ = self.tol_n;
            n_ = self.n;
            N_k = length(self.kx)*length(self.ky);
            E_k_matrix = self.E_kxkymesh(:,:,1);
            G0_hf_ = self.G0_hf;
            
            % Check for half-filling
            if n_ == 1
                mu_out = 0;
                [G_out, Sigma_out, D_out] = self.converge_G(G0_hf_,mu_out);
                Pi_out = self.Pi(G_out);
                self.Pi_vn0 = Pi_out(:,:,ceil(end/2));
                return
            end
            
            %First Guess
            mu_1 = mu_1_guess;
            n0_1 = (2/N_k)*sum(1./(exp(beta.*(E_k_matrix - mu_1)) + 1), [1 2]);
            G0_1 = 1./(1./G0_hf_ + mu_1);
            [G_1, Sigma_mk, D_nq] = self.converge_G(G0_1,mu_1);
            if isnan(G_1)
                G_out = NaN;
                Sigma_out = NaN;
                D_out = NaN;
                mu_out = NaN;
                return
            end
            n_1 = n0_1 + (2/(beta*N_k))*sum(G_1-G0_1, 1:3);
            delta_n_1 = real(n_ - n_1);
            
            mu_counter = 1;
            disp([num2str(mu_counter),' Iterations Completed to Converge mu'])
            disp(['|delta_n|=',num2str(abs(delta_n_1))])
            
            if abs(delta_n_1) < tol_n_
                G_out = G_1;
                Sigma_out = Sigma_mk;
                D_out = D_nq;
                mu_out = mu_1;
                Pi_out = self.Pi(G_out);
                self.Pi_vn0 = Pi_out(:,:,ceil(end/2));
                return
            end
            
            % Second Guess
            mu_2 = mu_2_guess;
            n0_2 = (2/N_k)*sum(1./(exp(beta.*(E_k_matrix - mu_2)) + 1), [1 2]);
            G0_2 = 1./(1./G0_hf_ + mu_2);
            [G_2, Sigma_mk, D_nq] = self.converge_G(G0_2,mu_2);
            if isnan(G_2)
                G_out = NaN;
                Sigma_out = NaN;
                D_out = NaN;
                mu_out = NaN;
                return
            end
            n_2 = n0_2 + (2/(beta*N_k))*sum(G_2-G0_2, 1:3);
            delta_n_2 = real(n_ - n_2);
            
            mu_counter = mu_counter + 1;
            disp([num2str(mu_counter),' Iterations Completed to Converge mu'])
            disp(['|delta_n|=',num2str(abs(delta_n_2))])
            
            if abs(delta_n_2) < tol_n_
                G_out = G_2;
                Sigma_out = Sigma_mk;
                D_out = D_nq;
                mu_out = mu_2;
                Pi_out = self.Pi(G_out);
                self.Pi_vn0 = Pi_out(:,:,ceil(end/2));
                return
            end
            
            iter_lim = 25;
            while abs(delta_n_2) > tol_n_ && mu_counter <= iter_lim
                if abs(delta_n_2 - delta_n_1) < 1e-14
                    delta_n_2 = delta_n_1 + 1;
                end
                % Apply Discrete Newton-Raphson to improve (3+j)rd guess for mu
                mu_3 = mu_2 - (mu_2-mu_1)*(delta_n_2)/(delta_n_2-delta_n_1);
                
                % (3+j)rd Guess
                n0_3 = (2/N_k)*sum((1./(exp(beta.*(E_k_matrix - mu_3)) + 1)), [1 2]);
                G0_3 = 1./(1./G0_hf_ + mu_3);
                [G_3, Sigma_mk, D_nq] = self.converge_G(G0_3,mu_3);
                if isnan(G_3)
                    G_out = NaN;
                    Sigma_out = NaN;
                    D_out = NaN;
                    mu_out = NaN;
                    return
                end
                n_3 = n0_3 + (2/(beta*N_k))*sum(G_3-G0_3, 1:3);
                delta_n_3 = real(n_ - n_3);
                
                % Recursively assign values to previous guesses
                mu_1 = mu_2;
                mu_2 = mu_3;
                delta_n_1 = delta_n_2;
                delta_n_2 = delta_n_3;
                
                mu_counter = mu_counter + 1;
                disp([num2str(mu_counter),' Iterations Completed to Converge mu'])
                disp(['|delta_n|=',num2str(abs(delta_n_2))])
            end
            
            if mu_counter > iter_lim               
                G_out = NaN;
                Sigma_out = NaN;
                D_out = NaN;
                mu_out = NaN;
                disp('-------------- mu Convergence Failed! Outputting NaNs ---------------')
                return
            end
            
            G_out = G_3;
            Sigma_out = Sigma_mk;
            D_out = D_nq;
            mu_out = mu_2;
            Pi_out = self.Pi(G_out);
            self.Pi_vn0 = Pi_out(:,:,ceil(end/2));
            disp('--------------------- End mu, G Convergence ----------------------')
            
        end
        
        function [G_mk, Sigma_mk, D_qn] = converge_G(self,G0,mu)
            disp('--------------------- Begin G Convergence ----------------------')
            % Utilizes Anderson mixing
            % Condition for CDW_bar used to identify CDW divergence
            % CDW_bar_cond = @(x)(max(-((self.lambda_0 - self.U)...
                %/self.gwE^2).*real(x(ceil(end/2),:))));
            
            % Iteration objects
            counter = 1;
            tol_ = self.tol_G;
            
            % Anderson mixing objects
            M = 3; % number of previous iterations to be mixed
            res_cell = cell(1,M);
            trial_cell = cell(1,M);
            
            % Non-int Phonon propagator
            D0_nq = self.D0;
            
            % 'Non-int' phonon self-energy
            %Pi0_nq = self.Pi0(mu);
            
            % 'Non-int' electron self-energy
            Sigma0_km = self.Sigma0(mu);
            
            % Construct Sigma0 as guess for first T iteration)
            Sigma_trial_1 = self.Sigma_trial_init;
            if isempty(Sigma_trial_1)
                Sigma_trial_1 = Sigma0_km;
            elseif ~isempty(Sigma_trial_1)
                % pad previous temperature iteration with Sigma0 if
                % number of matsubara frequencies are different
                N_w_prev = length(Sigma_trial_1(1,1,:));
                N_w = length(self.w_m);
                %zero_matrix = zeros(size(D0_nq));
                Sigma_trial_1 = cat(3,Sigma0_km(:,:,1:round((N_w - N_w_prev)/2)), Sigma_trial_1,...
                    Sigma0_km(:,:,(round((N_w + N_w_prev)/2 + 1)):N_w));
            end
            
            % First Guess
            self.D_divergence_counter = 0;
            Sigma_trial = Sigma_trial_1;
            trial_cell{1} = Sigma_trial;
            G_trial = 1./(1./G0 - Sigma_trial);
            [Sigma_calc, D_kk_mm] = self.Sigma(G_trial,G0,D0_nq,Sigma0_km);
            res_1 = Sigma_calc - Sigma_trial;
            res_cell{1} = res_1;
            trial_diff_1 = max(res_1,[],'all');
            disp([num2str(counter),' Iterations Completed to Converge G'])
            disp(['|trial_diff|=',num2str(abs(trial_diff_1))])
            
            % Initial damping factor
            b = 0.1;
            
            iter_lim = 500;
            % (2 + j)th Guess
            while abs(trial_diff_1) > tol_ && counter <= iter_lim
                counter = counter + 1;
                if counter <= M
                    Sigma_trial = Sigma_trial + b.*res_1;
                    trial_cell{counter} = Sigma_trial;
                    G_trial = 1./(1./G0 - Sigma_trial);
                    [Sigma_calc, D_kk_mm] = self.Sigma(G_trial,G0,D0_nq,Sigma0_km);
                    res_2 = Sigma_calc - Sigma_trial;
                    res_cell{counter} = res_2;
                    trial_diff_2 = max(res_2,[],'all');
                    b = b/(1 - trial_diff_2/trial_diff_1);
                    trial_diff_1 = trial_diff_2;
                    res_1 = res_2;
                    
                elseif counter == M + 1
                    Sigma_trial = Sigma_trial + b.*res_1;
                    G_trial = 1./(1./G0 - Sigma_trial);
                    [Sigma_calc, D_kk_mm] = self.Sigma(G_trial,G0,D0_nq,Sigma0_km);
                    res_2 = Sigma_calc - Sigma_trial;
                    trial_diff_2 = max(res_2,[],'all');
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
                    G_trial = 1./(1./G0 - Sigma_trial_1);
                    [Sigma_calc, D_kk_mm] = self.Sigma(G_trial,G0,D0_nq,Sigma0_km);
                    res_2 = Sigma_calc - Sigma_trial_1;
                    trial_diff_2 = max(res_2,[],'all');
                    trial_diff_1 = trial_diff_2;
                    res_1 = res_2;
                    Sigma_trial = Sigma_trial_1;
                    if abs(trial_diff_1) < 1e-2
                        iter_lim = iter_lim + 1;
                    end
                end
                
                disp([num2str(counter),' Iterations Completed to Converge G'])
                disp(['|trial_diff|=',num2str(abs(trial_diff_1))])
                
            end
            
            if counter > iter_lim
                disp('Convergence Failed! Electron Self-Energy Too Large!')
                disp(['max(D_nq) = ',num2str(max(abs(D_kk_mm),[],'all'))]);

                G_trial = NaN;
                Sigma_trial = NaN;
                D_kk_mm = NaN;
            end
            
            G_mk = G_trial;
            Sigma_mk = Sigma_trial;
            D_qn = D_kk_mm;
            disp('--------------------- End G Convergence ----------------------')
            
        end
        
        function out = anderson_mixing(self, b, prev_res_cell, new_res, prev_trial_cell, new_trial)
            M = length(prev_res_cell);
            N_kx = length(self.kx);
            N_ky = length(self.ky);
            
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
            
            out = reshape(new_trial_bar + b.*new_res_bar,N_ky,N_kx,[]);
            
        end
        
        function [Sigma_km, D_out]= Sigma(self,G_km,G0_km,D0_qn,Sigma0_km)
            % FUNCTION DESCRIPTION ----------------------------------------------------
            % Calculates electron self-energy including or excluding phonon-self energy
            % ARGUMENTS ---------------------------------------------------------------
            % G_mk is the single-electron propagator under trial
            % OUTPUTS -----------------------------------------------------------------
            % Sigma_mk is the calculated electron self-energy
            % D_qn is the effective phonon propagator
            

            N_k = length(self.kx)*length(self.ky);
            beta = self.Beta;
            gw_E = self.gwE;
            lambda0 = self.lambda_0;
            
            rx_phase = self.rx_phase_factor;
            ry_phase = self.ry_phase_factor;
            jtau_phase = self.jtau_phase_factor;
            
            Pi_qn = (2*gw_E^2/(beta*N_k)).*ifftn(jtau_phase.*rx_phase.*ry_phase.*fftn(G_km).*fftn(conj(G_km)));
            % Calculate phonon self-energy
            
            chi_bar = -1/self.gwE^2*real(Pi_qn(:,:,ceil(end/2)));
            if ~isempty(chi_bar*lambda0 >= 1) && self.D_divergence_counter <= 15
                chi_bar(chi_bar*lambda0 >= 1) = 1 - 1e-6;
                Pi_qn(:,:,ceil(end/2)) = -self.gwE^2*chi_bar;
                self.D_divergence_counter = self.D_divergence_counter + 1;
            end
            % Remove accidental divergences
            
            D_qn = 1./(1./D0_qn - Pi_qn);
            D_out = D_qn;
            % Calculate single phonon propagator
            Sigma_km = Sigma0_km + -(gw_E^2/(beta*N_k)).*(ifftn(conj(jtau_phase).*rx_phase.*ry_phase.*...
                (fftn(D_qn).*fftn(G_km) - fftn(D0_qn).*fftn(G0_km))));
            % Calculate electron self-energy
            
            
        end
        
        function out = Pi(self,G_km)
            N_k = length(self.kx)*length(self.ky);
            beta = self.Beta;
            gw_E = self.gwE;
            U_ = self.U;
            
            rx_phase = self.rx_phase_factor;
            ry_phase = self.ry_phase_factor;
            jtau_phase = self.jtau_phase_factor;
            
            Chi_qn = (1/(beta*N_k)).*ifftn(jtau_phase.*rx_phase.*ry_phase.*fftn(G_km).*fftn(conj(G_km)));
            % Calculate electron-electron interaction function
            
            out = 2*gw_E^2.*Chi_qn./(1-U_.*Chi_qn);
            % Calculate phonon self-energy
        end
        
        function out = Sigma0(self,mu)
            N_kx = length(self.kx);
            N_ky = length(self.ky);
            N_k = N_kx*N_ky;
            [~,~,w_m_zmesh] = meshgrid(self.kx,self.ky,self.w_m);
            beta = self.Beta;
            Z_k = self.E_kxkymesh - mu;
            w_E_ = self.w_E;
            BE_w_E = 1./(exp(beta.*w_E_) - 1);
            FD_Z_k = 1./(exp(beta.*Z_k) + 1);
            Sigma0_mkk = (1 + BE_w_E - FD_Z_k)./(1i.*w_m_zmesh - w_E_ - Z_k)...
                + (BE_w_E + FD_Z_k)./(1i.*w_m_zmesh + w_E_ - Z_k);
            Sigma0_m = sum(Sigma0_mkk,[1 2]);
            out = (self.gwE^2/N_k).*repmat(Sigma0_m,[N_ky N_kx 1]);
        end
        

        
    end
end

