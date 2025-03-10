function ionizationpressure()
    kB = 1.380649e-23;       % Boltzmann constant (J/K)
    h  = 6.62607015e-34;     % Planck's constant (J*s)
    me = 9.10938356e-31;     % Electron mass (kg)
    eV_to_J = 1.602176634e-19;
    J_ion = 15.6 * eV_to_J;  % Ionization energy for Argon (J)
    e_charge = 1.602176634e-19; % Elementary charge (C)
    epsilon0 = 8.854187817e-12; % Vacuum permittivity (F/m)
    
    g_i = 6; g_e = 2; g_a = 1; % Statistical weights
    pressures_Pa = [0.1, 1, 10] * 101325; % Pressures in Pascals
    Tvals = linspace(1000, 40000, 200);
    
    gamma_analytic_vals = zeros(length(pressures_Pa), length(Tvals));
    gamma_numeric_vals = zeros(length(pressures_Pa), length(Tvals));
    
    figure;
    subplot(2,1,1); hold on; grid on;
    
    for P_idx = 1:length(pressures_Pa)
        P = pressures_Pa(P_idx);
        alpha_analytic = zeros(size(Tvals));
        alpha_numeric = zeros(size(Tvals));
        
        for iT = 1:length(Tvals)
            T = Tvals(iT);
            ST = (g_i * g_e / g_a) * (2 * pi * me * kB * T / h^2)^(3/2) * exp(-J_ion / (kB * T));
            
            ni_analytic = (-2 * ST + sqrt((2 * ST)^2 + 4 * ST * (P / (kB * T)))) / 2;
            n0_analytic = P / (kB * T) - 2 * ni_analytic;
            ne_analytic = ni_analytic;
            alpha_analytic(iT) = ni_analytic / (n0_analytic + ni_analytic);
            
            r_d_analytic = sqrt(epsilon0 * kB * T / (ne_analytic * e_charge^2));
            gamma_analytic_vals(P_idx, iT) = e_charge^2 / (4 * pi * epsilon0 * r_d_analytic * kB * T);
            
            guess = log([n0_analytic, ni_analytic, ni_analytic]);
            
            obj_fun = @(vars) sum(([exp(vars(3)) - exp(vars(2));
                                    (exp(vars(2)) * exp(vars(3))) / exp(vars(1)) - ST;
                                    exp(vars(1)) + exp(vars(2)) + exp(vars(3)) - P / (kB * T)]).^2);
            
            options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off');
            sol_numeric = fmincon(obj_fun, guess, [], [], [], [], [], [], [], options);
            
            n0_numeric = exp(sol_numeric(1));
            ni_numeric = exp(sol_numeric(2));
            ne_numeric = exp(sol_numeric(3));
            
            alpha_numeric(iT) = ni_numeric / (n0_numeric + ni_numeric);
            
            r_d_numeric = sqrt(epsilon0 * kB * T / (ne_numeric * e_charge^2));
            gamma_numeric_vals(P_idx, iT) = e_charge^2 / (4 * pi * epsilon0 * r_d_numeric * kB * T);
        end
        
        plot(Tvals, alpha_analytic, '-', 'DisplayName', sprintf('Analytic \alpha, P=%.1f atm', P/101325));
        plot(Tvals, alpha_numeric, 'o', 'DisplayName', sprintf('Numeric \alpha, P=%.1f atm', P/101325));
        xlabel('Temperature (K)'); ylabel('Ionization Degree \alpha');
        title('Ionization Degree Comparison: Analytic vs Numeric');
        legend('Location', 'best');
    end

    
    subplot(2,1,2); hold on; grid on;
    for P_idx = 1:length(pressures_Pa)
        plot(Tvals, gamma_analytic_vals(P_idx, :), '--', 'DisplayName', sprintf('Gamma Analytic, P=%.1f atm', pressures_Pa(P_idx)/101325));
        plot(Tvals, gamma_numeric_vals(P_idx, :), 'x', 'DisplayName', sprintf('Gamma Numeric, P=%.1f atm', pressures_Pa(P_idx)/101325));
    end
    xlabel('Temperature (K)'); ylabel('\Gamma');
    title('Gamma Comparison: Analytic vs Numeric');
    legend('Location', 'best');
    hold off;
end
