clear; clc; close all;

gi = 6; ge = 2; ga = 1;
m = 9.1093837e-31;
h = 6.62607015e-34;
h_ = h / (2 * pi);
k_B = 1.380649e-23;
e = 1.60218e-19;
J1 = 15.6 * e;
J2 = 12.13 * e;
P_atm = 101325;
p = [1, 10];
R = 8.314;
Na = 6.02214076e23;
epsilon_0 = 8.8542e-12;

Temp = linspace(30000, 7500, 100);

f1 = @(T, J) (gi * ge / ga) * (m * k_B * T / (2 * pi * h_^2))^(3/2) .* exp(-J ./ (k_B * T)) / Na;
f2 = @(T, J) (gi * ge / ga) * (m * k_B * T / (2 * pi * h_^2))^(3/2) .* exp(-J ./ (k_B * T)) / Na;

debye_length = @(T, n_e) sqrt(epsilon_0 * k_B * T ./ (n_e * e^2));
delta_I = @(r, lambda_D) (r + 1) * e^2 / (4 * pi * epsilon_0 * lambda_D);
effective_ionization_energy = @(J, r, lambda_D) J - delta_I(r, lambda_D);
G = @(T, P, n) (e^2 / (4 * pi * epsilon_0 * k_B * T)) * (4 * pi * n * Na / 3)^(1/3);

n = [1e-6, 1e-6, 1e-6, 1e-6, 1e-6];

dat_ideal = cell(1, length(p));
dat_non_ideal = cell(1, length(p));
Temp_valid = cell(1, length(p));

for i = 1:length(p)
    for j = 1:length(Temp)
        T = Temp(j);
        try
            sol = fsolve(@(vars) [
                vars(1) * vars(2) - vars(4) * f1(T, J1);
                vars(1) * vars(3) - vars(5) * f2(T, J2);
                (p(i) * P_atm) / (T * R) - sum(vars);
                vars(1) - vars(2) - vars(3);
                vars(2) + vars(4) - 4 * vars(3) - 4 * vars(5)
            ], n, optimoptions('fsolve', 'Display', 'off'));
            
            n_ideal = sol;
            dat_ideal{i}(j) = G(T, p(i) * P_atm, n_ideal(1));
            
            lambda_D = debye_length(T, n_ideal(1) * Na);
            J1_eff = effective_ionization_energy(J1, 1, lambda_D);
            J2_eff = effective_ionization_energy(J2, 1, lambda_D);

            sol = fsolve(@(vars) [
                vars(1) * vars(2) - vars(4) * f1(T, J1_eff);
                vars(1) * vars(3) - vars(5) * f2(T, J2_eff);
                (p(i) * P_atm) / (T * R) - sum(vars);
                vars(1) - vars(2) - vars(3);
                vars(2) + vars(4) - 4 * vars(3) - 4 * vars(5)
            ], n_ideal, optimoptions('fsolve', 'Display', 'off'));

            n_non_ideal = sol;
            dat_non_ideal{i}(j) = G(T, p(i) * P_atm, n_non_ideal(1));
            Temp_valid{i}(j) = T;
        catch
            fprintf("⚠️ Solver failed at T = %.2f K, p = %.1f atm\n", T, p(i));
        end
    end
end

colors = {'b', 'r'};
markers = {'o', 's'};

figure; hold on;
for i = 1:length(p)
    plot(Temp_valid{i}, dat_ideal{i}, 'LineWidth', 1.5, 'Color', colors{i}, 'Marker', markers{i}, 'MarkerSize', 5, ...
         'DisplayName', sprintf('Ideal: p = %.1f atm', p(i)));
    plot(Temp_valid{i}, dat_non_ideal{i}, '--', 'LineWidth', 1.5, 'Color', colors{i}, 'Marker', markers{i}, 'MarkerSize', 5, ...
         'DisplayName', sprintf('Non-Ideal: p = %.1f atm', p(i)));
end

xlabel('Temperature (K)', 'FontSize', 12);
ylabel('Coulomb Coupling Parameter \Gamma', 'FontSize', 12);
title('Coulomb Coupling Parameter (Γ) vs Temperature', 'FontSize', 14);
legend('FontSize', 10);
grid on;
hold off;
