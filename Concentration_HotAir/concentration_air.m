clear; clc; close all;

% Constants
gi_N = 4; gi_O = 9; ga_N = 9; ga_O = 4; ge = 2;
m = 9.1093837e-31;
h_ = 1.0545718e-34;
k_B = 1.380649e-23;
J_N = 12.067 * 1.60218e-19;
J_O = 15.58 * 1.60218e-19;
P_atm = 101325;
p = 0.8;
R = 8.314;
Na = 6.02214076e23;

% Temperature range
Temp = linspace(40000, 7500, 100);

% Saha functions
f1 = @(T) (gi_N * ge / ga_N) * (m * k_B * T / (2 * pi * h_^2))^(3/2) .* exp(-J_N ./ (k_B * T)) / Na;
f2 = @(T) (gi_O * ge / ga_O) * (m * k_B * T / (2 * pi * h_^2))^(3/2) .* exp(-J_O ./ (k_B * T)) / Na;

% Initial guess
n = [1e-10, 1e-10, 1e-10, 1e-10, 1e-10];
dat = cell(1, length(n));

% System of equations
equations = @(vars, T) [
    vars(1) * vars(2) - vars(4) * f1(T);
    vars(1) * vars(3) - vars(5) * f2(T);
    (p * P_atm) / (T * R) - sum(vars);
    vars(1) - vars(2) - vars(3);
    vars(2) + vars(4) - 4 * vars(3) - 4 * vars(5)
];

% Solve numerically
for i = 1:length(Temp)
    T = Temp(i);
    sol = fsolve(@(vars) equations(vars, T), n, optimoptions('fsolve', 'Display', 'off'));
    if all(isfinite(sol))
        n = sol;
        for j = 1:length(n)
            dat{j}(i) = n(j) * Na;
        end
    else
        fprintf("⚠️ Warning: Solution did not converge at T = %.2f K\n", T);
    end
end

% Labels and colors
species_labels = {'Electrons (n_e)', 'Neutral Nitrogen (n_a_N)', 'Neutral Oxygen (n_a_O)', 'Ionized Nitrogen (n_i_N)', 'Ionized Oxygen (n_i_O)'};
colors = {'b', 'r', 'g', 'm', 'c'};
markers = {'o', 's', 'd', '^', 'v'};

% Plot results
figure; hold on;
for j = 1:length(n)
    plot(Temp, dat{j}, 'LineWidth', 1.5, 'Color', colors{j}, 'Marker', markers{j}, 'MarkerSize', 5, 'DisplayName', species_labels{j});
end

xlabel('Temperature (K)', 'FontSize', 12);
ylabel('Particle Concentration (molecules/m³)', 'FontSize', 12);
title('Ionization Equilibrium of Nitrogen and Oxygen in Air', 'FontSize', 14);
legend('FontSize', 10);
grid on;
hold off;
