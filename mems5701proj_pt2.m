%% Part 2
% pi_fan_initial = 0.7 % Initial guess for compressor pressure ratio
% A0_initial = 0
% net_thrust = 0
%
% pi_fan = pi_fan_initial
% A0 = A0_initial

% Cp_turbine = ???

% rho0 = P0 / (T0 * Rair)
% u0 = M0 * sqrt(gamma * Rair * T0)

% eta_fan = 0.95
% eta_turbine = 0.9
% eta_shaft = 0.95

% Iterate over fan pressure ratio

    % Pt0 = P0 / PrixM(M0, gamma)
    % Pt2 = pi_fan * Pt0

    % Pt15 = Pt2
    % P15 = PrixM(M15, gamma) * Pt15

    % Px = P15
    % Ptx = Px / PrixM(1, gamma)
    % pi_turbine = Ptx / Ptc
    % tau_turbine = 1 - eta_turbine * (1 - pow(pi_turbine, (gamma-1) / gamma ))

    % Ttx = tau_turbine * Ttc
    % specific_power = Cpt * (Ttc - Ttx)

    % Iterate over freestream area

        % m_dot0 = rho0 * u0 * A0

        % tau_fan = 1 + (1/eta_fan) * (pow(pi_fan, (gamma-1)/gamma) - 1)
        % Tt12 = Tt0
        % Tt13 = tau_fan * Tt12

        % power_fan = m_dot0 * Cpfan * (Tt13 - Tt12)
        % m_dot_turbine = power_fan / (eta_shaft * specific_power)