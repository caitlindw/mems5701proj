%% Part 2

% Cp_turbine = ???

% rho0 = P0 / (T0 * Rair)
% u0 = M0 * sqrt(gamma * Rair * T0)

%% Constants
Ptc = 2e5;
Ttc = 540; % todo: is this right for the chamber temperature?

MWT_air = MWT_yHyOyN(0, 0.232, 0.768);
R0 = 5.97994*(8314.4598/MWT_air); %gas constant / mixture avg molecular weight (slug/Kmol)

eta_fan = 0.95;
eta_turbine = 0.9;
eta_shaft = 0.95;


gamma = 1.4;

%% Cruise Conditions
M0_cruise       = 4.0;
rho0_cruise     = 2.26e-4;      % slug / cu_ft - Atmosphere Table
T0_cruise       = 389.97;       % R - Atmosphere Table
P0_cruise       = 1.049 * 144;  % lbf / sq_ft
F_np_min_cruise = 5000;     %lbf


%% Initialize pi_fan loop
M0 = M0_cruise;
P0 = P0_cruise;
T0 = T0_cruise;
Rho0 = rho0_cruise;
M15 = 0.2;

FnMin = F_np_min_cruise;

pi_fan_index = 0;

%% Iterate over fan pressure ratio
for pi_fan = 1.1:4.1:0.1

    Pt0 = P0 / PrixM(M0, gamma);
    Pt2 = pi_fan * Pt0;

    Pt15 = Pt2;
    P15 = PrixM(M15, gamma) * Pt15;

    Px = P15;
    Ptx = Px / PrixM(1, gamma);
    pi_turbine = Ptx / Ptc;

    tau_turbine = 1 - eta_turbine * (1 - pow(pi_turbine, (gamma-1) / gamma ));
    Ttx = tau_turbine * Ttc;
    specific_power = Cpt * (Ttc - Ttx);

    % Initialize Fn loop
    A0 = 0.01;
    counter = 0;

    while Fn < FnMin && counter < 10000

        u0 = M0 * sqrt(gamma*R0*T0);
        m_dot_0 = Rho0 * u0 * A0;

        tau_fan = 1 + (1/eta_fan) * (pow(pi_fan, (gamma-1)/gamma) - 1);
        Tt12 = Tt0;
        Tt13 = tau_fan * Tt12;

        power_fan = m_dot0 * Cpfan * (Tt13 - Tt12);
        m_dot_turbine = power_fan / (eta_shaft * specific_power);

        A0 = A0 + 0.01;
        counter = counter + 1;
        % end while loop

    end

end

    % Pt0 = P0 / PrixM(M0, gamma);
    % Pt2 = pi_fan * Pt0;

    % Pt15 = Pt2;
    % P15 = PrixM(M15, gamma) * Pt15;

    % Px = P15;
    % Ptx = Px / PrixM(1, gamma);
    % pi_turbine = Ptx / Ptc;
    % tau_turbine = 1 - eta_turbine * (1 - pow(pi_turbine, (gamma-1) /
    % gamma ));

    % Ttx = tau_turbine * Ttc
    % specific_power = Cpt * (Ttc - Ttx)

    % Iterate over freestream area

        % m_dot0 = rho0 * u0 * A0

        % tau_fan = 1 + (1/eta_fan) * (pow(pi_fan, (gamma-1)/gamma) - 1)
        % Tt12 = Tt0
        % Tt13 = tau_fan * Tt12

        % power_fan = m_dot0 * Cpfan * (Tt13 - Tt12)
        % m_dot_turbine = power_fan / (eta_shaft * specific_power)

