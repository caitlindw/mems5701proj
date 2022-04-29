%% Part 2


%% Constants
Ptc = 2e5;
Ttc = 540; 

sonic_mft = 0.6847; % Mft parameter at M=1 for gamma=1.4

MWT_air = MWT_yHyOyN(0, 0.232, 0.768);
R0 = 5.97994*(8314.4598/MWT_air); %gas constant / mixture avg molecular weight (slug/Kmol)

%Cp_H2 = (14.31*1000)*5.97994; %lbf/slug-R
Cp_H2 = 87010;
%Cp_O2 = (0.918*1000)*5.97994; %lbf/slug-R
Cp_O2 = 5438;

eta_fan = 0.95;
eta_turbine = 0.9;
eta_shaft = 0.95;

Beta = 0;
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

M15 = 0.2; % todo: make this separate for cruise and launch

Fn = 0;
FnMin = F_np_min_cruise;

phi_index = 0;
pi_fan_index = 0;


%% Iterate over equivalence ratio
for phi = [1, 2, 10]

    % Iterate over fan pressure ratio
    for pi_fan = 1.1:0.1:4.1
    
        Pt0 = P0 / PrixM(M0, gamma);
        Tt0 = T0 / TrixM(M0, gamma);
        Pt2 = pi_fan * Pt0;
    
        Pt15 = Pt2;
        P15 = PrixM(M15, gamma) * Pt15;

        Tt15 = Tt0;
        
        Px = P15;
        Ptx = Px / PrixM(1, gamma);
        Mx = MxPri(Px/Ptx, gamma);

        pi_turbine = Ptx / Ptc;
    
        %tau_turbine = 1 - eta_turbine * (1 - pow(pi_turbine, (gamma-1) / gamma ));
        tau_turbine = 1 - eta_turbine * (pi_turbine^( (gamma-1) / gamma ));
        Ttx = tau_turbine * Ttc;
        Tx = Ttx * TrixM(Mx, gamma);

        % Initialize Fn loop
        A0 = 0.01;
        counter = 0;
    
        while Fn < FnMin && counter < 10000

            fprintf('hello there');
    
            u0 = M0 * sqrt(gamma*R0*T0);
            m_dot_0 = Rho0 * u0 * A0;
            
            % Station 15
            A15 = (m_dot_0 * sqrt(R0 * Tt15)) / (Pt15 * MftxM(M15, gamma));
    
            %tau_fan = 1 + (1/eta_fan) * (pow(pi_fan, (gamma-1)/gamma) - 1);
            tau_fan = 1 + (1/eta_fan) * ( (pi_fan^(gamma-1)/gamma) - 1);
            Tt12 = Tt0;
            Tt13 = tau_fan * Tt12;

            % Station t/x
            m_dot_x = m_dot_turbine;

            % todo: fix circular reference in this section
            [yHON] = MassFracs(phi,Beta);
            yH_p = yHON(1);
            yO_p = yHON(2);
            yN_p = yHON(3);
            m_dot_H2_turbine = m_dot_x*yH_p;
            m_dot_O2_turbine = m_dot_x*yO_p;
            MWT_c = MWT_yHyOyN(yH_p, yO_p, yN_p);
            Rc = 5.97994*(8314.4598/MWT_c); % (slug/Kmol)
            Rx = Rc;
            Cp_c = ((gamma)/(gamma-1))*Rc; %lbf/slug-R 

            specific_power = Cp_c * (Ttc - Ttx);
            power_fan = m_dot0 * Cpfan * (Tt13 - Tt12);
            m_dot_turbine = power_fan / (eta_shaft * specific_power);
            % end todo

            ht_H2_turbine = Cp_H2 * Ttc * m_dot_H2_turbine;
            ht_O2_turbine = Cp_O2 * Ttc * m_dot_O2_turbine;

            ht_initial_turbine = (ht_H2_turbine + ht_O2_turbine) / m_dot_x;

            % Find the total temp and total pressure at the chamber exit
            % nozzle
            Ttx_turbine = Ttbrn_yHyOyNhi(yH_p, yO_p, yN_p, ht_initial_turbine, gamma);
            Ptx_turbine = Ptc;

            Ax = (m_dot_x * sqrt(Rx * Ttx_turbine)) / (Ptx_turbine * sonic_mft);
            Vx = Mx * sqrt(Rx * gamma * Tx);

            Ix = m_dot_x * Vx + Px * Ax;

            
            % With t/x and 15 fully defined, calculate station 7 properties
            % as an impulsive mixer
            A7 = Ax + A15;
            
            
            
            A0 = A0 + 0.01;
            counter = counter + 1;
            % end while loop
    
        end

        pi_fan_index = pi_fan_index + 1;
    
    end

    phi_index = phi_index + 1;

end


