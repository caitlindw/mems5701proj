%% Part 2


%% Constants
Ptc = 2e5;
Ttc = 540; 

Beta = 0;
gamma = 1.4;
Cfg = 0.96;

sonic_mft = 0.6847; % Mft parameter at M=1 for gamma=1.4

MWT_air = MWT_yHyOyN(0, 0.232, 0.768);
R0 = 5.97994*(8314.4598/MWT_air); %gas constant / mixture avg molecular weight (slug/Kmol)
Cp0 = R0 * (gamma / (gamma-1));

%Cp_H2 = (14.31*1000)*5.97994; %lbf/slug-R
Cp_H2 = 87010;
%Cp_O2 = (0.918*1000)*5.97994; %lbf/slug-R
Cp_O2 = 5438;

eta_fan = 0.95;
eta_turbine = 0.9;
eta_shaft = 0.95;



%% Cruise Conditions
M0_cruise       = 4.0;
rho0_cruise     = 2.26e-4;      % slug / cu_ft - Atmosphere Table
T0_cruise       = 389.97;       % R - Atmosphere Table
P0_cruise       = 1.049 * 144;  % lbf / sq_ft
F_np_min_cruise = 5000;         %lbf
M15_cruise      = 0.2;
Cd_cruise       = 0.4;
pi_inlet_cruise = 0.669;


%% Initialize pi_fan loop
M0 = M0_cruise;
P0 = P0_cruise;
T0 = T0_cruise;
Rho0     = rho0_cruise;
M15      = M15_cruise;
Cd       = Cd_cruise;
pi_inlet = pi_inlet_cruise;

Fnp = 0;
FnpMin = F_np_min_cruise;

phi_index = 0;
pi_fan_index = 0;


%% Iterate over equivalence ratio
for phi = [1, 2, 10]

    % Iterate over fan pressure ratio
    for pi_fan = 1.1:0.1:4.1
    
        Pt0 = P0 / PrixM(M0, gamma);
        Tt0 = T0 / TrixM(M0, gamma);

        Pt12 = pi_inlet * Pt0;
        Tt12 = Tt0;

        tau_fan = 1 + (1/eta_fan) * (pi_fan^((gamma-1)/gamma) - 1);
        Tt13 = tau_fan * Tt12;
        Pt13 = pi_fan * Pt12;

        

        Pt15 = Pt13;
        P15  = PrixM(M15, gamma) * Pt15;

        Tt15 = Tt13;
        
        Px = P15;
        Ptx = Px / PrixM(1, gamma);
        Mx = MxPri(Px/Ptx, gamma);

        pi_turbine = Ptx / Ptc;
    
        %tau_turbine = 1 - eta_turbine * (1 - pow(pi_turbine, (gamma-1) / gamma ));
        tau_turbine = 1 - eta_turbine * (pi_turbine^( (gamma-1) / gamma ));
        Ttx = tau_turbine * Ttc;
        Tx = Ttx * TrixM(Mx, gamma);

        % Initialize Fn loop
        A0 = 1;
        counter = 0;
    
        while Fnp < FnpMin && counter < 1000

            % todo: remove this line
            A0 = 2.11;
    
            % Freestream 
            u0      = M0 * sqrt(gamma*R0*T0);
            m_dot_0 = Rho0 * u0 * A0;
            ht_0    = Cp0 * Tt0 * m_dot_0;
            q0      = 0.5 * Rho0 * u0^2;
            
            % Station 15
            A15 = (m_dot_0 * sqrt(R0 * Tt15)) / (Pt15 * MftxM(M15, gamma));
            T15 = TrixM(M15, gamma) * Tt15;
            fprintf(num2str(Tt15));
            fprintf('\n');
            fprintf(num2str(TrixM(M15, gamma) * Tt15));
            v15 = M15 * sqrt(R0 * gamma * T15);
            I15 = m_dot_0 * v15 + P15 * A15;
    
            %tau_fan = 1 + (1/eta_fan) * ( (pi_fan^((gamma-1)/gamma)) - 1);
            %Tt12 = Tt0;
            %Tt13 = tau_fan * Tt12;

            % Station t/x
            [yHON] = MassFracs(phi,Beta);
            yH_p = yHON(1);
            yO_p = yHON(2);
            yN_p = yHON(3);
            MWT_c = MWT_yHyOyN(yH_p, yO_p, yN_p);
            Rc = 5.97994*(8314.4598/MWT_c); % (slug/Kmol)
            Rx = Rc;
            Cp_c = ((gamma)/(gamma-1))*Rc; %lbf/slug-R 

            power_turbine = Cp_c * (Ttc - Ttx);
            power_fan = m_dot_0 * Cp0 * (Tt13 - Tt12);
            m_dot_turbine = power_fan / (eta_shaft * power_turbine);

            % Find enthalpy of constituent species
            m_dot_x = m_dot_turbine;
            m_dot_H2_turbine = m_dot_x*yH_p;
            m_dot_O2_turbine = m_dot_x*yO_p;
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
            I7 = Ix + I15;
            m_dot_7      = m_dot_x + m_dot_0;
            ht_initial_7 = (ht_0 * ht_H2_turbine * ht_O2_turbine) / m_dot_7;

            % Find the mass fractions of oxygen, hydrogen, and nitrogen in
            % the combined stream
            yH2_7 = m_dot_H2_turbine / m_dot_7;
            yO2_7 = (0.768 * m_dot_0 + m_dot_O2_turbine) / m_dot_7;
            yN2_7 = (0.232 * m_dot_0) / m_dot_7;
            MWT_7 = MWT_yHyOyN(yH2_7, yO2_7, yN2_7);
            R7    = 5.97994*(8314.4598/MWT_7);
            Cp7   = R7 * (gamma / (gamma-1));

            % Find TtBurned
            Tt7 = Ttbrn_yHyOyNhi(yH2_7, yO2_7, yN2_7, ht_initial_7, gamma);
            ht7 = Tt7 * Cp7;

            % Find Mfi parameter and Mach
            Mfi_7 = ht7 * (m_dot_7 / I7)^2;
            M7 = MbxMfi(Mfi_7, gamma);

            % Find Pt7 using Mft
            Pt7 = m_dot_7 * sqrt(R7 * Tt7) / (MftxM(M7, gamma) * A7);


            % Station 9

            % Assume these are true for isentropic, ideal expansion
            m_dot_9 = m_dot_7;
            Tt9 = Tt7;
            P9  = P0;
            Pt9_ideal = Pt7;
            R9_ideal  = R7;

            % Find ideal station 9 properties
            M9_ideal = MxPri(P9/Pt9_ideal, gamma);
            T9_ideal = TrixM(M9_ideal, gamma) * Tt9;
            V9_ideal = M9_ideal * sqrt(gamma * R9_ideal * T9_ideal);

            % Find actual gross thrust
            Fg_ideal = m_dot_9 * V9_ideal;
            Fg = Cfg * Fg_ideal;

            % Find ram drag and drag increments
            Dram = m_dot_0 * u0;

            Acowl = 0.1 * A0;
            Dcowl = q0 * Cd * Acowl;

            % Calc net propulsive force
            Fnp = Fg - Dram - Dcowl;

            
            A0 = A0 + 0.01;
            counter = counter + 1;
            % end while loop
    
        end

        if counter < 1000

            % todo: calc Fs, Isp, and SFC

            if pi_fan == 2.0 && phi == 2 

                fprintf('\nEngine: \n');
                fprintf(append('Fnp:            ', num2str(Fnp), '\n'));

                fprintf('\Fan: \n');
                fprintf(append('Pt12:           ', num2str(Pt12), '\n'));
                fprintf(append('Pt13:           ', num2str(Pt13), '\n'));
                fprintf(append('m_dot_0:        ', num2str(m_dot_0), '\n'));
                fprintf(append('power_fan:      ', num2str(power_fan), '\n'));
                fprintf(append('tau_fan:        ', num2str(tau_fan), '\n'));
                
                fprintf('\nStation 15: \n');
                fprintf(append('A15:           ', num2str(A15), '\n'));
                fprintf(append('M15:           ', num2str(M15), '\n'));
                fprintf(append('T15:           ', num2str(T15), '\n'));
                fprintf(append('Tt15:          ', num2str(Tt15), '\n'));
                fprintf(append('P15:           ', num2str(P15), '\n'));
                fprintf(append('Pt15:          ', num2str(Pt15), '\n'));
                fprintf(append('v15:           ', num2str(v15), '\n'));
                fprintf(append('I15:           ', num2str(I15), '\n'));


                fprintf('\nTurbine: \n');
                fprintf(append('power_turbine:  ', num2str(power_turbine), '\n'));
                fprintf(append('Tt13:           ', num2str(Tt13), '\n'));
                fprintf(append('m_dot_turbine:  ', num2str(m_dot_turbine), '\n'));
                fprintf(append('ht_H2_turbine:  ', num2str(ht_H2_turbine), '\n'));
                fprintf(append('ht_O2_turbine:  ', num2str(ht_O2_turbine), '\n'));

            end

        else

            fprintf('hit max iter');

        end

        pi_fan_index = pi_fan_index + 1;
    
    end

    phi_index = phi_index + 1;

end


