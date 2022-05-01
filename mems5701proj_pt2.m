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

eta_fan     = 0.95;
eta_turbine = 0.9;
eta_shaft   = 0.99;



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
    
        %Flight Condition
        Pt0 = P0 / PrixM(M0, gamma);
        Tt0 = T0 / TrixM(M0, gamma);
        
        %Staton 12
        Pt12 = pi_inlet * Pt0;
        Tt12 = Tt0;
        
        %Station 13
        tau_fan = 1 + (1/eta_fan) * (pi_fan^((gamma-1)/gamma) - 1);
        Tt13 = tau_fan * Tt12;
        Tt15 = Tt13;
%         if pi_fan == 2.0 && phi == 2 
%             fprintf(['\n Tt12: ' num2str(Tt12)]);
%             fprintf(['\n Tt13: ' num2str(Tt13)]);
%             fprintf(['\n Tt15: ' num2str(Tt15)]);
%         end
        Pt13 = pi_fan * Pt12;

        %Station 15
        Pt15 = Pt13;
        P15  = PrixM(M15, gamma) * Pt15;
        
        Px = P15;
        Mx = 1; %given
        Ptx = Px / PrixM(Mx, gamma);

        pi_turbine = Ptx / Ptc; %Pt5/Pt4
    
        tau_turbine = 1 - eta_turbine * (1-pi_turbine^( (gamma-1) / gamma ));
        
        % Station t/x
        [yHON] = MassFracs(phi,Beta);
        yH_p = yHON(1);
        yO_p = yHON(2);
        yN_p = yHON(3);
        MWT_c = MWT_yHyOyN(yH_p, yO_p, yN_p);
        Rc = 5.97994*(8314.4598/MWT_c); % (slug/Kmol)
        Rx = Rc;
        Cp_c = ((gamma)/(gamma-1))*Rc; %lbf/slug-R 

        ht_initial_turbine = Ttc * Cp_c;
        Tt4 = Ttbrn_yHyOyNhi(yH_p, yO_p, yN_p, ht_initial_turbine, gamma);
        Tt5 = Tt4 * tau_turbine;
        Ttx = Tt5;

        Tx = Ttx * TrixM(Mx, gamma);

        power_turbine = Cp_c * (Tt4 - Tt5);       

        % Initialize Fn loop
        A0 = 0.01;
        counter = 0;
        Fnp = 0;
    
        while Fnp < FnpMin && counter < 10000
            
            %A0 = 2.11; % todo: remove this line
    
            % Freestream 
            u0      = M0 * sqrt(gamma*R0*T0);
            m_dot_0 = Rho0 * u0 * A0;
            ht_0    = Cp0 * Tt0 * m_dot_0;
            q0      = 0.5 * Rho0 * u0^2;
            
            % Station 15
            A15 = (m_dot_0 * sqrt(R0 * Tt15)) / (Pt15 * MftxM(M15, gamma));
            T15 = TrixM(M15, gamma) * Tt15;
            v15 = M15 * sqrt(R0 * gamma * T15);
            I15 = m_dot_0 * v15 + P15 * A15;
            
            power_fan = m_dot_0 * Cp0 * (Tt13 - Tt12);
            m_dot_turbine = abs(power_fan / (eta_shaft * power_turbine));

            % Find enthalpy of constituent species
            m_dot_x = m_dot_turbine;
            m_dot_H2_turbine = m_dot_x*yH_p;
            m_dot_O2_turbine = m_dot_x*yO_p;
            ht_H2_turbine = Cp_H2 * Ttc * m_dot_H2_turbine;
            ht_O2_turbine = Cp_O2 * Ttc * m_dot_O2_turbine;
    
            % Find the total temp and total pressure at the chamber exit
            % nozzle
            Ax = (m_dot_x * sqrt(Rx * Ttx)) / (Ptx * sonic_mft);
            Vx = Mx * sqrt(Rx * gamma * Tx);
            Ix = m_dot_x * Vx + Px * Ax;

            
            % With t/x and 15 fully defined, calculate station 7 properties
            % as an impulsive mixer
            A7 = Ax + A15;
            I7 = Ix + I15;
            m_dot_7      = m_dot_x + m_dot_0;
            ht_initial_7 = (ht_0 + ht_H2_turbine + ht_O2_turbine) / m_dot_7;

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

        if counter < 10000

            % Calculate performance parameters
            SFC = (m_dot_turbine * 32.17*3600) / Fnp;
            Isp = Fnp / (32.2 * m_dot_turbine);
            Fs  = Fnp / m_dot_0;

            if pi_fan == 2.0 && phi == 2 

                fprintf('\nFan: \n');
                fprintf(append('Pt12:           ', num2str(Pt12), '\n'));
                fprintf(append('Pt13:           ', num2str(Pt13), '\n'));
                fprintf(append('A0:             ', num2str(A0), '\n'));
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

                fprintf('\nTtb: \n');
                fprintf(append('Tt4:           ', num2str(Tt4), '\n'));

                fprintf('\nTurbine: \n');
                fprintf(append('tau_turbine:      ', num2str(tau_turbine), '\n'));
                fprintf(append('power_turbine:    ', num2str(power_turbine), '\n'));
                fprintf(append('Ttx:              ', num2str(Ttx), '\n'));
                fprintf(append('Tx:               ', num2str(Tx), '\n'));
                fprintf(append('Ptx:              ', num2str(Ptx), '\n'));
                fprintf(append('Px:               ', num2str(Px), '\n'));
                fprintf(append('Cp:               ', num2str(Cp_c), '\n'));
                fprintf(append('R:                ', num2str(Rc), '\n'));
                fprintf(append('m_dot_turbine:    ', num2str(m_dot_turbine), '\n'));
                fprintf(append('m_dot_H2_turbine: ', num2str(m_dot_H2_turbine), '\n'));
                fprintf(append('m_dot_O2_turbine: ', num2str(m_dot_O2_turbine), '\n'));
                fprintf(append('ht_H2_turbine:    ', num2str(ht_H2_turbine), '\n'));
                fprintf(append('ht_O2_turbine:    ', num2str(ht_O2_turbine), '\n'));

                fprintf('\nStation 7: \n');
                fprintf(append('m_dot_7:          ', num2str(m_dot_7), '\n'));
                fprintf(append('Pt7:              ', num2str(Pt7), '\n'));
                fprintf(append('Tt7:              ', num2str(Tt7), '\n'));
                fprintf(append('M7:               ', num2str(M7), '\n'));
                fprintf(append('ht7:              ', num2str(ht7), '\n'));
                fprintf(append('I7:               ', num2str(I7), '\n'));
                fprintf(append('A7:               ', num2str(A7), '\n'));
                fprintf(append('ht_initial_7:     ', num2str(ht_initial_7), '\n'));
                fprintf(append('Cp7:              ', num2str(Cp7), '\n'));
                fprintf(append('R7:               ', num2str(R7), '\n'));

                fprintf('\nStation 9: \n');
                fprintf(append('M9_ideal:              ', num2str(M9_ideal), '\n'));
                fprintf(append('T9_ideal:              ', num2str(T9_ideal), '\n'));
                fprintf(append('V9_ideal:              ', num2str(V9_ideal), '\n'));

                fprintf('\nPerformance: \n');
                fprintf(append('Fnp:            ', num2str(Fnp), '\n'));
                fprintf(append('SFC:            ', num2str(SFC), '\n'));
                fprintf(append('Isp:            ', num2str(Isp), '\n'));
                fprintf(append('Fs:             ', num2str(Fs), '\n'));

            end

        else

            fprintf('hit max iter');

        end

        pi_fan_index = pi_fan_index + 1;
    
    end

    phi_index = phi_index + 1;

end


