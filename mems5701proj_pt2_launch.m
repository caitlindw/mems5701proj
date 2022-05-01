%% Part 2
close all;

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

Cp_H2 = 87010;
Cp_O2 = 5438;

eta_fan     = 0.95;
eta_turbine = 0.9;
eta_shaft   = 0.99;

%% Launch Conditions
M0_launch       = 0.8;
rho0_launch     = 8.91e-4;      % slug / cu_ft - Atmosphere Table
T0_launch       = 411.84;       % R - Atmosphere Table
P0_launch       = 4.373 * 144;  % lbf / sq_ft
F_np_min_launch = 4000;
M15_launch      = 0.2;
Cd_launch       = 0.2;
pi_inlet_launch = 0.96;

%% Initialize pi_fan loop
M0 = M0_launch;
P0 = P0_launch;
T0 = T0_launch;

Rho0     = rho0_launch;
M15      = M15_launch;
Cd       = Cd_launch;
pi_inlet = pi_inlet_launch;
FnpMin  = F_np_min_launch;

phi_index = 1;
pi_fan_index = 1;
pi_fan_list = 1.1:0.1:4.1;

%Flight Condition
Pt0 = P0 / PrixM(M0, gamma);
Tt0 = T0 / TrixM(M0, gamma);

%Staton 12
Pt12 = pi_inlet * Pt0;
Tt12 = Tt0;

%% Initialize Matrices for Flight Values
tau_fan = zeros(3,31);
Tt13 = zeros(3,31);
Tt15 = zeros(3,31);
Pt13 = zeros(3,31);
Pt15 = zeros(3,31);
P15 = zeros(3,31);
Ptx = zeros(3,31);
Px = zeros(3,31); 
pi_turbine = zeros(3,31); 
tau_turbine = zeros(3,31); 
Rc = zeros(3,31); 
Rx = zeros(3,31); 
Cp_c = zeros(3,31); 
ht_initial_turbine = zeros(3,31); 
Tt4 = zeros(3,31);
Tt5 = zeros(3,31);
Ttx = zeros(3,31);
Tx = zeros(3,31);
power_turbine = zeros(3,31); 
m_dot_0 = zeros(3,31); 
ht_0 = zeros(3,31); 
A15 = zeros(3,31); 
T15 = zeros(3,31); 
v15 = zeros(3,31); 
I15 = zeros(3,31); 
power_fan = zeros(3,31);
m_dot_turbine = zeros(3,31);
m_dot_x = zeros(3,31);
m_dot_H2_turbine = zeros(3,31);
m_dot_O2_turbine = zeros(3,31);
ht_H2_turbine = zeros(3,31);
ht_O2_turbine = zeros(3,31);
Ax = zeros(3,31);
Vx = zeros(3,31);
Ix = zeros(3,31);
A7 = zeros(3,31);
I7 = zeros(3,31);
m_dot_7 = zeros(3,31);
ht_initial_7 = zeros(3,31);
yH2_7 = zeros(3,31);
yO2_7 = zeros(3,31);
yN2_7 = zeros(3,31);
MWT_7 = zeros(3,31);
R7 = zeros(3,31);
Cp7 = zeros(3,31);
Tt7 = zeros(3,31);
ht7 = zeros(3,31);
Mfi_7 = zeros(3,31);
M7 = zeros(3,31);
Pt7 = zeros(3,31);
m_dot_9 = zeros(3,31);
Tt9 = zeros(3,31);
Pt9_ideal = zeros(3,31);
R9_ideal = zeros(3,31);
M9_ideal = zeros(3,31);
T9_ideal = zeros(3,31);
V9_ideal = zeros(3,31);
Fg_ideal = zeros(3,31);
Fg = zeros(3,31);
Dram = zeros(3,31);
Acowl = zeros(3,31);
Dcowl = zeros(3,31);
Fnp = zeros(3,31);
SFC = zeros(3,31);
Isp = zeros(3,31);
Fs = zeros(3,31);
A0 = zeros(3,31);


%% Iterate over equivalence ratio
for phi = [1, 2, 10]

    % Iterate over fan pressure ratio
    for pi_fan = 1.1:0.1:4.1
        
        %Station 13
        tau_fan(phi_index,pi_fan_index) = 1 + (1/eta_fan) * (pi_fan^((gamma-1)/gamma) - 1);
        Tt13(phi_index,pi_fan_index) = tau_fan(phi_index,pi_fan_index) * Tt12;
        Tt15(phi_index,pi_fan_index) = Tt13(phi_index,pi_fan_index);
%         if pi_fan == 2.0 && phi == 2 
%             fprintf(['\n Tt12: ' num2str(Tt12)]);
%             fprintf(['\n Tt13: ' num2str(Tt13)]);
%             fprintf(['\n Tt15: ' num2str(Tt15)]);
%         end
        Pt13(phi_index,pi_fan_index) = pi_fan * Pt12;

        %Station 15
        Pt15(phi_index,pi_fan_index) = Pt13(phi_index,pi_fan_index);
        P15(phi_index,pi_fan_index)  = PrixM(M15, gamma) * Pt15(phi_index,pi_fan_index);
        
        Px(phi_index,pi_fan_index) = P15(phi_index,pi_fan_index);
        Mx = 1; %given
        Ptx(phi_index,pi_fan_index) = Px(phi_index,pi_fan_index) / PrixM(Mx, gamma);

        pi_turbine(phi_index,pi_fan_index) = Ptx(phi_index,pi_fan_index) / Ptc; %Pt5/Pt4
    
        tau_turbine(phi_index,pi_fan_index) = 1 - eta_turbine * (1-pi_turbine(phi_index,pi_fan_index)^( (gamma-1) / gamma ));
        
        % Station t/x
        [yHON] = MassFracs(phi_index,Beta);
        yH_p = yHON(1);
        yO_p = yHON(2);
        yN_p = yHON(3);
        MWT_c = MWT_yHyOyN(yH_p, yO_p, yN_p);
        Rc(phi_index,pi_fan_index) = 5.97994*(8314.4598/MWT_c); % (slug/Kmol)
        Rx(phi_index,pi_fan_index) = Rc(phi_index,pi_fan_index);
        Cp_c(phi_index,pi_fan_index) = ((gamma)/(gamma-1))*Rc(phi_index,pi_fan_index); %lbf/slug-R 

        ht_initial_turbine(phi_index,pi_fan_index) = Ttc * Cp_c(phi_index,pi_fan_index);
        Tt4(phi_index,pi_fan_index) = Ttbrn_yHyOyNhi(yH_p, yO_p, yN_p, ht_initial_turbine(phi_index,pi_fan_index), gamma);
        Tt5(phi_index,pi_fan_index) = Tt4(phi_index,pi_fan_index) * tau_turbine(phi_index,pi_fan_index);
        Ttx(phi_index,pi_fan_index) = Tt5(phi_index,pi_fan_index);

        Tx(phi_index,pi_fan_index) = Ttx(phi_index,pi_fan_index) * TrixM(Mx, gamma);

        power_turbine(phi_index,pi_fan_index) = Cp_c(phi_index,pi_fan_index) * (Tt4(phi_index,pi_fan_index) - Tt5(phi_index,pi_fan_index));       

        % Initialize Fn loop
        A0(phi_index,pi_fan_index) = 0.01;
        counter = 0;
        %Fnp = 0;
    
        while Fnp(phi_index,pi_fan_index) < FnpMin && counter < 10000

            % Freestream 
            u0      = M0 * sqrt(gamma*R0*T0);
            m_dot_0(phi_index,pi_fan_index) = Rho0 * u0 * A0(phi_index,pi_fan_index);
            ht_0(phi_index,pi_fan_index)    = Cp0 * Tt0 * m_dot_0(phi_index,pi_fan_index);
            q0      = 0.5 * Rho0 * u0^2;
            
            % Station 15
            A15(phi_index,pi_fan_index) = (m_dot_0(phi_index,pi_fan_index) * sqrt(R0 * Tt15(phi_index,pi_fan_index))) / (Pt15(phi_index,pi_fan_index) * MftxM(M15, gamma));
            T15(phi_index,pi_fan_index) = TrixM(M15, gamma) * Tt15(phi_index,pi_fan_index);
            v15(phi_index,pi_fan_index) = M15 * sqrt(R0 * gamma * T15(phi_index,pi_fan_index));
            I15(phi_index,pi_fan_index) = m_dot_0(phi_index,pi_fan_index) * v15(phi_index,pi_fan_index) + P15(phi_index,pi_fan_index) * A15(phi_index,pi_fan_index);
            
            power_fan(phi_index,pi_fan_index) = m_dot_0(phi_index,pi_fan_index) * Cp0 * (Tt13(phi_index,pi_fan_index) - Tt12);
            m_dot_turbine(phi_index,pi_fan_index) = abs(power_fan(phi_index,pi_fan_index) / (eta_shaft * power_turbine(phi_index,pi_fan_index)));

            % Find enthalpy of constituent species
            m_dot_x(phi_index,pi_fan_index) = m_dot_turbine(phi_index,pi_fan_index);
            m_dot_H2_turbine(phi_index,pi_fan_index) = m_dot_x(phi_index,pi_fan_index)*yH_p;
            m_dot_O2_turbine(phi_index,pi_fan_index) = m_dot_x(phi_index,pi_fan_index)*yO_p;
            ht_H2_turbine(phi_index,pi_fan_index) = Cp_H2 * Ttc * m_dot_H2_turbine(phi_index,pi_fan_index);
            ht_O2_turbine(phi_index,pi_fan_index) = Cp_O2 * Ttc * m_dot_O2_turbine(phi_index,pi_fan_index);
    
            % Find the total temp and total pressure at the chamber exit
            % nozzle
            Ax(phi_index,pi_fan_index) = (m_dot_x(phi_index,pi_fan_index) * sqrt(Rx(phi_index,pi_fan_index) * Ttx(phi_index,pi_fan_index))) / (Ptx(phi_index,pi_fan_index) * sonic_mft);
            Vx(phi_index,pi_fan_index) = Mx * sqrt(Rx(phi_index,pi_fan_index) * gamma * Tx(phi_index,pi_fan_index));
            Ix(phi_index,pi_fan_index) = m_dot_x(phi_index,pi_fan_index) * Vx(phi_index,pi_fan_index) + Px(phi_index,pi_fan_index) * Ax(phi_index,pi_fan_index);

            
            % With t/x and 15 fully defined, calculate station 7 properties
            % as an impulsive mixer
            A7(phi_index,pi_fan_index) = Ax(phi_index,pi_fan_index) + A15(phi_index,pi_fan_index);
            I7(phi_index,pi_fan_index) = Ix(phi_index,pi_fan_index) + I15(phi_index,pi_fan_index);
            m_dot_7(phi_index,pi_fan_index)      = m_dot_x(phi_index,pi_fan_index) + m_dot_0(phi_index,pi_fan_index);
            ht_initial_7(phi_index,pi_fan_index) = (ht_0(phi_index,pi_fan_index) + ht_H2_turbine(phi_index,pi_fan_index) + ht_O2_turbine(phi_index,pi_fan_index)) / m_dot_7(phi_index,pi_fan_index);

            % Find the mass fractions of oxygen, hydrogen, and nitrogen in
            % the combined stream
            yH2_7(phi_index,pi_fan_index) = m_dot_H2_turbine(phi_index,pi_fan_index) / m_dot_7(phi_index,pi_fan_index);
            yO2_7(phi_index,pi_fan_index) = (0.768 * m_dot_0(phi_index,pi_fan_index) + m_dot_O2_turbine(phi_index,pi_fan_index)) / m_dot_7(phi_index,pi_fan_index);
            yN2_7(phi_index,pi_fan_index) = (0.232 * m_dot_0(phi_index,pi_fan_index)) / m_dot_7(phi_index,pi_fan_index);
            MWT_7(phi_index,pi_fan_index) = MWT_yHyOyN(yH2_7(phi_index,pi_fan_index), yO2_7(phi_index,pi_fan_index), yN2_7(phi_index,pi_fan_index));
            R7(phi_index,pi_fan_index)    = 5.97994*(8314.4598/MWT_7(phi_index,pi_fan_index));
            Cp7(phi_index,pi_fan_index)   = R7(phi_index,pi_fan_index) * (gamma / (gamma-1));

            % Find TtBurned
            Tt7(phi_index,pi_fan_index) = Ttbrn_yHyOyNhi(yH2_7(phi_index,pi_fan_index), yO2_7(phi_index,pi_fan_index), yN2_7(phi_index,pi_fan_index), ht_initial_7(phi_index,pi_fan_index), gamma);
            ht7(phi_index,pi_fan_index) = Tt7(phi_index,pi_fan_index) * Cp7(phi_index,pi_fan_index);

            % Find Mfi parameter and Mach
            Mfi_7(phi_index,pi_fan_index) = ht7(phi_index,pi_fan_index) * (m_dot_7(phi_index,pi_fan_index) / I7(phi_index,pi_fan_index))^2;
            M7(phi_index,pi_fan_index) = MbxMfi(Mfi_7(phi_index,pi_fan_index), gamma);

            % Find Pt7 using Mft
            Pt7(phi_index,pi_fan_index) = m_dot_7(phi_index,pi_fan_index) * sqrt(R7(phi_index,pi_fan_index) * Tt7(phi_index,pi_fan_index)) / (MftxM(M7(phi_index,pi_fan_index), gamma) * A7(phi_index,pi_fan_index));


            % Station 9

            % Assume these are true for isentropic, ideal expansion
            m_dot_9(phi_index,pi_fan_index) = m_dot_7(phi_index,pi_fan_index);
            Tt9(phi_index,pi_fan_index) = Tt7(phi_index,pi_fan_index);
            P9  = P0;
            Pt9_ideal(phi_index,pi_fan_index) = Pt7(phi_index,pi_fan_index);
            R9_ideal(phi_index,pi_fan_index)  = R7(phi_index,pi_fan_index);

            % Find ideal station 9 properties
            M9_ideal(phi_index,pi_fan_index) = MxPri(P9/Pt9_ideal(phi_index,pi_fan_index), gamma);
            T9_ideal(phi_index,pi_fan_index) = TrixM(M9_ideal(phi_index,pi_fan_index), gamma) * Tt9(phi_index,pi_fan_index);
            V9_ideal(phi_index,pi_fan_index) = M9_ideal(phi_index,pi_fan_index) * sqrt(gamma * R9_ideal(phi_index,pi_fan_index) * T9_ideal(phi_index,pi_fan_index));

            % Find actual gross thrust
            Fg_ideal(phi_index,pi_fan_index) = m_dot_9(phi_index,pi_fan_index) * V9_ideal(phi_index,pi_fan_index);
            Fg(phi_index,pi_fan_index) = Cfg * Fg_ideal(phi_index,pi_fan_index);

            % Find ram drag and drag increments
            Dram(phi_index,pi_fan_index) = m_dot_0(phi_index,pi_fan_index) * u0;

            Acowl(phi_index,pi_fan_index) = 0.1 * A0(phi_index,pi_fan_index);
            Dcowl(phi_index,pi_fan_index) = q0 * Cd * Acowl(phi_index,pi_fan_index);

            % Calc net propulsive force
            Fnp(phi_index,pi_fan_index) = Fg(phi_index,pi_fan_index) - Dram(phi_index,pi_fan_index) - Dcowl(phi_index,pi_fan_index);

            A0(phi_index,pi_fan_index) = A0(phi_index,pi_fan_index) + 0.01;
            counter = counter + 1;
            % end while loop
    
        end

        if counter < 10000

            % Calculate performance parameters
            SFC(phi_index,pi_fan_index) = (m_dot_turbine(phi_index,pi_fan_index) * 32.17*3600) / Fnp(phi_index,pi_fan_index);
            Isp(phi_index,pi_fan_index) = Fnp(phi_index,pi_fan_index) / (32.2 * m_dot_turbine(phi_index,pi_fan_index));
            Fs(phi_index,pi_fan_index)  = Fnp(phi_index,pi_fan_index) / m_dot_0(phi_index,pi_fan_index);

            if pi_fan == 2.0 && phi == 2 

                fprintf('\nFan: \n');
                fprintf(append('Pt12:           ', num2str(Pt12), '\n'));
                fprintf(append('Pt13:           ', num2str(Pt13(2,10)), '\n'));
                fprintf(append('A0:             ', num2str(A0(2,10)), '\n'));
                fprintf(append('m_dot_0:        ', num2str(m_dot_0(2,10)), '\n'));
                fprintf(append('power_fan:      ', num2str(power_fan(2,10)), '\n'));
                fprintf(append('tau_fan:        ', num2str(tau_fan(2,10)), '\n'));
                
                fprintf('\nStation 15: \n');
                fprintf(append('A15:           ', num2str(A15(2,10)), '\n'));
                fprintf(append('M15:           ', num2str(M15), '\n'));
                fprintf(append('T15:           ', num2str(T15(2,10)), '\n'));
                fprintf(append('Tt15:          ', num2str(Tt15(2,10)), '\n'));
                fprintf(append('P15:           ', num2str(P15(2,10)), '\n'));
                fprintf(append('Pt15:          ', num2str(Pt15(2,10)), '\n'));
                fprintf(append('v15:           ', num2str(v15(2,10)), '\n'));
                fprintf(append('I15:           ', num2str(I15(2,10)), '\n'));

                fprintf('\nTtb: \n');
                fprintf(append('Tt4:           ', num2str(Tt4(2,10)), '\n'));

                fprintf('\nTurbine: \n');
                fprintf(append('tau_turbine:      ', num2str(tau_turbine(2,10)), '\n'));
                fprintf(append('power_turbine:    ', num2str(power_turbine(2,10)), '\n'));
                fprintf(append('Ttx:              ', num2str(Ttx(2,10)), '\n'));
                fprintf(append('Tx:               ', num2str(Tx(2,10)), '\n'));
                fprintf(append('Ptx:              ', num2str(Ptx(2,10)), '\n'));
                fprintf(append('Px:               ', num2str(Px(2,10)), '\n'));
                fprintf(append('Cp:               ', num2str(Cp_c(2,10)), '\n'));
                fprintf(append('R:                ', num2str(Rc(2,10)), '\n'));
                fprintf(append('m_dot_turbine:    ', num2str(m_dot_turbine(2,10)), '\n'));
                fprintf(append('m_dot_H2_turbine: ', num2str(m_dot_H2_turbine(2,10)), '\n'));
                fprintf(append('m_dot_O2_turbine: ', num2str(m_dot_O2_turbine(2,10)), '\n'));
                fprintf(append('ht_H2_turbine:    ', num2str(ht_H2_turbine(2,10)), '\n'));
                fprintf(append('ht_O2_turbine:    ', num2str(ht_O2_turbine(2,10)), '\n'));

                fprintf('\nStation 7: \n');
                fprintf(append('m_dot_7:          ', num2str(m_dot_7(2,10)), '\n'));
                fprintf(append('Pt7:              ', num2str(Pt7(2,10)), '\n'));
                fprintf(append('Tt7:              ', num2str(Tt7(2,10)), '\n'));
                fprintf(append('M7:               ', num2str(M7(2,10)), '\n'));
                fprintf(append('ht7:              ', num2str(ht7(2,10)), '\n'));
                fprintf(append('I7:               ', num2str(I7(2,10)), '\n'));
                fprintf(append('A7:               ', num2str(A7(2,10)), '\n'));
                fprintf(append('ht_initial_7:     ', num2str(ht_initial_7(2,10)), '\n'));
                fprintf(append('Cp7:              ', num2str(Cp7(2,10)), '\n'));
                fprintf(append('R7:               ', num2str(R7(2,10)), '\n'));

                fprintf('\nStation 9: \n');
                fprintf(append('M9_ideal:              ', num2str(M9_ideal(2,10)), '\n'));
                fprintf(append('T9_ideal:              ', num2str(T9_ideal(2,10)), '\n'));
                fprintf(append('V9_ideal:              ', num2str(V9_ideal(2,10)), '\n'));

                fprintf('\nPerformance: \n');
                fprintf(append('Fnp:            ', num2str(Fnp(2,10)), '\n'));
                fprintf(append('SFC:            ', num2str(SFC(2,10)), '\n'));
                fprintf(append('Fs:             ', num2str(Fs(2,10)), '\n'));
                fprintf(append('Isp:            ', num2str(Isp(2,10)), '\n'));                

            end

        else

            fprintf('hit max iter\n');
            fprintf(append('phi:    ', num2str(phi), '\n'));
            fprintf(append('pi_fan: ', num2str(pi_fan), '\n'));

        end

        pi_fan_index = pi_fan_index + 1;
    
    end
    pi_fan_index = 1;
    phi_index = phi_index + 1;

end

%% Plotting for Cruise, Max Isp calculation

figure(1);
hold on;
plot(pi_fan_list,SFC(1,:));
plot(pi_fan_list,SFC(2,:));
plot(pi_fan_list,SFC(3,:));
xlabel('\pi_{fan} Sweep');
ylabel('TSFC');
title('TSFC Values over \pi_{fan} Sweep for Launch');
legend('\phi = 1','\phi = 2','\phi = 10');
hold off;

figure(2);
hold on;
plot(pi_fan_list,Isp(1,:));
plot(pi_fan_list,Isp(2,:));
plot(pi_fan_list,Isp(3,:));
xlabel('\pi_{fan} Sweep');
ylabel('Isp (s)');
title('Isp Values over \pi_{fan} Sweep for Launch');
legend('\phi = 1','\phi = 2','\phi = 10');
hold off;

figure(3);
hold on;
plot(pi_fan_list,Fs(1,:));
plot(pi_fan_list,Fs(2,:));
plot(pi_fan_list,Fs(3,:));
xlabel('\pi_{fan} Sweep');
ylabel('Fs (lb)');
title('Specific Thrust Values over \pi_{fan} Sweep for Launch');
legend('\phi = 1','\phi = 2','\phi = 10');
hold off;

%maximum Isp
[Isp_max_cruise,Isp_max_indices_cruise] = max(Isp,[],2);
