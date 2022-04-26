%General values
close all;
Tt_p = 540; %Rankine
Pt_c = 2*10^5; %psf
Pt_x = Pt_c;
Cfg = 0.96;
gamma = 1.4;
MWT_air = MWT_yHyOyN(0, 0.232, 0.768);
R0 = 5.97994*(8314.4598/MWT_air); %gas constant / mixture avg molecular weight (slug/Kmol)
Cp0 = ((gamma)/(gamma-1))*R0; %lbf/slug-R
%Cp_H2 = (14.31*1000)*5.97994; %lbf/slug-R
Cp_H2 = 87010;
%Cp_O2 = (0.918*1000)*5.97994; %lbf/slug-R
Cp_O2 = 5438;
A0_list = 1*10^-4:0.1:15.0001;
hpr = 1300000000;
a = 0.000088;
b = 540;

Beta = 0;


%% Launch condition: M=0.8, 30,000 ft altitude, min F_np = 4000 lbf
% Freestream Conditions
M0_launch   = 0.8;
rho0_launch = 8.91e-4;      % slug / cu_ft - Atmosphere Table
T0_launch   = 411.84;       % R - Atmosphere Table
P0_launch   = 4.373 * 144;  % lbf / sq_ft
Tt0_launch = (1/TrixM(M0_launch, gamma))*T0_launch;
Pt0_launch = (1/PrixM(M0_launch, gamma))*P0_launch;

F_np_min_launch = 4000; %lbf
V0_launch = M0_launch*sqrt(gamma*R0*T0_launch);


%cowl drag, launch
C_D_launch = 0.2;
pi_inlet_launch = 0.96;
q_launch = 0.5*rho0_launch*V0_launch^2;
A_cowl_launch = 0.10*A0_list;
D_cowl_launch = q_launch*C_D_launch*A_cowl_launch;

%section 15 & X to 7 , mixer
M15_launch = 0.1;

Tt2_launch = 1/(TrixM(M0_launch, gamma))*T0_launch;

%% Function

%intialize matrices for calculations
m_dot_0_launch = zeros(4,151);
m_dot_15_launch = zeros(4,151);
Pt15_launch = zeros(4,151);
Tt15_launch = zeros(4,151);
A15_launch = zeros(4,151);
m_dot_H2_sol_launch = zeros(4,151);
m_dot_O2_sol_launch = zeros(4,151);
ht_H2_launch = zeros(4,151);
ht_O2_launch = zeros(4,151);
h_int_p_launch = zeros(4,151);
Ttx_launch = zeros(4,151);
P15_launch = zeros(4,151);
Px_launch = zeros(4,151);
Mx_launch = zeros(4,151);
Mft_x_launch = zeros(4,151);
Tx_launch = zeros(4,151);
Ax_launch = zeros(4,151);
m_dot_7_launch = zeros(4,151);
m_dot_2_launch = zeros(4,151);
ht_2_launch = zeros(4,151);
ht_int_7_launch = zeros(4,151);
ht_7_launch = zeros(4,151);
Vx_launch = zeros(4,151);
Ix_launch = zeros(4,151);
V15_launch = zeros(4,151);
I15_launch = zeros(4,151);
mfi_7_launch = zeros(4,151);
M7_launch = zeros(4,151);
m_dot_O2_f_launch = zeros(4,151);
m_dot_N2_launch = zeros(4,151);
yH2_7_launch = zeros(4,151);
yN2_7_launch = zeros(4,151);
yO2_7_launch = zeros(4,151);
MWT_7_launch = zeros(4,151);
R7_launch = zeros(4,151);
Cp_7_launch = zeros(4,151);
Tt7_launch = zeros(4,151);
A7_launch = zeros(4,151);
Pt7_launch = zeros(4,151);
I7_launch = zeros(4,151);
m_dot_9_launch = zeros(4,151);
Tt9i_launch = zeros(4,151);
Pt9i_launch = zeros(4,151);
M9i_launch  = zeros(4,151);
R9i_launch = zeros(4,151);
T9i_launch = zeros(4,151);
V9i_launch = zeros(4,151);
Fgi_launch  = zeros(4,151);
Fg_launch  = zeros(4,151);
Dram_launch  = zeros(4,151);
A_cowl_launch = zeros(4,151);
D_cowl_launch = zeros(4,151);
F_N_launch  = zeros(4,151);
m_dot_c_launch = zeros(4,151);
Isp_launch = zeros(4,151);
TSFC_launch = zeros(4,151);
Fs_launch = zeros(4,151);
phi_i = 1;
A0_i = 1;
yH_p = zeros(4,1);
yO_p = zeros(4,1);
yN_p = zeros(4,1);

%equivalence = outer loop
for phi = [1 2 10 1000]
    
    %area = inner loop
    for A0 = 1*10^-4:0.1:15.0001
        %compute m_dot_0 = m_dot_15 from A0
        m_dot_0_launch(phi_i,A0_i) = Pt0_launch*A0*(MftxM(M0_launch, gamma))/sqrt(R0*Tt0_launch);
        m_dot_15_launch(phi_i,A0_i) = m_dot_0_launch(phi_i,A0_i);
        
        %find Pt15 from flight condition + recovery
        %can get Pt0
        %get Pt15 = Pt2 = Pt0*pi_inlet
        Pt15_launch(phi_i,A0_i) = Pt0_launch*pi_inlet_launch;
        Tt15_launch(phi_i,A0_i) = Tt0_launch;
        A15_launch(phi_i,A0_i) = m_dot_15_launch(phi_i,A0_i)*sqrt(R0*Tt15_launch(phi_i,A0_i))/(Pt15_launch(phi_i,A0_i)*(MftxM(M15_launch, gamma)));
        
        %innermost loop: guess m_dot_c until the launch condition is met
        m_dot_c = 0;
        counter = 0;
        %stop condition = given F_N 
        while (F_N_launch(phi_i,A0_i) < F_np_min_launch-0.001) || (counter >= 10000)
            m_dot_c = m_dot_c+0.001;


            %finding mass fractions in chamber to find Rc
            [yHON] = MassFracs(phi,Beta);
            yH_p(phi_i) = yHON(1);
            yO_p(phi_i) = yHON(2);
            yN_p(phi_i) = yHON(3);
            m_dot_H2_sol_launch(phi_i,A0_i) = m_dot_c*yH_p(phi_i);
            m_dot_O2_sol_launch(phi_i,A0_i) = m_dot_c*yO_p(phi_i);
            MWT_c = MWT_yHyOyN(yH_p(phi_i), yO_p(phi_i), yN_p(phi_i));
            Rc = 5.97994*(8314.4598/MWT_c); % (slug/Kmol)
            Rx = Rc;
            Cp_c = ((gamma)/(gamma-1))*Rc; %lbf/slug-R                             

            %Station x: Ttx
            m_dot_x = m_dot_c;
            ht_H2_launch(phi_i,A0_i) = Cp_H2*Tt_p*m_dot_H2_sol_launch(phi_i,A0_i);
            ht_O2_launch(phi_i,A0_i) = Cp_O2*Tt_p*m_dot_O2_sol_launch(phi_i,A0_i);
            h_int_p_launch(phi_i,A0_i) = (ht_H2_launch(phi_i,A0_i) + ht_O2_launch(phi_i,A0_i))/m_dot_x;
            Ttx_launch(phi_i,A0_i) = Ttbrn_yHyOyNhi(yH_p(phi_i), yO_p(phi_i), yN_p(phi_i), h_int_p_launch(phi_i,A0_i),gamma);
            
            %Station X
            P15_launch(phi_i,A0_i) = PrixM(M15_launch, gamma)*Pt15_launch(phi_i,A0_i);
            Px_launch(phi_i,A0_i) = P15_launch(phi_i,A0_i);
            Mx_launch(phi_i,A0_i) = MxPri(Px_launch(phi_i,A0_i)/Pt_x, gamma);
            Mft_x_launch(phi_i,A0_i) = MftxM(Mx_launch(phi_i,A0_i), gamma);
            Tx_launch(phi_i,A0_i) = TrixM(Mx_launch(phi_i,A0_i), gamma)*Ttx_launch(phi_i,A0_i);
            Ax_launch(phi_i,A0_i) = m_dot_x*sqrt(Rx*Ttx_launch(phi_i,A0_i))/(Mft_x_launch(phi_i,A0_i)*Pt_x);

            %station 7
            m_dot_7_launch(phi_i,A0_i) = m_dot_H2_sol_launch(phi_i,A0_i) + m_dot_O2_sol_launch(phi_i,A0_i) + m_dot_15_launch(phi_i,A0_i);
            m_dot_2_launch(phi_i,A0_i) = m_dot_15_launch(phi_i,A0_i);
            ht_2_launch(phi_i,A0_i) = m_dot_2_launch(phi_i,A0_i)*Cp0*(T0_launch/TrixM(M0_launch, gamma));
            ht_int_7_launch(phi_i,A0_i) = (1/m_dot_7_launch(phi_i,A0_i))*(ht_2_launch(phi_i,A0_i) + ht_H2_launch(phi_i,A0_i) + ht_O2_launch(phi_i,A0_i));
            
            
            
            Vx_launch(phi_i,A0_i) = Mx_launch(phi_i,A0_i)*sqrt(gamma*Rx*Tx_launch(phi_i,A0_i));
            Ix_launch(phi_i,A0_i) = m_dot_x*Vx_launch(phi_i,A0_i)+Px_launch(phi_i,A0_i)*Ax_launch(phi_i,A0_i);
            V15_launch(phi_i,A0_i) = M15_launch*sqrt(gamma*R0*(TrixM(M15_launch, gamma))*Tt15_launch(phi_i,A0_i));
            I15_launch(phi_i,A0_i) = m_dot_15_launch(phi_i,A0_i)*V15_launch(phi_i,A0_i)+P15_launch(phi_i,A0_i)*A15_launch(phi_i,A0_i);            

            m_dot_O2_f_launch(phi_i,A0_i) = m_dot_0_launch(phi_i,A0_i)*0.232;
            m_dot_N2_launch(phi_i,A0_i) = m_dot_0_launch(phi_i,A0_i)*0.768;
            yH2_7_launch(phi_i,A0_i) = m_dot_H2_sol_launch(phi_i,A0_i)/m_dot_7_launch(phi_i,A0_i);
            yN2_7_launch(phi_i,A0_i) = m_dot_N2_launch(phi_i,A0_i)/m_dot_7_launch(phi_i,A0_i);
            yO2_7_launch(phi_i,A0_i) = (m_dot_O2_f_launch(phi_i,A0_i)+m_dot_O2_sol_launch(phi_i,A0_i))/m_dot_7_launch(phi_i,A0_i);
            MWT_7_launch(phi_i,A0_i) = MWT_yHyOyN(yH2_7_launch(phi_i,A0_i), yO2_7_launch(phi_i,A0_i), yN2_7_launch(phi_i,A0_i));
            R7_launch(phi_i,A0_i) = 5.97994*(8314.4598/MWT_7_launch(phi_i,A0_i));
            Cp_7_launch(phi_i,A0_i) = gamma/(gamma-1)*R7_launch(phi_i,A0_i);

            Tt7_launch(phi_i,A0_i) = Ttbrn_yHyOyNhi(yH2_7_launch(phi_i,A0_i), yO2_7_launch(phi_i,A0_i), yN2_7_launch(phi_i,A0_i), ht_int_7_launch(phi_i,A0_i),gamma);
            ht_7_launch(phi_i,A0_i) = Tt7_launch(phi_i,A0_i) * Cp_7_launch(phi_i,A0_i);

            %Station 7 Mach and Other Properties
            I7_launch(phi_i,A0_i) = I15_launch(phi_i,A0_i)+Ix_launch(phi_i,A0_i);
            mfi_7_launch(phi_i,A0_i) = ht_7_launch(phi_i,A0_i)*(m_dot_7_launch(phi_i,A0_i)/I7_launch(phi_i,A0_i))^2;
            M7_launch(phi_i,A0_i) = MbxMfi(mfi_7_launch(phi_i,A0_i), gamma);


            A7_launch(phi_i,A0_i) = A15_launch(phi_i,A0_i)+Ax_launch(phi_i,A0_i);
            Pt7_launch(phi_i,A0_i) = m_dot_7_launch(phi_i,A0_i)*sqrt(R7_launch(phi_i,A0_i)*Tt7_launch(phi_i,A0_i))/(MftxM(M7_launch(phi_i,A0_i), gamma)*A7_launch(phi_i,A0_i));

            %station 9
            m_dot_9_launch(phi_i,A0_i) = m_dot_7_launch(phi_i,A0_i);
            Tt9i_launch(phi_i,A0_i) = Tt7_launch(phi_i,A0_i);
            P9i_launch = P0_launch;
            Pt9i_launch(phi_i,A0_i) = Pt7_launch(phi_i,A0_i);
            M9i_launch(phi_i,A0_i) = MxPri(P9i_launch/Pt9i_launch(phi_i,A0_i), gamma);
            R9i_launch(phi_i,A0_i) = R7_launch(phi_i,A0_i);
            T9i_launch(phi_i,A0_i) = TrixM(M9i_launch(phi_i,A0_i), gamma)*Tt9i_launch(phi_i,A0_i);
            V9i_launch(phi_i,A0_i) = M9i_launch(phi_i,A0_i)*sqrt(gamma*R9i_launch(phi_i,A0_i)*T9i_launch(phi_i,A0_i));

            Fgi_launch(phi_i,A0_i) = m_dot_9_launch(phi_i,A0_i)*V9i_launch(phi_i,A0_i);
            Fg_launch(phi_i,A0_i) = Cfg*Fgi_launch(phi_i,A0_i);
            Dram_launch(phi_i,A0_i) = m_dot_0_launch(phi_i,A0_i)*V0_launch;
            A_cowl_launch(phi_i,A0_i) = 0.10*A0;
            D_cowl_launch(phi_i,A0_i) = q_launch*C_D_launch*A_cowl_launch(phi_i,A0_i);
            F_N_launch(phi_i,A0_i) = Fg_launch(phi_i,A0_i) - Dram_launch(phi_i,A0_i) - D_cowl_launch(phi_i,A0_i);

            counter = counter + 1;

        end
        if counter >= 10000
           fprintf('Hit max iterations :(');
        else
            m_dot_c_launch(phi_i,A0_i) = m_dot_c;
            TSFC_launch(phi_i,A0_i) = (m_dot_c*32.17*3600)/F_N_launch(phi_i,A0_i); %convert m_dot_c from slug/s to lbm/hr
            Isp_launch(phi_i,A0_i) = F_N_launch(phi_i,A0_i)/(32.2*m_dot_c);
            Fs_launch(phi_i,A0_i) = F_N_launch(phi_i,A0_i)/(m_dot_0_launch(phi_i,A0_i)); %F/m_dot_0
        end
        
        A0_i = A0_i+1;    
    end
    A0_i = 1;
    phi_i = phi_i+1;
end

%% Printing for case of phi=10, Ao=4
fprintf(['\n \n Freesteam: \n Alt (kft): 30' '\n M: ' num2str(M0_launch)]);
fprintf(['\n P (psf): ' num2str(P0_launch) '\n Pt (psf): ' num2str(Pt0_launch) '\n T (R): ' num2str(T0_launch) '\n Tt (R): ' num2str(Tt0_launch)]);
fprintf(['\n q (psf): ' num2str(q_launch) '\n Cd0 (ft^2): ' num2str(A0_list(41)) '\n U (ft/s): ' num2str(V0_launch) '\n md (slug/s): ' num2str(m_dot_0_launch(3,41))]);
fprintf(['\n An (ft^2): ' num2str(A_cowl_launch(3,41)) '\n Dn (lbf): ' num2str(D_cowl_launch(3,41)) '\n Dram (lbf): ' num2str(Dram_launch(3,41)) '\n ']);

fprintf(['\n \n Constants: \n gamma: 1.4' '\n MWT_O2: 32 \n MWT_N2: 28.16 \n MWT_H2: 2']);
fprintf(['\n delhf: 1300000000' '\n a: 0.000088 \n b: 540 \n']);

fprintf(['\n \n Chamber: ' '\n Tt1: ' num2str(Tt_p) '\n Pt (psf): ' num2str(Pt_c) '\n md_c (slug/s): ' num2str(m_dot_c_launch(3,41))]);
fprintf(['\n phi: 10 ' '\n beta: 0' '\n Cp_H2: ' num2str(Cp_H2) '\n Cp_O2: ' num2str(Cp_O2)]);
fprintf(['\n md_H2: ' num2str(m_dot_H2_sol_launch(3,41)) '\n md_O2: ' num2str(m_dot_O2_sol_launch(3,41)) '\n md_N2: 0']);
fprintf(['\n R (ft^2/s^2): ' num2str(Rc) '\n Cp (ft^2/s^2): ' num2str(Cp_c) '\n ht_i (ft^2/s^2): ' num2str(h_int_p_launch(3,41))]);
fprintf(['\n Ttb: ' num2str(Ttx_launch(3,41)) '\n ']);

fprintf(['\n \n Primary Exit: ' '\n Mx: ' num2str(Mx_launch(3,41)) '\n Px: ' num2str(Px_launch(3,41)) '\n Tx: ' num2str(Tx_launch(3,41)) '\n Ux: ' num2str(Vx_launch(3,41))]);
fprintf(['\n Ae/A*: ' num2str(ArixM(Mx_launch(3,41), gamma)) '\n Ax: ' num2str(Ax_launch(3,41)) '\n I (lbf): ' num2str(Ix_launch(3,41))]);

fprintf(['\n \n Station 15: ' '\n M: ' num2str(M15_launch) '\n Pt: ' num2str(Pt15_launch(3,41)) '\n P: ' num2str(P15_launch(3,41)) '\n Tt: ' num2str(Tt15_launch(3,41))]);
fprintf(['\n T: ' num2str(Tt15_launch(3,41)*TrixM(M15_launch, gamma)) '\n md: ' num2str(m_dot_15_launch(3,41)) '\n R: ' num2str(R0) '\n Cp: ' num2str(Cp0)]);
fprintf(['\n A: ' num2str(A15_launch(3,41)) '\n V: ' num2str(V15_launch(3,41)) '\n I: ' num2str(I15_launch(3,41))]);

fprintf(['\n \n Station 7: ' '\n md: ' num2str(m_dot_7_launch(3,41)) '\n md_H2: ' num2str(m_dot_H2_sol_launch(3,41)) '\n md_O2: ' num2str(m_dot_O2_f_launch(3,41)+m_dot_O2_sol_launch(3,41)) '\n md_Ns: ' num2str(m_dot_N2_launch(3,41))]);
fprintf(['\n R: ' num2str(R7_launch(3,41)) '\n Cp: ' num2str(Cp_7_launch(3,41)) '\n ht: ' num2str(ht_7_launch(3,41)) '\n Tt: ' num2str(Tt7_launch(3,41)) '\n I: ' num2str(I7_launch(3,41))]);
fprintf(['\n A: ' num2str(A7_launch(3,41)) '\n M: ' num2str(M7_launch(3,41)) '\n Pt: ' num2str(Pt7_launch(3,41))]);

fprintf(['\n \n Nozzle: ' '\n R: ' num2str(R9i_launch(3,41)) '\n Tt: ' num2str(Tt9i_launch(3,41)) '\n T: ' num2str(T9i_launch(3,41)) '\n Pt: ' num2str(Pt9i_launch(3,41)) '\n P: ' num2str(P9i_launch)]);
fprintf(['\n M: ' num2str(M9i_launch(3,41)) '\n U: ' num2str(V9i_launch(3,41)) '\n F: ' num2str(Fg_launch(3,41))]);

fprintf(['\n \n Peformance: ' '\n F_N: ' num2str(F_N_launch(3,41)) '\n TSFC_eng: ' num2str(TSFC_launch(3,41)) '\n Isp: ' num2str(Isp_launch(3,41))]);

%% Plotting for launch, Max Isp calculation

figure(1);
hold on;
plot(A0_list,TSFC_launch(1,:));
plot(A0_list,TSFC_launch(2,:));
plot(A0_list,TSFC_launch(3,:));
plot(A0_list,TSFC_launch(4,:));
xlabel('A0 Sweep (ft^2)');
ylabel('TSFC');
title('TSFC Values over A0 Sweep for launch');
legend('\phi = 1','\phi = 2','\phi = 10','\phi = 1000');
hold off;

figure(2);
hold on;
plot(A0_list,Isp_launch(1,:));
plot(A0_list,Isp_launch(2,:));
plot(A0_list,Isp_launch(3,:));
plot(A0_list,Isp_launch(4,:));
xlabel('A0 Sweep (ft^2)');
ylabel('Isp (s)');
title('Isp Values over A0 Sweep for launch');
legend('\phi = 1','\phi = 2','\phi = 10','\phi = 1000');
hold off;

figure(3);
hold on;
plot(A0_list,Fs_launch(1,:));
plot(A0_list,Fs_launch(2,:));
plot(A0_list,Fs_launch(3,:));
plot(A0_list,Fs_launch(4,:));
xlabel('A0 Sweep (ft^2)');
ylabel('Fs (lb)');
title('Specific Thrust Values over A0 Sweep for launch');
legend('\phi = 1','\phi = 2','\phi = 10','\phi = 1000');
hold off;

figure(4);
hold on;
plot(A0_list(2:151),Fs_launch(1,2:151));
plot(A0_list(2:151),Fs_launch(2,2:151));
plot(A0_list(2:151),Fs_launch(3,2:151));
plot(A0_list(2:151),Fs_launch(4,2:151));
xlabel('A0 Sweep (ft^2)');
ylabel('Fs (lb)');
title('Specific Thrust Values over A0 Sweep for launch without Pure Rocket');
legend('\phi = 1','\phi = 2','\phi = 10','\phi = 1000');
hold off;

%maximum Isp
[Isp_max_launch,Isp_max_indices_launch] = max(Isp_launch,[],2);


