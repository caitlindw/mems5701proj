%General values
Tt_p = 540; %Rankine
Pt_c = 2*10^5; %psf
Pt_x = Pt_c;
Cfg = 0.96;
gamma = 1.4;
R0 = 49710/(28.97*0.06852177); %gas constant / mixture avg molecular weight (slug/Kmol)
Cp0 = ((gamma)/(gamma-1))*R0; %lbf/slug-R
Cp_H2 = 82.65; %lbf/slug-R
Cp_O2 = 5.2925; %lbf/slug-R
A0_list = 1*10^-4:0.1:15.0001;

%finding mass fractions in chamber to find Rc
yN_p = 0;
yO_p = 32/34;
yH_p = 2/34;
MWT_c = MWT_yHyOyN(yH_p, yO_p, 0);
Rc = 49710/MWT_c;  
Rx = Rc;


%% Launch condition: M=0.8, 30,000 ft altitude, min F_np = 4000 lbf
% Freestream Conditions
M0_launch   = 0.8;
rho0_launch = 8.91e-4;      % slug / cu_ft - Atmosphere Table
T0_launch   = 411.84;       % R - Atmosphere Table
P0_launch   = 4.373 * 144;  % lbf / sq_ft

F_np_min_launch = 4000; %lbf
V0_launch = M0_launch*sqrt(gamma*R0*T0_launch);


%cowl drag, launch
C_D_launch = 0.2;
pi_inlet_launch = 0.96;
A0_launch = [1*10^(-4):0.1:15];
q_launch = 0.5*rho0_launch*V0_launch^2;
A_cowl_launch = 0.10*A0_launch;
D_cowl_launch = q_launch*C_D_launch*A_cowl_launch;

%section 15 & X to 7 , mixer
M15_launch = 0.1;

Tt2_launch = 1/(TrixM(M0_launch, gamma))*T0_launch;

%% Cruise condition: M=4, 60,000 ft altitude, min F_np = 5000 lbf
% Freestream Conditions
M0_cruise   = 4.0;
rho0_cruise = 2.26e-4;      % slug / cu_ft - Atmosphere Table
T0_cruise   = 389.97;       % R - Atmosphere Table
Tt0_cruise = (1/TrixM(M0_cruise, gamma))*T0_cruise;
P0_cruise   = 1.049 * 144;  % lbf / sq_ft
Pt0_cruise = (1/PrixM(M0_cruise, gamma))*P0_cruise;

F_np_min_cruise = 5000; %lbf

V0_cruise = M0_cruise*sqrt(gamma*R0*T0_cruise);

%cowl drag cruise
C_D_cruise = 0.4;
pi_inlet_cruise = 0.669;
q_cruise = 0.5*rho0_cruise*V0_cruise^2;

%section 15 & X to 7 
M15_cruise = 0.2;

%intialize matrices for claculations
m_dot_0_cruise = zeros(4,151);
m_dot_15_cruise = zeros(4,151);
Pt15_cruise = zeros(4,151);
Tt15_cruise = zeros(4,151);
A15_cruise = zeros(4,151);
m_dot_H2_sol_cruise = zeros(4,151);
m_dot_O2_sol_cruise = zeros(4,151);
ht_H2_cruise = zeros(4,151);
ht_O2_cruise = zeros(4,151);
h_int_p_cruise = zeros(4,151);
Ttx_cruise = zeros(4,151);
P15_cruise = zeros(4,151);
Px_cruise = zeros(4,151);
Mx_cruise = zeros(4,151);
Mft_x_cruise = zeros(4,151);
Tx_cruise = zeros(4,151);
Ax_cruise = zeros(4,151);
m_dot_7_cruise = zeros(4,151);
m_dot_2_cruise = zeros(4,151);
ht_2_cruise = zeros(4,151);
ht_7_cruise = zeros(4,151);
Vx_cruise = zeros(4,151);
Ix_cruise = zeros(4,151);
V15_cruise = zeros(4,151);
I15_cruise = zeros(4,151);
mfi_7_cruise = zeros(4,151);
M7_cruise = zeros(4,151);
m_dot_O2_f_cruise = zeros(4,151);
m_dot_N2_cruise = zeros(4,151);
yH2_7_cruise = zeros(4,151);
yN2_7_cruise = zeros(4,151);
yO2_7_cruise = zeros(4,151);
MWT_7_cruise = zeros(4,151);
R7_cruise = zeros(4,151);
Cp_7_cruise = zeros(4,151);
Tt7_cruise = zeros(4,151);
A7_cruise = zeros(4,151);
Pt7_cruise = zeros(4,151);
 m_dot_9_cruise = zeros(4,151);
Tt9i_cruise = zeros(4,151);
Pt9i_cruise = zeros(4,151);
M9i_cruise  = zeros(4,151);
R9i_cruise = zeros(4,151);
T9i_cruise = zeros(4,151);
V9i_cruise = zeros(4,151);
Fgi_cruise  = zeros(4,151);
Fg_cruise  = zeros(4,151);
Dram_cruise  = zeros(4,151);
F_N_cruise  = zeros(4,151);
m_dot_c_cruise = zeros(4,151);
Isp_cruise = zeros(4,151);
TSFC_cruise = zeros(4,151);
Fs_cruise = zeros(4,151);
phi_i = 1;
A0_i = 1;

%equivalence = outer loop
for phi = [1 2 10 1000]
    
    %area = inner loop
    for A0 = 1*10^-4:0.1:15.0001
        %compute m_dot_0 = m_dot_15 from A0
        m_dot_0_cruise(phi_i,A0_i) = Pt0_cruise*A0*(MftxM(M0_cruise, gamma))/sqrt(R0*Tt0_cruise);
        m_dot_15_cruise(phi_i,A0_i) = m_dot_0_cruise(phi_i,A0_i);
        
        %find Pt15 from flight condition + recovery
        %can get Pt0
        %get Pt15 = Pt2 = Pt0*pi_inlet
        Pt15_cruise(phi_i,A0_i) = Pt0_cruise*pi_inlet_cruise;
        Tt15_cruise(phi_i,A0_i) = Tt0_cruise;
        A15_cruise(phi_i,A0_i) = m_dot_15_cruise(phi_i,A0_i)*sqrt(R0*Tt15_cruise(phi_i,A0_i))/(Pt15_cruise(phi_i,A0_i)*(MftxM(M15_cruise, gamma)));
        
        %innermost loop: guess m_dot_c until the cruise condition is met
        m_dot_c = 0.9;
        counter = 0;
        %stop condition = F_N given -> make a while loop
        while F_N_cruise(phi_i,A0_i) < F_np_min_cruise-0.1 || counter >= 10000
            m_dot_c = m_dot_c+0.1;
            %for selected ph and beta (bypass ratio) = 0, compute mass flows of H, O

            %Station C
            %Mk_O2 = 32;
            %Mk_H2 = 2;

            %finding mass flow of H2, O2
            %mH2_st = 2
            %mO2_st = 32
            syms m_dotH2 m_dotO2
            eq1 = phi == (m_dotH2/m_dotO2)/(2/32);
            eq2 = m_dot_c == m_dotH2 + m_dotO2;
            [m_dot_H2_sol_cruise(phi_i,A0_i),m_dot_O2_sol_cruise(phi_i,A0_i)] = solve([eq1,eq2],[m_dotH2,m_dotO2]);
                             

            %Station x: Ttx
            m_dot_x = m_dot_c;
            ht_H2_cruise(phi_i,A0_i) = Cp_H2*Tt_p*m_dot_H2_sol_cruise(phi_i,A0_i);
            ht_O2_cruise(phi_i,A0_i) = Cp_O2*Tt_p*m_dot_O2_sol_cruise(phi_i,A0_i);
            h_int_p_cruise(phi_i,A0_i) = (ht_H2_cruise(phi_i,A0_i) + ht_O2_cruise(phi_i,A0_i))/m_dot_x;
            Ttx_cruise(phi_i,A0_i) = Ttbrn_yHyOyNhi(yH_p, yO_p, yN_p, h_int_p_cruise(phi_i,A0_i),gamma);
            
            %Station 5
            P15_cruise(phi_i,A0_i) = PrixM(M15_cruise, gamma)*Pt15_cruise(phi_i,A0_i);
            Px_cruise(phi_i,A0_i) = P15_cruise(phi_i,A0_i);
            Mx_cruise(phi_i,A0_i) = MxPri(Px_cruise(phi_i,A0_i)/Pt_x, gamma);
            Mft_x_cruise(phi_i,A0_i) = MftxM(Mx_cruise(phi_i,A0_i), gamma);
            Tx_cruise(phi_i,A0_i) = TrixM(Mx_cruise(phi_i,A0_i), gamma)*Ttx_cruise(phi_i,A0_i);
            Ax_cruise(phi_i,A0_i) = m_dot_x*sqrt(Rx*Ttx_cruise(phi_i,A0_i))/(Mft_x_cruise(phi_i,A0_i)*Pt_x);

            %station 7
            m_dot_7_cruise(phi_i,A0_i) = m_dot_H2_sol_cruise(phi_i,A0_i) + m_dot_O2_sol_cruise(phi_i,A0_i) + m_dot_15_cruise(phi_i,A0_i);
            m_dot_2_cruise(phi_i,A0_i) = m_dot_15_cruise(phi_i,A0_i);
            ht_2_cruise(phi_i,A0_i) = m_dot_2_cruise(phi_i,A0_i)*Cp0*(T0_cruise/TrixM(M0_cruise, gamma));
            ht_7_cruise(phi_i,A0_i) = (1/m_dot_7_cruise(phi_i,A0_i))*(ht_2_cruise(phi_i,A0_i) + ht_H2_cruise(phi_i,A0_i) + ht_O2_cruise(phi_i,A0_i));
            Vx_cruise(phi_i,A0_i) = Mx_cruise(phi_i,A0_i)*sqrt(gamma*Rx*Tx_cruise(phi_i,A0_i));
            Ix_cruise(phi_i,A0_i) = m_dot_x*Vx_cruise(phi_i,A0_i)+Px_cruise(phi_i,A0_i)*Ax_cruise(phi_i,A0_i);
            V15_cruise(phi_i,A0_i) = M15_cruise*sqrt(gamma*R0*(TrixM(M15_cruise, gamma))*Tt15_cruise(phi_i,A0_i));
            I15_cruise(phi_i,A0_i) = m_dot_15_cruise(phi_i,A0_i)*V15_cruise(phi_i,A0_i)+P15_cruise(phi_i,A0_i)*A15_cruise(phi_i,A0_i);            
            
            %Station 7 Mach and Other Properties
            mfi_7_cruise(phi_i,A0_i) = ht_7_cruise(phi_i,A0_i)*(m_dot_7_cruise(phi_i,A0_i)/(I15_cruise(phi_i,A0_i)+Ix_cruise(phi_i,A0_i)))^2;
            M7_cruise(phi_i,A0_i) = MbxMfi(mfi_7_cruise(phi_i,A0_i), gamma);

            m_dot_O2_f_cruise(phi_i,A0_i) = m_dot_0_cruise(phi_i,A0_i)*0.232;
            m_dot_N2_cruise(phi_i,A0_i) = m_dot_0_cruise(phi_i,A0_i)*0.768;
            yH2_7_cruise(phi_i,A0_i) = m_dot_H2_sol_cruise(phi_i,A0_i)/m_dot_7_cruise(phi_i,A0_i);
            yN2_7_cruise(phi_i,A0_i) = m_dot_N2_cruise(phi_i,A0_i)/m_dot_7_cruise(phi_i,A0_i);
            yO2_7_cruise(phi_i,A0_i) = (m_dot_O2_f_cruise(phi_i,A0_i)+m_dot_O2_sol_cruise(phi_i,A0_i))/m_dot_7_cruise(phi_i,A0_i);
            MWT_7_cruise(phi_i,A0_i) = MWT_yHyOyN(yH2_7_cruise(phi_i,A0_i), yO2_7_cruise(phi_i,A0_i), yN2_7_cruise(phi_i,A0_i));
            R7_cruise(phi_i,A0_i) = 49710/MWT_7_cruise(phi_i,A0_i);
            Cp_7_cruise(phi_i,A0_i) = gamma/(gamma-1)*R7_cruise(phi_i,A0_i);
            Tt7_cruise(phi_i,A0_i) = ht_7_cruise(phi_i,A0_i)/Cp_7_cruise(phi_i,A0_i);
            A7_cruise(phi_i,A0_i) = A15_cruise(phi_i,A0_i)+Ax_cruise(phi_i,A0_i);
            Pt7_cruise(phi_i,A0_i) = m_dot_7_cruise(phi_i,A0_i)*sqrt(R7_cruise(phi_i,A0_i)*Tt7_cruise(phi_i,A0_i))/(MftxM(M7_cruise(phi_i,A0_i), gamma)*A7_cruise(phi_i,A0_i));

            %station 9
            m_dot_9_cruise(phi_i,A0_i) = m_dot_7_cruise(phi_i,A0_i);
            Tt9i_cruise(phi_i,A0_i) = Tt7_cruise(phi_i,A0_i);
            P9i_cruise = P0_cruise;
            Pt9i_cruise(phi_i,A0_i) = Pt7_cruise(phi_i,A0_i);
            M9i_cruise(phi_i,A0_i) = MxPri(P9i_cruise/Pt9i_cruise(phi_i,A0_i), gamma);
            R9i_cruise(phi_i,A0_i) = R7_cruise(phi_i,A0_i);
            T9i_cruise(phi_i,A0_i) = TrixM(M9i_cruise(phi_i,A0_i), gamma)*Tt9i_cruise(phi_i,A0_i);
            V9i_cruise(phi_i,A0_i) = M9i_cruise(phi_i,A0_i)*sqrt(gamma*R9i_cruise(phi_i,A0_i)*T9i_cruise(phi_i,A0_i));

            Fgi_cruise(phi_i,A0_i) = m_dot_7_cruise(phi_i,A0_i)*V9i_cruise(phi_i,A0_i);
            Fg_cruise(phi_i,A0_i) = Cfg*Fgi_cruise(phi_i,A0_i);
            Dram_cruise(phi_i,A0_i) = m_dot_0_cruise(phi_i,A0_i)*V0_cruise;
            A_cowl_cruise = 0.10*A0;
            D_cowl_cruise = q_cruise*C_D_cruise*A_cowl_cruise;
            F_N_cruise(phi_i,A0_i) = Fg_cruise(phi_i,A0_i) - Dram_cruise(phi_i,A0_i) - D_cowl_cruise;

            %check if F_N requirement for cruise/launch is met
            counter = counter + 1;
        end
        if counter >= 10000
           fprintf('Hit max iterations :(');
        else
            m_dot_c_cruise(phi_i,A0_i) = m_dot_c;
            TSFC_cruise(phi_i,A0_i) = (m_dot_c*3600)/F_N_cruise(phi_i,A0_i);
            Isp_cruise(phi_i,A0_i) = F_N_cruise(phi_i,A0_i)/(32.2*m_dot_c);
            Fs_cruise(phi_i,A0_i) = F_N_cruise(phi_i,A0_i)/(m_dot_0_cruise(phi_i,A0_i)); %F/m_dot_0
            %fprintf(['phi = ',num2str(phi),' A0 = ',num2str(A0),' m_dot_c = ',num2str(m_dot_c_cruise(phi_i,A0_i)) '\n']);
        end

        A0_i = A0_i+1;      
    end
    A0_i = 1;
    phi_i = phi_i+1;
end

%% Plotting for Cruise, Max Isp calculation

figure(1);
hold on;
plot(A0_list,TSFC_cruise(1,:));
plot(A0_list,TSFC_cruise(2,:));
plot(A0_list,TSFC_cruise(3,:));
plot(A0_list,TSFC_cruise(4,:));
xlabel('A0 Sweep (ft^2)');
ylabel('TSFC');
title('TSFC Values over A0 Sweep for Cruise');
legend('\phi = 1','\phi = 2','\phi = 10','\phi = 1000');
hold off;

figure(2);
hold on;
plot(A0_list,Isp_cruise(1,:));
plot(A0_list,Isp_cruise(2,:));
plot(A0_list,Isp_cruise(3,:));
plot(A0_list,Isp_cruise(4,:));
xlabel('A0 Sweep (ft^2)');
ylabel('Isp (s)');
title('Isp Values over A0 Sweep for Cruise');
legend('\phi = 1','\phi = 2','\phi = 10','\phi = 1000');
hold off;

figure(3);
hold on;
plot(A0_list,Fs_cruise(1,:));
plot(A0_list,Fs_cruise(2,:));
plot(A0_list,Fs_cruise(3,:));
plot(A0_list,Fs_cruise(4,:));
xlabel('A0 Sweep (ft^2)');
ylabel('Fs (lb)');
title('Specific Thrust Values over A0 Sweep for Cruise');
legend('\phi = 1','\phi = 2','\phi = 10','\phi = 1000');
hold off;

%maximum Isp
[Isp_max,Isp_max_indices] = max(Isp_cruise,[],2);