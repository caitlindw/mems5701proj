%General values
Tt_p = 540; %Rankine
Pt_c = 2*10^5; %psf
Cfg = 0.96;
gamma = 1.4;
R0 = 49710/(28.97*0.06852177); %gas constant / mixture avg molecular weight (slug/Kmol)
Cp0 = ((gamma)/(gamma-1))*R0; %lbf/slug-R

%TO DO
% compute R, C_p for all species


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
P0_cruise   = 1.049 * 144;  % lbf / sq_ft

F_np_min_cruise = 5000; %lbf

V0_cruise = M0_cruise*sqrt(gamma*R0*T0_cruise);

%cowl drag cruise
C_D_cruise = 0.4;
pi_inlet_cruise = 0.669;
A0_cruise = [1*10^(-4):0.1:15];
q_cruise = 0.5*rho0_cruise*V0_cruise^2;
A_cowl_cruise = 0.10*A0_cruise;
D_cowl_cruise = q_cruise*C_D_cruise*A_cowl_cruise;

%section 15 & X to 7 
M15_cruise = 0.2;

%section 7 -> 9 isentropic flow, mixer

%select a0 using 2 nested for-loops
%equivalence = outer loop
    

    %area = inner loop
        %compute m_dot_0 = m_dot_15 from A0
        %m_dot_0 = V0*rho0*A0
        
        %find Pt15 from flight condition + recovery
            %can get Pt0
            %get Pt15 = Pt2 = Pt0*pi_inlet
            %M15 = 0.2
            %could find Tt15, then A15
            
            %innermost loop: guess m_dotc = 200lb/s
                %stop condition = F_N given -> make a while loop
                %for selected ph and beta (bypass ratio) = 0, compute mass flows of H, O

                    %Station C
                    %Mk_O2 = 32
                    %Mk_H2 = 2
                    %phi = (m_dotH2/m_dotO2)/(mH2_st/mO2_st)
                    %mH2_st = 2*2;
                    %mO2_st = 1*32;
                    %syms m_dotH2 m_dotO2
                    %eq1 = phi == (m_dotH2/m_dotO2)/(mH2_st/mO2_st);
                    %eq2 = m_dot_c = m_dotH2 + m_dotO2;
                    %[m_dotH2_sol,m_dotO2_sol] = solve([eq1,eq2],[m_dotH2,m_dotO2])
                    
                    %ht_H2 = Cp_H2*Tt_p*m_dot_H2_sol;
                    %ht_O2 = Cp_O2*Tt_p*m_dot_O2_sol;
                    %yN_p = 0;
                    %yO_p = 32/34;
                    %yH_p = 2/34;
                    %MWT_c = MWT_yHyOyN(yH_p, yO_p, 0);
                    %Rc = 49710/MWT_t;                    

                    %Station x
                    %m_dot_x = m_dot_c;
                    %h_int_p = (ht_H2 + ht_O2)/m_dot_x;
                    %Ttx = Ttbrn_yHyOyNhi(yH_p, yO_p, yN_p, h_int_p,
                    %gamma);
                    
                    %Pt_x = Pt_c;
                    %Pr_15 = PrixM(M15, gamma); %P/Pt
                    %P15 = Pr_15*Pt15;
                    %Px = P15;
                    %Mx = MxPri(Px/Pt_x, gamma);
                    %Mft_x = MftxM(Mx, gamma);
                    %Tx = TrixM(Mx, gamma)*Ttx;
                    %Rx = Rc;
                    %Ax = m_dot_x*sqrt(Rx*Ttx)/(Mft_x*Pt_x);

                    %station 7
                    %m_dot_7 = m_dot_H2_sol + m_dot_O2_sol + m_dot_15;
                    %m_dot_2 = m_dot_15;
                    %h_2 = m_dot_2*Cp0*(T0/TrixM(M0, gamma));
                    %ht_7 = (1/m_dot_7)*(h_2 + h_H2 + h_O2);
                    %Vx = Mx*sqrt(gamma*Rx*Tx);
                    %Ix = m_dot_x*Vx+Px*Ax;
                    %V15 = M15*sqrt(gamma*R0*(TrixM(M15, gamma))*Tt15)
                    %I15 = m_dot_15*V15+P15*A15;
                    
                    
                    %mfi_7 = ht_7*(m_dot_7/(I15+Ix))^2;
                    %M7 = MbxMfi(mfi_7, gamma);

                    %m_dot_O2_f = m_dot_0*0.232;
                    %m_dot_N2 = m_dot_0*0.768;
                    %yH2_7 = m_dot_H2_sol/m_dot_7;
                    %yN2_7 = m_dot_N2/m_dot_7;
                    %yO2_7 = (m_dot_O2_f+m_dot_O2_sol)/m_dot_7;
                    %MWT_7 = MWT_yHyOyN(yH2_7, yO2_7, yN2_7);
                    %R7 = 49710/MWT_7;
                    %Cp_7 = gamma/(gamma-1)*R7;
                    %Tt7 = ht_7/Cp_7;

                    %Pt7 = (1/PrixM(M7, gamma))*P7

                    %station 9
                    %m_dot_9 = m_dot_7;
                    %Tt9i = Tt7;
                    %P9i = P0;
                    %Pt9i = Pt7;
                    %M9i = MxPri(P9i/Pt9i, gamma);
                    %R9i = R7;
                    %T9i = TrixM(M9i, gamma)*Tt9i;
                    %V9i = M9i*sqrt(gamma*R9i*T9i);

                    %Fgi = m_dot_7*V9i;
                    %Fg = Cfg*Fgi;
                    %Dram = m_dot_0*V0;
                    %q = 0.5*rho0*V0^2;
                    %D_cowl = q*C_D*A_cowl;
                    %F_N = Fg - Dram - D_cowl;

                    %check if F_N requirement for cruise/launch is met
                    

                    
    


