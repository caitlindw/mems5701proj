%General values
Tt_p = 540; %Rankine
Pt_chamber = 2*10^5; %psf
Cfg = 0.96;
gamma = 1.4;
R0 = 49710/(28.97*0.06852177); %gas constant / mixture avg molecular weight (slug/Kmol)
C_P0 = ((gamma)/(gamma-1))*R0; %lbf/slug-R

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

%h_init_7_launch = (1/m_dot_7_launch)*(m_dot_H2_launch*C_p_H2_launch*Tf_H2_feed_launch+m_dot_O2_launch*C_p_O2_launch*Tf_O2_feed_launch+m_dot_2_launch*C_p_air_launch*Tt2_launch);

%section 7 -> 9 isentropic flow
%F = m_dot*V_exhaust;

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
        %compute m_dot0 = m_dot15 from A0
        %mdot_0 = V0*rho0*A0
        
        %find Pt15 from flight condition + recovery
            %can get Pt0
            %get Pt15 = Pt2 = Pt0*pi_inlet
            %M15 = 0.2
            %could find Tt15, then A15
            
            %innermost loop: guess m_dotc = 200lb/s
                %stop condition = F_N given -> make a while loop
                %for selected ph and beta (bypass ratio) = 0, compute mass flows of H, O
                    %Mk_O2 = 32
                    %Mk_H2 = 4
                    %phi = (m_dotH2/m_dotO2)/(mH2_st/mO2_st)
                    %mH2_st = 2*2;
                    %mO2_st = 1*32;
                    %syms m_dotH2 m_dotO2
                    %eq1 = phi == (m_dotH2/m_dotO2)/(mH2_st/mO2_st);
                    %eq2 = m_dotc = m_dotH2 + m_dotO2;
                    %[m_dotH2_sol,m_dotO2_sol] = solve([eq1,eq2],[m_dotH2,m_dotO2])
                    
                    %MAYBE
                    %find At using C* = sqrt(R*Tt)/mft_M=1 and m_dot =
                    %a*P_c*Ab*rho_s and At = m_dot*C*/Pt

                    
    


