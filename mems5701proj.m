%General values
Tt_p = 540; %Rankin
Pt_chamber = 2*10^5; %psf
Cfg = 0.96;
gamma = 1.4;
R0_air = 49710/(28.97*0.06852177); %gas constant / mixture avg molecular weight (slug/Kmol)

%TO DO
% get standard atmosphere table values
% compute R, C_p for all species

%% Launch condition: M=0.8, 30,000 ft altitude, min F_np = 4000 lbf
M0_launch = 0.8;
rho0_launch = 0; %find using altitude
T0_launch = 0; %find using altitude
F_np_min_launch = 4000; %lbf
V0_launch = M0_launch*sqrt(gamma*R0_air*T0_launch);

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
M0_cruise = 4;
rho0_cruise = 0; %find using altitude
T0_cruise = 0; %find using altitude
F_np_min_cruise = 5000; %lbf

V0_cruise = M0_cruise*sqrt(gamma*R0_air*T0_cruise);

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

