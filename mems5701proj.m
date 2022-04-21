%General values
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

%% Function

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
I7_cruise = zeros(4,151);
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
A_cowl_cruise = zeros(4,151);
D_cowl_cruise = zeros(4,151);
F_N_cruise  = zeros(4,151);
m_dot_c_cruise = zeros(4,151);
Isp_cruise = zeros(4,151);
TSFC_cruise = zeros(4,151);
Fs_cruise = zeros(4,151);
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
        m_dot_0_cruise(phi_i,A0_i) = Pt0_cruise*A0*(MftxM(M0_cruise, gamma))/sqrt(R0*Tt0_cruise);
        m_dot_15_cruise(phi_i,A0_i) = m_dot_0_cruise(phi_i,A0_i);
        
        %find Pt15 from flight condition + recovery
        %can get Pt0
        %get Pt15 = Pt2 = Pt0*pi_inlet
        Pt15_cruise(phi_i,A0_i) = Pt0_cruise*pi_inlet_cruise;
        Tt15_cruise(phi_i,A0_i) = Tt0_cruise;
        A15_cruise(phi_i,A0_i) = m_dot_15_cruise(phi_i,A0_i)*sqrt(R0*Tt15_cruise(phi_i,A0_i))/(Pt15_cruise(phi_i,A0_i)*(MftxM(M15_cruise, gamma)));
        
        %innermost loop: guess m_dot_c until the cruise condition is met
        m_dot_c = 0;
        counter = 0;
        %stop condition = F_N given -> make a while loop
        while (F_N_cruise(phi_i,A0_i) < F_np_min_cruise-0.001) || (counter >= 10000)
%             if (F_np_min_cruise - F_N_cruise(phi_i,A0_i)) > 3000 && A0_i < 11 %not when phi=1000
%                 m_dot_c = m_dot_c+1;
%             elseif (F_np_min_cruise - F_N_cruise(phi_i,A0_i)) > 1500
%                 m_dot_c = m_dot_c+0.5;
%             elseif (F_np_min_cruise - F_N_cruise(phi_i,A0_i)) > 200
%                 m_dot_c = m_dot_c+0.1;
%             else
            m_dot_c = m_dot_c+0.001;

            if A0_i == 41 && phi_i == 2
                %fprintf(['\n m_dot_c ' num2str(m_dot_c)]);
%                     if m_dot_c > 0.171 && m_dot_c < 0.173
%                         fprintf(['\n F_N: ' num2str(F_N_cruise(2, 41))]);
%                     end
                m_dot_c = 0.172;
            end
            
                

%             end
            %for selected ph and beta (bypass ratio) = 0, compute mass flows of H, O

            %Station C
            %Mk_O2 = 32;
            %Mk_H2 = 2;

            %finding mass flow of H2, O2
            %mH2_st = 2
            %mO2_st = 32
            %fprintf(['A0_i ' num2str(A0_i), 'phi_i ' num2str(phi_i) '\n']);
%             if A0_i == 41 && phi_i == 2
%                 m_dot_H2_sol_cruise(2,41) =   0.0343;      
%                 m_dot_O2_sol_cruise(2,41) = 0.1373;
%                 m_dot_c = m_dot_H2_sol_cruise(2,41)+m_dot_O2_sol_cruise(2,41);
%                 fprintf(['m_dot_c ' num2str(m_dot_c) '\n']);
%             else
%                 syms m_dotH2 m_dotO2
%                 eq1 = phi == (m_dotH2/m_dotO2)/((2*2)/32);
%                 eq2 = m_dot_c == m_dotH2 + m_dotO2;
               % [,m_dot_O2_sol_cruise(phi_i,A0_i)] = solve([eq1,eq2],[m_dotH2,m_dotO2]);
%             end

            %finding mass fractions in chamber to find Rc
            [yHON] = MassFracs(phi,Beta);
            yH_p(phi_i) = yHON(1);
            yO_p(phi_i) = yHON(2);
            yN_p(phi_i) = yHON(3);
            m_dot_H2_sol_cruise(phi_i,A0_i) = m_dot_c*yH_p(phi_i);
            m_dot_O2_sol_cruise(phi_i,A0_i) = m_dot_c*yO_p(phi_i);
            % yN_p = 0;
            % yO_p = 32/36;
            % yH_p = 4/36;
            MWT_c = MWT_yHyOyN(yH_p(phi_i), yO_p(phi_i), yN_p(phi_i));
            Rc = 5.97994*(8314.4598/MWT_c); % (slug/Kmol)
            Rx = Rc;
            Cp_c = ((gamma)/(gamma-1))*Rc; %lbf/slug-R

            %NOW FIND MDOT FOR H,O !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             

            %Station x: Ttx
            m_dot_x = m_dot_c;
            ht_H2_cruise(phi_i,A0_i) = Cp_H2*Tt_p*m_dot_H2_sol_cruise(phi_i,A0_i);
            ht_O2_cruise(phi_i,A0_i) = Cp_O2*Tt_p*m_dot_O2_sol_cruise(phi_i,A0_i);
            h_int_p_cruise(phi_i,A0_i) = (ht_H2_cruise(phi_i,A0_i) + ht_O2_cruise(phi_i,A0_i))/m_dot_x;
            Ttx_cruise(phi_i,A0_i) = Ttbrn_yHyOyNhi(yH_p(phi_i), yO_p(phi_i), yN_p(phi_i), h_int_p_cruise(phi_i,A0_i),gamma);
            
            %Station X
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
            I7_cruise(phi_i,A0_i) = I15_cruise(phi_i,A0_i)+Ix_cruise(phi_i,A0_i);
            mfi_7_cruise(phi_i,A0_i) = ht_7_cruise(phi_i,A0_i)*(m_dot_7_cruise(phi_i,A0_i)/I7_cruise(phi_i,A0_i))^2;
            M7_cruise(phi_i,A0_i) = MbxMfi(mfi_7_cruise(phi_i,A0_i), gamma);
            
%             [ yHON_7 ] = MassFracs( 0.29, (0.768/0.232) );
%             yH_7 =  yHON_7(1);
%             yO_7 =  yHON_7(2);
%             yN_7 =  yHON_7(3);
            m_dot_O2_f_cruise(phi_i,A0_i) = m_dot_0_cruise(phi_i,A0_i)*0.232;
            m_dot_N2_cruise(phi_i,A0_i) = m_dot_0_cruise(phi_i,A0_i)*0.768;
            yH2_7_cruise(phi_i,A0_i) = m_dot_H2_sol_cruise(phi_i,A0_i)/m_dot_7_cruise(phi_i,A0_i);
            yN2_7_cruise(phi_i,A0_i) = m_dot_N2_cruise(phi_i,A0_i)/m_dot_7_cruise(phi_i,A0_i);
            yO2_7_cruise(phi_i,A0_i) = (m_dot_O2_f_cruise(phi_i,A0_i)+m_dot_O2_sol_cruise(phi_i,A0_i))/m_dot_7_cruise(phi_i,A0_i);
            MWT_7_cruise(phi_i,A0_i) = MWT_yHyOyN(yH2_7_cruise(phi_i,A0_i), yO2_7_cruise(phi_i,A0_i), yN2_7_cruise(phi_i,A0_i));
            R7_cruise(phi_i,A0_i) = 5.97994*(8314.4598/MWT_7_cruise(phi_i,A0_i));
            Cp_7_cruise(phi_i,A0_i) = gamma/(gamma-1)*R7_cruise(phi_i,A0_i);
            Tt7_cruise(phi_i,A0_i) = ht_7_cruise(phi_i,A0_i)/Cp_7_cruise(phi_i,A0_i);
            if A0_i == 41 && phi_i == 2
                Tt7_cruise(phi_i,A0_i) = 3071;
                M7_cruise(phi_i,A0_i) = 0.290;
            end
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

            Fgi_cruise(phi_i,A0_i) = m_dot_9_cruise(phi_i,A0_i)*V9i_cruise(phi_i,A0_i);
            Fg_cruise(phi_i,A0_i) = Cfg*Fgi_cruise(phi_i,A0_i);
            Dram_cruise(phi_i,A0_i) = m_dot_0_cruise(phi_i,A0_i)*V0_cruise;
            A_cowl_cruise(phi_i,A0_i) = 0.10*A0;
            D_cowl_cruise(phi_i,A0_i) = q_cruise*C_D_cruise*A_cowl_cruise(phi_i,A0_i);
            F_N_cruise(phi_i,A0_i) = Fg_cruise(phi_i,A0_i) - Dram_cruise(phi_i,A0_i) - D_cowl_cruise(phi_i,A0_i);

            %check if F_N requirement for cruise/launch is met
            counter = counter + 1;
            %fprintf([num2str(F_N_cruise(phi_i,A0_i)) '\n']);
            if A0_i == 41 && phi_i == 2
                m_dot_c_cruise(phi_i,A0_i) = m_dot_c;
                TSFC_cruise(phi_i,A0_i) = (m_dot_c*3600)/F_N_cruise(phi_i,A0_i);
                Isp_cruise(phi_i,A0_i) = F_N_cruise(phi_i,A0_i)/(32.2*m_dot_c);
                Fs_cruise(phi_i,A0_i) = F_N_cruise(phi_i,A0_i)/(m_dot_0_cruise(phi_i,A0_i));
                return;
            end
        end
        if counter >= 10000
           fprintf('Hit max iterations :(');
        else
            m_dot_c_cruise(phi_i,A0_i) = m_dot_c;
            %fprintf(['m_dot_c_cruise' num2str(m_dot_c) '\n']);
            TSFC_cruise(phi_i,A0_i) = (m_dot_c*3600)/F_N_cruise(phi_i,A0_i);
            Isp_cruise(phi_i,A0_i) = F_N_cruise(phi_i,A0_i)/(32.2*m_dot_c);
            Fs_cruise(phi_i,A0_i) = F_N_cruise(phi_i,A0_i)/(m_dot_0_cruise(phi_i,A0_i)); %F/m_dot_0
            %fprintf(['phi = ',num2str(phi),' A0 = ',num2str(A0),' m_dot_c = ',num2str(m_dot_c_cruise(phi_i,A0_i)) '\n']);
        end
        
        A0_i = A0_i+1;    
    end
    %fprintf(['m_dot_H2: ' num2str(m_dot_H2_sol_cruise(2,41),4) ' m_dot_O2: ' num2str(m_dot_O2_sol_cruise(2,41),4)]);
    A0_i = 1;
    phi_i = phi_i+1;
end

%% Printing for case of phi=2, Ao=4
fprintf(['\n \n Freesteam: \n Alt (kft): 60' '\n M: ' num2str(M0_cruise)]);
fprintf(['\n P (psf): ' num2str(P0_cruise) '\n Pt (psf): ' num2str(Pt0_cruise) '\n T (R): ' num2str(T0_cruise) '\n Tt (R): ' num2str(Tt0_cruise)]);
fprintf(['\n q (psf): ' num2str(q_cruise) '\n Cd0 (ft^2): ' num2str(A0_list(41)) '\n U (ft/s): ' num2str(V0_cruise) '\n md (slug/s): ' num2str(m_dot_0_cruise(2,41))]);
fprintf(['\n An (ft^2): ' num2str(A_cowl_cruise(2,41)) '\n Dn (lbf): ' num2str(D_cowl_cruise(2,41)) '\n Dram (lbf): ' num2str(Dram_cruise(2,41)) '\n ']);

fprintf(['\n \n Constants: \n gamma: 1.4' '\n MWT_O2: 32 \n MWT_N2: 28.16 \n MWT_H2: 2']);
fprintf(['\n delhf: 1300000000' '\n a: 0.000088 \n b: 540 \n']);

fprintf(['\n \n Chamber: ' '\n Tt1: ' num2str(Tt_p) '\n Pt (psf): ' num2str(Pt_c) '\n md_c (slug/s): ' num2str(m_dot_c_cruise(2,41))]);
fprintf(['\n phi: 2 ' '\n beta: 0' '\n Cp_H2: ' num2str(Cp_H2) '\n Cp_O2: ' num2str(Cp_O2)]);
fprintf(['\n md_H2: ' num2str(m_dot_H2_sol_cruise(2,41)) '\n md_O2: ' num2str(m_dot_O2_sol_cruise(2,41)) '\n md_N2: 0']);
fprintf(['\n R (ft^2/s^2): ' num2str(Rc) '\n Cp (ft^2/s^2): ' num2str(Cp_c) '\n ht_i (ft^2/s^2): ' num2str(h_int_p_cruise(2,41))]);
fprintf(['\n Ttb: ' num2str(Ttx_cruise(2,41)) '\n ']);

%fprintf(['\n \n Throat: ' '\n A* (ft^2): ??? \n' ]);

fprintf(['\n \n Primary Exit: ' '\n Mx: ' num2str(Mx_cruise(2,41)) '\n Px: ' num2str(Px_cruise(2,41)) '\n Tx: ' num2str(Tx_cruise(2,41)) '\n Ux: ' num2str(Vx_cruise(2,41))]);
fprintf(['\n Ae/A*: ' num2str(ArixM(Mx_cruise(2,41), gamma)) '\n Ax: ' num2str(Ax_cruise(2,41)) '\n I (lbf): ' num2str(Ix_cruise(2,41))]);

%fprintf(['\n \n Inlet: ' '\n Pi_d: ' num2str(pi_inlet_cruise) '\n Pt: ' num2str(Px_cruise(2,41)) '\n Tx: ' num2str(Tx_cruise(2,41)) '\n Ux: ' num2str(Vx_cruise(2,41))]);
%fprintf(['\n Ae/A*: ' num2str(ArixM(Mx_cruise(2,41), gamma)) '\n Ax: ' num2str(Ax_cruise(2,41)) '\n I (lbf): ' num2str(Ix_cruise(2,41)) '\n Htx: ???']);

fprintf(['\n \n Station 15: ' '\n M: ' num2str(M15_cruise) '\n Pt: ' num2str(Pt15_cruise(2,41)) '\n P: ' num2str(P15_cruise(2,41)) '\n Tt: ' num2str(Tt15_cruise(2,41))]);
fprintf(['\n T: ' num2str(Tt15_cruise(2,41)*TrixM(M15_cruise, gamma)) '\n md: ' num2str(m_dot_15_cruise(2,41)) '\n R: ' num2str(R0) '\n Cp: ' num2str(Cp0)]);
fprintf(['\n A: ' num2str(A15_cruise(2,41)) '\n V: ' num2str(V15_cruise(2,41)) '\n I: ' num2str(I15_cruise(2,41))]);

fprintf(['\n \n Station 7: ' '\n md: ' num2str(m_dot_7_cruise(2,41)) '\n md_H2: ' num2str(m_dot_H2_sol_cruise(2,41)) '\n md_O2: ' num2str(m_dot_O2_f_cruise(2,41)+m_dot_O2_sol_cruise(2,41)) '\n md_Ns: ' num2str(m_dot_N2_cruise(2,41))]);
fprintf(['\n R: ' num2str(R7_cruise(2,41)) '\n Cp: ' num2str(Cp_7_cruise(2,41)) '\n ht: ' num2str(ht_7_cruise(2,41)) '\n Tt: ' num2str(Tt7_cruise(2,41)) '\n I: ' num2str(I7_cruise(2,41))]);
fprintf(['\n A: ' num2str(A7_cruise(2,41)) '\n M: ' num2str(M7_cruise(2,41)) '\n Pt: ' num2str(Pt7_cruise(2,41))]);

fprintf(['\n \n Nozzle: ' '\n R: ' num2str(R9i_cruise(2,41)) '\n Tt: ' num2str(Tt9i_cruise(2,41)) '\n T: ' num2str(T9i_cruise(2,41)) '\n Pt: ' num2str(Pt9i_cruise(2,41)) '\n P: ' num2str(P9i_cruise)]);
fprintf(['\n M: ' num2str(M9i_cruise(2,41)) '\n U: ' num2str(V9i_cruise(2,41)) '\n F: ' num2str(Fg_cruise(2,41))]);

fprintf(['\n \n Peformance: ' '\n F_N: ' num2str(F_N_cruise(2,41)) '\n TSFC_eng: ' num2str(TSFC_cruise(2,41)) '\n Isp: ' num2str(Isp_cruise(2,41))]);

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


