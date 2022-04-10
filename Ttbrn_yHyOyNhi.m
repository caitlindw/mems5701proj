function [ Ttbrn ] = Ttbrn_yHyOyNhi(yH, yO, yN, hi, g)
% Approximates the gas total temperature after combution of hydrogen
% including a crude estimate to account for dissociation.
% Inputs: yH  = mass fraction of hydrogen atoms
%         yO  = mass fraction of oxygen atoms
%         yN  = mass fraction of inert atoms
%         hi  = specific enthalpy of reactants, not including heats of formation [ft^2/s^2]
% 
% 'Combustion of hydrogen approximated by:
%         hpr = heat of combustion = 1.3e9 [ft^2/s^2]
%         a   = combustion efficiency/dissociation sensitivity = 8.8e-5
%         b   = combustion efficiency/dissociation constant = 540

% Set constants
  hpr = 1300000000;
  a = 0.000088;
  b = 540;
 
% Ensure mass fractions add to one and get the combusted gas properties.
  Ytot = yH + yO + yN;
  yH = yH / Ytot;
  yO = yO / Ytot;
  yN = yN / Ytot;

  MWT = MWT_yHyOyN(yH, yO, yN);
  R = 49710. / MWT;
  Cp = R * g / (g - 1);

% Compute the mass fraction of burned hydrogen and the mixture temperature.
  Yh2b = min(yH, yO / 8);
  
  Ttbrn = (hi + Yh2b * hpr * (1 + a * b)) / (Cp + a * Yh2b * hpr);

end

