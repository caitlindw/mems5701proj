function [Mp] = MpxMfi(MFIin, g)
% Computes the supersonic Mach number for the given mass flow function, ht*(mdot/I)^2

  if nargin < 2, g = 1.4; end
  eps = 0.0000001;
  itmax = 1000;

  Ma = 1;
  Fa = MfixM(Ma, g);
  Mb = 10 ^ 5;
  Fb = MfixM(Mb, g);
    
  for iter = 1: itmax
    if mod(iter, 2) == 1
      Mg = 0.5 * (Ma + Mb);
    else
      Mg = 1 / sqrt(Ma ^ -2 + (MFIin - Fa) * (Mb ^ -2 - Ma ^ -2) / (Fb - Fa));
    end
    Fg = MfixM(Mg, g);
    
    if abs(Fg - MFIin) <= eps
      break;
    elseif Fg < MFIin
      Mb = Mg;
      Fb = Fg;
    else
      Ma = Mg;
      Fa = Fg;
    end
  end

  Mp = Mg;
end
