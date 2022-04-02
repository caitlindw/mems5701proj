function [Mp] = MpxAri(ARin, g)
% Computes the supersonic Mach number for the given A/A*

  if nargin < 2, g = 1.4; end
  eps = 0.00000001;
  itmax = 1000;

  Ma = 1;
  Aa = ArixM(Ma, g);
  Mb = 10 ^ 5;
  Ab = ArixM(Mb, g);
  
  n = g / (g - 1);
    
  for iter = 1: itmax
    if mod(iter, 2) == 1
      Mg = 0.5 * (Ma + Mb);
    else
      Mg = ((ARin - Aa) * (Mb ^ n - Ma ^ n) / (Ab - Aa) + Ma ^ n) ^ (1 / n);
    end
    Ag = ArixM(Mg, g);
    
    if abs(Ag - ARin) <= eps
      break;
    elseif Ag > ARin
      Mb = Mg;
      Ab = Ag;
    else
      Ma = Mg;
      Aa = Ag;
    end
  end

  Mp = Mg;
end
