function [Mb] = MbxAri(ARin, g)
% Computes the subsonic Mach number for the given A/A*

  if nargin < 2, g = 1.4; end
  eps = 0.0000001;
  itmax = 1000;

  Ma = eps;
  Aa = ArixM(Ma, g);
  Mb = 1;
  Ab = ArixM(Mb, g);
    
  for iter = 1: itmax;
    if mod(iter, 2) == 1
      Mg = 0.5 * (Ma + Mb);
    else
      Mg = 1 / ((ARin - Aa) * (1 / Mb - 1 / Ma) / (Ab - Aa) + 1 / Ma);
    end
    Ag = ArixM(Mg, g);
    
    if abs(Ag - ARin) <= eps
      break;
    elseif Ag < ARin
      Mb = Mg;
      Ab = Ag;
    else
      Ma = Mg;
      Aa = Ag;
    end
  end

  Mb = Mg;
end
