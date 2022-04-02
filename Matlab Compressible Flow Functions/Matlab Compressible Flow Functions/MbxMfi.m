function [Mb] = MbxMfi(MFIin, g)
% Computes the subsonic Mach number for the given mass flow function, ht*(mdot/I)^2

  if nargin < 2, g = 1.4; end
  eps = 0.0000001;
  itmax = 1000;

  Ma = 0;
  Fa = MfixM(Ma, g);
  Mb = 1;
  Fb = MfixM(Mb, g);
    
  for iter = 1: itmax
    if mod(iter, 2) == 1
      Mg = 0.5 * (Ma + Mb);
    else
      Mg = Ma + (sqrt(MFIin) - sqrt(Fa)) * (Mb - Ma) / (sqrt(Fb) - sqrt(Fa));
    end
    Fg = MfixM(Mg, g);
    
    if abs(Fg - MFIin) <= eps
      break;
    elseif Fg > MFIin
      Mb = Mg;
      Fb = Fg;
    else
      Ma = Mg;
      Fa = Fg;
    end
  end

  Mb = Mg;
end
