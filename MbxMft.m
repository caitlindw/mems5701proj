function [Mb] = MbxMft(MFtin, g)
% Computes the subsonic Mach number for the given mass flow function, mdot(RTt)^.5/(PtA)

  if nargin < 2, g = 1.4; end
  eps = 0.0000001;
  itmax = 1000;

  Ma = 0;
  Fa = MftxM(Ma, g);
  Mb = 1;
  Fb = MftxM(Mb, g);
    
  for iter = 1: itmax
    if mod(iter, 2) == 1
      Mg = 0.5 * (Ma + Mb);
    else
      Mg = Ma + (MFtin - Fa) * (Mb - Ma) / (Fb - Fa);
    end
    Fg = MftxM(Mg, g);
    
    if abs(Fg - MFtin) <= eps
      break;
    elseif Fg > MFtin
      Mb = Mg;
      Fb = Fg;
    else
      Ma = Mg;
      Fa = Fg;
    end
  end

  Mb = Mg;
end
