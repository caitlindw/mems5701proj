function [Mft] = MftxM(M, g)
% Computes mdot(RTt)^.5/(PtA) for the given M.
  if nargin < 2, g = 1.4; end
  Mft = g ^ 0.5 * M * (1 + (g - 1) / 2 * M ^ 2) ^ (-(g + 1) / (2 * (g - 1)));
end
