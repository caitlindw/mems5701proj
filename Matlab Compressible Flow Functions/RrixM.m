function [Rri] = RrixM(M, g)
% Computes Rho/Rhot for the given M.
  if nargin < 2, g = 1.4; end
  Rri = (1 + (g - 1) / 2 * M ^ 2) ^ (1 / (1 - g));
end
