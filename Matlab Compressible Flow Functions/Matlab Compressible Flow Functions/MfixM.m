function [Mfi] = MfixM(M, g)
% Computes ht*(mdot/I)^2 for the given M.
  if nargin < 2, g = 1.4; end
  Mfi = g ^ 2 / (g - 1) * M ^ 2 * (1 + (g - 1) / 2 * M ^ 2) / (1 + g * M ^ 2) ^ 2;
end
