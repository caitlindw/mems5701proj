function [Ari] = ArixM(M, g)
% Computes A/A* for the given M.
  if nargin < 2, g = 1.4; end
  Ari = (2 * (1 + (g - 1) / 2 * M ^ 2) / (g + 1)) ^ ((g + 1) / (2 * (g - 1))) / M;
end
