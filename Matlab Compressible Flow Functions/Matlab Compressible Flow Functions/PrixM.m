function [Pri] = PrixM(M, g)
% Computes P/Pt for the given M.
  if nargin < 2, g = 1.4; end
  Pri = (1 + (g - 1) / 2 * M ^ 2) ^ (g / (1 - g));
end
