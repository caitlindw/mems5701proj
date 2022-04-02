function [M] = MxPri(PR, g)
% Computes M for the given P/Pt.
  if nargin < 2, g = 1.4; end
  M = (2 / (g - 1) * (PR ^ ((1 - g) / g) - 1)) ^ 0.5;
end
