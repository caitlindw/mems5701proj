function [M] = MxTri(TR, g)
% Computes M for the given T/Tt.
  if nargin < 2, g = 1.4; end
  M = sqrt((2 / (g - 1)) * (1 / TR - 1));
end
