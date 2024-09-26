function u = unitNamer(u)
% A helper function to AnPSD_units and eyeAnalyses for constructing unit
% name string with trailing zeroes.

u = num2str(u);
if length(u) < 2
  u = ['000' u];
elseif length(u) < 3
  u = ['00' u];
elseif length(u) < 4
  u = ['0' u];
end