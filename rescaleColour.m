function rescaleColour(scaling_factor)

if nargin < 1
  scaling_factor = 1.15;
end

c = colormap(parula);
rescaled_c = c;

for i = 1:size(c, 1)
  rescaled_c(i,:) = c(min(round(i^scaling_factor), size(c, 1)), :);
end

colormap(rescaled_c);