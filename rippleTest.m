function passClause = rippleTest(series, fullSeries)
% Test if series contains ripples

passClause = false;
for s = 1:numel(fullSeries)
  if strcmpi(series, fullSeries{s})
    passClause = true;
    break
  end
end