function condition = determineCondition(seriesName, awakeSeries, anaesthesiaSeries)
% Determine whether the series belongs to awake or anaesthesia conditions.
% A helper function to globalPCA, globalPCA2, globalPCA_noS1_noRSC, and
% globalUnitsFigs.

condition = 0;
for entry = 1:numel(awakeSeries)
  if (numel(seriesName) >= 14 && strcmpi(awakeSeries{entry}, seriesName(1:14))) ||...
      (numel(seriesName) < 14 && strcmpi(awakeSeries{entry}, seriesName))
    condition = 1;
  end
end
for entry = 1:numel(anaesthesiaSeries)
  if numel(seriesName) >= 14 && strcmpi(anaesthesiaSeries{entry}, seriesName(1:14)) ||...
      (numel(seriesName) < 14 && strcmpi(anaesthesiaSeries{entry}, seriesName))
    condition = 2;
  end
end
end