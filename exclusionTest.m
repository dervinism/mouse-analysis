function breakClause = exclusionTest(seriesName, exclusionsSeries)
% Check if a series should be excluded. A helper function to globalPCA,
% globalPCA2, globalPCA_noS1_noRSC, and globalUnitsFigs.

for iExclude = 1:numel(exclusionsSeries)
  if strcmpi(seriesName, exclusionsSeries{iExclude})
    breakClause = true;
    break
  else
    breakClause = false;
  end
end
end