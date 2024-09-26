for iCond = 1:numel(areaBetaWindowsIndividual)
  for iArea = 1:numel(areaBetaWindowsIndividual{iCond})
    betas = [];
    for iBeta = 1:numel(areaBetaWindowsIndividual{iCond}{iArea})
      betas = [betas; areaBetaWindowsIndividual{iCond}{iArea}{iBeta}]; %#ok<*AGROW>
    end
    areaBetaWindowsIndividual{iCond}{iArea} = betas; %#ok<*SAGROW>
  end
end