function legendParams = updateLegendsSeries(conditions, legendParams)
% Update series figure legends. A helper function to globalPCA2 and globalFR.

for iCond = 1:numel(conditions)
  legendParams{iCond}.update = true;
end
end