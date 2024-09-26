function legendParams = initLegendsSeries(conditions)
% Initialise series figure legends. A helper function to globalPCA2 and globalFR.

for iCond = 1:numel(conditions)
  legendParams{iCond}.lines = []; %#ok<*AGROW>
  legendParams{iCond}.text = {};
  legendParams{iCond}.update = true;
end
end