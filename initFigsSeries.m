function seriesFigs = initFigsSeries(conditions)
% Initialise series figures. A helper function to globalPCA2 and globalFR.

for iCond = 1:numel(conditions)
  seriesFigs{iCond} = figure; %#ok<*AGROW>
end
end