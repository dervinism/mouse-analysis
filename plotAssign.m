function [fH, areaVar, fLegend] = plotAssign(fH, dataVec, animalColour, fLegend, animal, areaVar)
% Storage variable assignment and plotting. A helper functon for globalPCA
% and globalPCA2.

figure(fH); hold on
p = plot(dataVec, 'color',animalColour); hold off
if fLegend.update
  fLegend.lines = [fLegend.lines p];
  if isempty(fLegend.text)
    fLegend.text = {animal};
  else
    fLegend.text{numel(fLegend.text)+1} = animal;
  end
  fLegend.update = false;
end
if isempty(areaVar)
  areaVar = dataVec;
else
  areaVarTemp = areaVar;
  areaVar = [areaVarTemp; dataVec];
end
end