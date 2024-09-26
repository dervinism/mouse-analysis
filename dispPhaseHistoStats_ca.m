function [fPEst, fWTest, pEst, U2, pObs, U2Obs] = dispPhaseHistoStats_ca(phaseData, cond, comp, frequency)

% Find the comparison index
vInd = strfind(comp, 'v');
comp1 = comp(1:vInd-1);
comp2 = comp(vInd+1:end);
for iComp1 = 1:numel(phaseData.comparisons.comparison1{cond})
  if strcmpi(phaseData.comparisons.comparison1{cond}{iComp1}, comp1)
    for iComp2 = 1:numel(phaseData.comparisons.comparison2{cond}{iComp1})
      if strcmpi(phaseData.comparisons.comparison2{cond}{iComp1}{iComp2}, comp2)
        break
      end
    end
    break
  end
end

% Extract the stats
fPEst = NaN(size(frequency));
fWTest = NaN(size(frequency));
pEst = NaN(size(frequency));
U2 = NaN(size(frequency));
pObs = NaN(size(frequency));
U2Obs = NaN(size(frequency));
for iF = 1:numel(frequency)
  if ~isempty(phaseData.fPEst_ca{cond}{iComp1}{iComp2}(iF))
    fPEst(iF) = phaseData.fPEst_ca{cond}{iComp1}{iComp2}(iF);
  end
  if ~isempty(phaseData.fWTest_ca{cond}{iComp1}{iComp2}(iF))
    fWTest(iF) = phaseData.fWTest_ca{cond}{iComp1}{iComp2}(iF);
  end
  if ~isempty(phaseData.pEst_ca{cond}{iComp1}{iComp2}(iF))
    pEst(iF) = phaseData.pEst_ca{cond}{iComp1}{iComp2}(iF);
  end
  if ~isempty(phaseData.U2_ca{cond}{iComp1}{iComp2}(iF))
    U2(iF) = phaseData.U2_ca{cond}{iComp1}{iComp2}(iF);
  end
  if ~isempty(phaseData.pObs_ca{cond}{iComp1}{iComp2}(iF))
    pObs(iF) = phaseData.pObs_ca{cond}{iComp1}{iComp2}(iF);
  end
  if ~isempty(phaseData.U2Obs_ca{cond}{iComp1}{iComp2}(iF))
    U2Obs(iF) = phaseData.U2Obs_ca{cond}{iComp1}{iComp2}(iF);
  end
end