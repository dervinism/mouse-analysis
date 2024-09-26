function [strOmnibus, strRao, strHR, strHartigan, strKDE] = summaryStatsPupilPhase(area, phaseData)


%% Fig01
[strOmnibus{1}, strRao{1}, strHR{1}, strHartigan{1}, strKDE{1}] = summaryStatsFOI(area, phaseData, 0.03);


%% Fig02
[strOmnibus{2}, strRao{2}, strHR{2}, strHartigan{2}, strKDE{2}] = summaryStatsFOI(area, phaseData, 0.3);



%% Local functions
function [strOmnibus, strRao, strHR, strHartigan, strKDE] = summaryStatsFOI(area, phaseData, FOI)
freq = phaseData.areaFreqFullInterpIndividual;
inds = freq == FOI;
mOmnibus = phaseData.distributionStats.mOmnibus{1}{determineArea(area)}{inds};
pOmnibus = phaseData.distributionStats.pOmnibus{1}{determineArea(area)}{inds};
U_Rao = phaseData.distributionStats.U_Rao{1}{determineArea(area)}{inds};
pRao = phaseData.distributionStats.pRao{1}{determineArea(area)}{inds};
T_HR = phaseData.distributionStats.T_HR{1}{determineArea(area)}{inds};
pHR = phaseData.distributionStats.pHR{1}{determineArea(area)}{inds};
dipHDT = phaseData.distributionStats.dipHDT{1}{determineArea(area)}{inds};
pHDT = phaseData.distributionStats.pHDT{1}{determineArea(area)}{inds};

for iFit = 1:numel(phaseData.distributionStats.nModes{1}{determineArea(area)}{inds})
  pKDE = phaseData.distributionStats.pKDE{1}{determineArea(area)}{inds}(iFit);
  if pKDE < 0.05
    break
  end
end
nModes = phaseData.distributionStats.nModes{1}{determineArea(area)}{inds}(iFit);
U2_KDE = phaseData.distributionStats.U2_KDE{1}{determineArea(area)}{inds}(iFit);

strOmnibus = sprintf('Omnibus M = %.3g p = %.3g \n', mOmnibus, pOmnibus);
strRao = sprintf('Rao U = %.3g p = %.3g \n', U_Rao, pRao);
strHR = sprintf('Hermans-Rasson T = %.3g p = %.3g \n', T_HR, pHR);
strHartigan = sprintf('Hartigan''s dip = %.3g p = %.3g \n', dipHDT, pHDT);
strKDE = sprintf('Kernel Density Excess number of modes = %.3g U2 = %.3g p = %.3g \n', nModes, U2_KDE, pKDE);