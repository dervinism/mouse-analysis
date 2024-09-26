function str = summaryStatsPupilHalfPhase(area, phaseData)


%% Fig01
[str{1}] = summaryStatsFOI(area, phaseData, 0.03);


%% Fig02
[str{2}] = summaryStatsFOI(area, phaseData, 0.3);



%% Local functions
function str = summaryStatsFOI(area, phaseData, FOI)
freq = phaseData.areaFreqFullInterpIndividual;
inds = freq == FOI;
rho = phaseData.rPhaseHalves{1}{determineArea(area)}(inds);
pval = phaseData.pvalPhaseHalves{1}{determineArea(area)}(inds);
str = sprintf('Spearman rho = %.3g p = %.3g \n', rho, pval);