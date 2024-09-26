% Run this script to draw figures for phase and coherence comparisons
% between series, as well as to perform corresponding correlation analyses.


% LOAD PRE-PROCESSED DATA
if ~exist('dataStruct','var')
  load(dataFile);
end


% INITIALISE PARAMETERS
visibility = 'on';


% DRAW FIGURES
[figPhase, figCoherence, phaseVec, coherenceVec, cohConfVec] = halfCorrPlot(phaseCohSeries, {'awake', 'anaesthesia'}, visibility);


% CORRELATION ANALYSES
[phaseCohCorr.rPhase, phaseCohCorr.rhoPhase, phaseCohCorr.rhoCircPhase, phaseCohCorr.rCoherence,...
  phaseCohCorr.rhoCoherence] = corrStates(figPhase, figCoherence, phaseVec, coherenceVec, cohConfVec);


% SAVE FIGURES
initFigName = [dataFileSeries(1:end-4) '_ephys'];
saveFigsState(initFigName, figPhase, figCoherence, FOI);


% SAVE DATA
save(dataFileSeries,'-append','phaseCohCorr');

clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca dataStruct