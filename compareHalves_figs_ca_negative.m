% Run this script to draw figures for phase and coherence comparisons
% between half files, as well as to perform corresponding correlation
% analyses for area comparisons (within negative subpopulation only).


%% LOAD PRE-PROCESSED DATA
if ~exist('dataStruct','var')
  load(dataFile);
end


%% INITIALISE PARAMETERS
params
lists
visibility = 'on';


%% DRAW FIGURES
fnsData = fieldnames(dataStruct.seriesData_negative);
fnsData_ca = fieldnames(dataStruct.seriesData_ca_negative);
if (~isempty(dbEntries_ca) && dbEntries_ca(1) == inf) ||...
    (~isempty(dbEntries_ca) && dbEntries_ca(end) > numel(fnsData_ca))
  dbEntries_ca_negative = 1:numel(fnsData_ca);
else
  dbEntries_ca_negative = dbEntries_ca;
end
for dbCount = dbEntries_ca_negative % Loop through db entries
  dbStruct = dataStruct.seriesData_ca_negative.(fnsData_ca{dbCount});
  if ~isfield(dbStruct.popData, 'phaseCohHalves')
    continue
  end
  phaseCohHalves = dbStruct.popData.phaseCohHalves;
  entryName = fnsData_ca{dbCount};
  FOI = dataStruct.seriesData_negative.(fnsData{1}).conf.FOI;
  if isempty(phaseCohHalves)
    continue
  end
  
  % Output folder
  strSep = strfind(entryName,'__');
  entryName1 = entryName(1:strSep-1);
  strSep1 = strfind(entryName1,'s');
  entryName1 = entryName1(strSep1+1:end);
  entryName2 = entryName(strSep+2:end);
  strSep2 = strfind(entryName2,'s');
  entryName2 = entryName2(strSep2+1:end);
  figsubdirname = ['negative' filesep entryName1 'v' entryName2];
  if ~exist(figsubdirname,'dir')
    mkdir(figsubdirname)
  end
  
  % Figures
  [figPhase, figCoherence, phaseVec, coherenceVec, cohConfVec] = halfCorrPlot(phaseCohHalves, {'half1', 'half2'}, visibility);
  
  
  % Correlation analyses for db
  [phaseCohHalvesCorr.rPhase, phaseCohHalvesCorr.rhoPhase, phaseCohHalvesCorr.rhoCircPhase, phaseCohHalvesCorr.rCoh,...
    phaseCohHalvesCorr.rhoCoh] = corrStates(figPhase, figCoherence, phaseVec, coherenceVec, cohConfVec);
  
  
  % Save figures
  initFigName = [figsubdirname filesep entryName '_ephys_correlation_halves'];
  saveFigsState(initFigName, figPhase, figCoherence, FOI);
  %close(figPhase); close(figCoherence);
  
  
  % Save data
  dataString = ['dataStruct.seriesData_ca_negative.' entryName '.popData.phaseCohHalves = phaseCohHalves;'];
  eval(dataString);
  dataString = ['dataStruct.seriesData_ca_negative.' entryName '.popData.phaseCohHalvesCorr = phaseCohHalvesCorr;'];
  eval(dataString);
end


%% CORRELATION ANALYSES FOR ANIMAL
phaseCohHalves = dataStruct.popData_ca_negative.phaseCohHalves;
if ~isempty(phaseCohHalves.half1.phase)
  figPhase = phaseVphasePlot(phaseCohHalves.half1.phase, phaseCohHalves.half2.phase, FOI, 'Half1', 'Half2', visibility);
  correctedCoh1 = phaseCohHalves.half1.coh;
  correctedCoh1(correctedCoh1 - phaseCohHalves.half1.coh_conf <= 0) = NaN;
  correctedCoh2 = phaseCohHalves.half2.coh;
  correctedCoh2(correctedCoh2 - phaseCohHalves.half2.coh_conf <= 0) = NaN;
  figCoherence = cohVcohPlot(correctedCoh1,correctedCoh2, FOI, 'Half1', 'Half2', visibility);
  for iF = 1:numel(FOI)
    phaseVec{iF} = [phaseCohHalves.half1.phase; phaseCohHalves.half2.phase];
    coherenceVec{iF} = [phaseCohHalves.half1.coh; phaseCohHalves.half2.coh];
    cohConfVec{iF} = [phaseCohHalves.half1.coh_conf; phaseCohHalves.half2.coh_conf];
  end
  [phaseCohHalves.rPhase, phaseCohHalves.rhoPhase, phaseCohHalves.rhoCircPhase, phaseCohHalves.rCoh,...
    phaseCohHalves.rhoCoh] = corrStates(figPhase, figCoherence, phaseVec, coherenceVec, cohConfVec);
  close(figPhase); close(figCoherence);
  dataStruct.popData_ca_negative.phaseCohHalves = phaseCohHalves;
  save(dataFile,'dataStruct','-v7.3');
end

clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca dataStruct
