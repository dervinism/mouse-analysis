% Run this script to draw figures for phase and coherence comparisons
% between half files, as well as to perform corresponding correlation
% analyses.


%% LOAD PRE-PROCESSED DATA
if ~exist('dataStruct','var')
  load(dataFile);
end


%% INITIALISE PARAMETERS
params
lists
visibility = 'on';


%% DRAW FIGURES
fnsData = fieldnames(dataStruct.seriesData);
if ~isempty(dbEntries) && dbEntries(1) == inf
  dbEntriesLocal = 1:numel(fnsData);
else
  dbEntriesLocal = dbEntries;
end
for dbCount = dbEntriesLocal % Loop through db entries
  
  % Load the contents of dbStruct
  [~, ~, ~, entryName, ~, ~, ~,...
    ~, ~, ~, ~, ~, ~, ~, FOI,...
    ~, ~, ~, ~, phaseCohHalves] = get_dbStruct(dataStruct, dbCount);
  if isempty(phaseCohHalves)
    continue
  end
  
  % Output folder
  strSep = strfind(entryName,'s');
  figsubdirname = entryName(strSep+1:end);
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
  dataString = ['dataStruct.seriesData.' entryName '.popData.phaseCohHalves = phaseCohHalves;'];
  eval(dataString);
  dataString = ['dataStruct.seriesData.' entryName '.popData.phaseCohHalvesCorr = phaseCohHalvesCorr;'];
  eval(dataString);
end

%% Correlation analyses for animal
phaseCohHalves = dataStruct.popData.phaseCohHalves;
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
  dataStruct.popData.phaseCohHalves = phaseCohHalves;
  save(dataFile,'dataStruct','-v7.3');
end

clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca dataStruct