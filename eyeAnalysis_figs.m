% Run this script to draw figures after eye analyses were performed.


%% LOAD PRE-PROCESSED DATA
if ~exist('dataStruct','var')
  load(dataFile);
end


%% INITIALISE PARAMETERS
params
lists
UPF = 5; % Units Per Figure
x_lim = [0.01 2];
visibility = 'on';
figType1 = 'pupilPop';
figType2 = 'pupilUnit';
powerUnitsPupil = '(area units)^2/Hz';
powerUnitsEphys = 'APs^2/Hz';


%% DRAW THE FIGURES
fnsData = fieldnames(dataStruct.seriesData);
if ~isempty(dbEntries) && dbEntries(1) == inf
  dbEntriesLocal = 1:numel(fnsData);
else
  dbEntriesLocal = dbEntries;
end
for dbCount = dbEntriesLocal % Loop through dbEntries
  
  % Load the contents of dbStruct
  [dbStruct, ~, ~, entryName, ~, ~, shankIDs,...
    ~, period, ~, ~, srData, ~, ~, FOI,...
    MUAsAll, ~, ~, phaseCoh] = get_dbStruct(dataStruct, dbCount);
  if isempty(MUAsAll) || ~sum(sum(MUAsAll))
    disp(['No spiking data for ' fnsData{dbCount} '. Skippig to the next db entry...']);
    continue
  end
  if isempty(phaseCoh)
    disp(['No coherence data for ' fnsData{dbCount} '. Skippig to the next db entry...']);
    continue
  end
  try
    pupilDB = dbStruct.popData.pupil;
  catch
    disp(['No pupil data for ' fnsData{dbCount} '. Skipping to the next db entry...']);
    continue
  end
  
  % Trim data
  [seriesName, animal] = seriesFromEntry(entryName);
  entryNameEye = [animal '_s' seriesName(1:14)];
  commonPeriod = combinePeriods(period, dataStruct.eyeData.(entryNameEye).period, srData);
  if isempty(commonPeriod)
    continue
  end
  [inds, MUAsAll] = determineInds(commonPeriod, srData, MUAsAll);
  
  % Output folder
  strSep = strfind(entryName,'s');
  figsubdirname = entryName(strSep+1:end);
  if ~exist(figsubdirname,'dir')
    mkdir(figsubdirname)
  end
  
% PLOT POPULATION DATA
  phaseCoherencePlots(0, 1, 1, pupilDB.popData, srData, 'mua', 1, powerUnitsPupil, x_lim, figsubdirname, entryName,...
    figType1, visibility, MUAsAll, numel(shankIDs), MUAsAll, true);

% LOOP THROUGH SHANKS
  for sh = 1:numel(shankIDs)
    shankStruct = eval(['dbStruct.shankData.' shankIDs{sh}]);
    fprintf('%s shank %d -------------------\n', fnsData{dbCount}, sh);
    
% LOAD THE CONTENTS OF THE SHANK STRUCTURE VARIABLE
    shank = shankStruct.shankID;
    MUAs = shankStruct.MUAs;
    if isfield(shankStruct,'pupil')
      pupil = shankStruct.pupil;
    else
      continue
    end
    
% LOOP THROUGH UNITS
    nUnits = numel(pupil.unitData);
    for u = 1:nUnits
      fprintf('Processing %s shank %d unit %d\n', fnsData{dbCount}, sh, u);
      unitData = pupil.unitData{u};
      
% PLOT UNIT DATA
      phaseCoherencePlots(sh, u, UPF, unitData, srData, unitData.unit, nUnits, powerUnitsEphys, x_lim, figsubdirname,...
        entryName, figType2, visibility, MUAsAll, shank, MUAs);
      
    end % loop over units
  end % loop over shanks
  
% PLOT CORRELATIONS
  if isfield(dbStruct.popData, 'phaseCoh')
    figPhase = phaseVphasePlot(phaseCoh.phaseFOI, pupilDB.unitData.phaseFOI, FOI, 'Spiking', 'Pupil', visibility);
    [rPhasePupVsPop, pvalPhasePupVsPop] = corrMulti(pupilDB.unitData.phaseFOI, dbStruct.popData.phaseCoh.phaseFOI, 'circular');
    correctedCohPup = pupilDB.unitData.cohFOI;
    correctedCohPup(correctedCohPup - pupilDB.unitData.coh_confFOI <= 0) = NaN;
    correctedCohPR = phaseCoh.cohFOI;
    correctedCohPR(correctedCohPR - phaseCoh.coh_confFOI <= 0) = NaN;
    figCoherence = cohVcohPlot(correctedCohPR, correctedCohPup, FOI, 'Spiking', 'Pupil', visibility);
    figName = [figsubdirname filesep entryName '_' figType2 'Full'];
    saveFigsState(figName, figPhase, figCoherence, FOI);
    %close(figPhase); close(figCoherence);
  end
  
end % loop over db entries

clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca dataStruct