% Run this script to draw figures after total motion analyses were performed.


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
powerUnitsMotion = '(area units)^2/Hz';
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
    motionDB = dbStruct.popData.motion;
  catch
    disp(['No pupil data for ' fnsData{dbCount} '. Skipping to the next db entry...']);
    continue
  end
  
  % Trim data
  [seriesName, animal] = seriesFromEntry(entryName);
  entryNameMotion = [animal '_s' seriesName(1:14)];
  commonPeriod = combinePeriods(period, dataStruct.motionData.(entryNameMotion).period, srData);
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
  phaseCoherencePlots(0, 1, 1, motionDB.popData, srData, 'mua', 1, powerUnitsMotion, x_lim, figsubdirname, entryName,...
    figType1, visibility, MUAsAll, numel(shankIDs), MUAsAll, true);

% LOOP THROUGH SHANKS
  for sh = 1:numel(shankIDs)
    shankStruct = eval(['dbStruct.shankData.' shankIDs{sh}]);
    fprintf('%s shank %d -------------------\n', fnsData{dbCount}, sh);
    
% LOAD THE CONTENTS OF THE SHANK STRUCTURE VARIABLE
    shank = shankStruct.shankID;
    MUAs = shankStruct.MUAs;
    if isfield(shankStruct,'motion')
      motion = shankStruct.motion;
    else
      continue
    end
    
% LOOP THROUGH UNITS
    nUnits = numel(motion.unitData);
    for u = 1:nUnits
      fprintf('Processing %s shank %d unit %d\n', fnsData{dbCount}, sh, u);
      unitData = motion.unitData{u};
      
% PLOT UNIT DATA
      phaseCoherencePlots(sh, u, UPF, unitData, srData, unitData.unit, nUnits, powerUnitsEphys, x_lim, figsubdirname,...
        entryName, figType2, visibility, MUAsAll, shank, MUAs);
      
    end % loop over units
  end % loop over shanks
  
% PLOT CORRELATIONS
  if isfield(dbStruct.popData, 'phaseCoh')
    figPhase = phaseVphasePlot(phaseCoh.phaseFOI, motionDB.unitData.phaseFOI, FOI, 'Spiking', 'Motion', visibility);
    [rPhaseMotVsPop, pvalPhaseMotVsPop] = corrMulti(motionDB.unitData.phaseFOI, dbStruct.popData.phaseCoh.phaseFOI, 'circular');
    correctedCohMot = motionDB.unitData.cohFOI;
    correctedCohMot(correctedCohMot - motionDB.unitData.coh_confFOI <= 0) = NaN;
    correctedCohPR = phaseCoh.cohFOI;
    correctedCohPR(correctedCohPR - phaseCoh.coh_confFOI <= 0) = NaN;
    figCoherence = cohVcohPlot(correctedCohPR, correctedCohMot, FOI, 'Spiking', 'Motion', visibility);
    figName = [figsubdirname filesep entryName '_' figType2 'Full'];
    saveFigsState(figName, figPhase, figCoherence, FOI);
    %close(figPhase); close(figCoherence);
  end
  
end % loop over db entries

clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca dataStruct
