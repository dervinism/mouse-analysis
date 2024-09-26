% Run this script to perform phase and coherence comparisons between halves
% of the same recording for area comparisons.


%% LOAD PRE-PROCESSED DATA
if ~exist('dataStruct','var')
  load(dataFile);
end


%% INITIALISE PARAMETERS
params
lists
intermediateSaving = false;


%% INITIALISE STORAGE VARIABLES
phaseCohHalves.half1.phase = [];
phaseCohHalves.half1.coh = [];
phaseCohHalves.half1.coh_conf = [];
phaseCohHalves.half2.phase = [];
phaseCohHalves.half2.coh = [];
phaseCohHalves.half2.coh_conf = [];


%% PERFORM CORRELATIONS
fnsData = fieldnames(dataStruct.seriesData);
fnsData_ca = fieldnames(dataStruct.seriesData_ca_muas);
if ~isempty(dbEntries_ca) && dbEntries_ca(1) == inf
  dbEntries_caLocal = 1:numel(fnsData_ca);
else
  dbEntries_caLocal = dbEntries_ca;
end
for dbCount = dbEntries_caLocal % Loop through db entries
  
  % Load the contents of dbStruct
  dbStruct = dataStruct.seriesData_ca_muas.(fnsData_ca{dbCount});
  if isempty(dbStruct)
    continue
  end
  seriesName1 = seriesNames(fnsData_ca{dbCount});
  [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, FOI] = get_dbStruct(dataStruct, find(endsWith(fnsData,seriesName1)));
  
  phaseCohHalvesDB = {};
  for half = 1:2 % Loop through halves
    
    units = dbStruct.units;
    phaseCoh = dbStruct.phaseCoh;
    if isempty(units)
      continue
    end
    if isempty(phaseCoh)
      continue
    end
    nUnits = numel(phaseCoh);
    for u = 1:nUnits % Loop through units
      fprintf('Processing %s unit %d\n', fnsData_ca{dbCount}, u);
      unitData = phaseCoh{u};
      
      if isfield(unitData,'coh_halves') && ~isempty(unitData.coh_halves)
        [state.phase, state.coh, state.coh_conf, state.FOI, state.fInds] = phaseCohFOI(FOI, unitData.coh_halves_freq,...
          unitData.phase_halves(half,:), unitData.coh_halves(half,:), unitData.coh_conf_halves(half,:),...
          unitData.rateadjust_kappa_halves(half,:));
      else
        state.phase = nan(size(FOI));
        state.coh = nan(size(FOI));
        state.coh_conf = nan(size(FOI));
        state.FOI = FOI;
        state.fInds = nan(size(FOI));
      end
      state.unit = units(u);
      state.series = fnsData_ca{dbCount};
      
      if half == 1
        phaseCohHalvesDB{u}.half1 = state; %#ok<*SAGROW>
        phaseCohHalves.half1.phase = [phaseCohHalves.half1.phase; state.phase];
        phaseCohHalves.half1.coh = [phaseCohHalves.half1.coh; state.coh];
        phaseCohHalves.half1.coh_conf = [phaseCohHalves.half1.coh_conf; state.coh_conf];
      elseif half == 2
        phaseCohHalvesDB{u}.half2 = state;
        phaseCohHalves.half2.phase = [phaseCohHalves.half2.phase; state.phase];
        phaseCohHalves.half2.coh = [phaseCohHalves.half2.coh; state.coh];
        phaseCohHalves.half2.coh_conf = [phaseCohHalves.half2.coh_conf; state.coh_conf];
      end
    end
  end
  
  % Save data
  dataString = ['dataStruct.seriesData_ca_muas.' fnsData_ca{dbCount} '.popData.phaseCohHalves = phaseCohHalvesDB;'];
  eval(dataString);
  dataStruct.popData_ca_muas.phaseCohHalves = phaseCohHalves;
  if intermediateSaving
    save(dataFile,'dataStruct','-v7.3'); %#ok<*UNRCH>
  end
end

if ~intermediateSaving
  save(dataFile,'dataStruct','-v7.3');
end
clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca dataStruct