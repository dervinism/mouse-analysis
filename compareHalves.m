% Run this script to perform phase and coherence comparisons between halves
% of the same recording. Not an actual analysis but rather extraction of
% relevant data that was already produced in AnPSD_units.


% LOAD PRE-PROCESSED DATA
if ~exist('dataStruct','var')
  load(dataFile);
end


%% INITIALISE PARAMETERS
params
lists
intermediateSaving = false; % Save everytime a db series is finished being analysed (bad idea since saving is a lengthy process)


%% INITIALISE STORAGE VARIABLES
phaseCohHalves.half1.phase = [];
phaseCohHalves.half1.coh = [];
phaseCohHalves.half1.coh_conf = [];
phaseCohHalves.half2.phase = [];
phaseCohHalves.half2.coh = [];
phaseCohHalves.half2.coh_conf = [];


%% PERFORM CORRELATIONS
fnsData = fieldnames(dataStruct.seriesData);
if ~isempty(dbEntries) && dbEntries(1) == inf % Here you can chose to execute only certain db entries
  dbEntriesLocal = 1:numel(fnsData);
else
  dbEntriesLocal = dbEntries;
end
for dbCount = dbEntriesLocal % Loop through db entries
  
  % Load the contents of dbStruct
  [dbStruct, ~, ~, entryName, ~, ~, shankIDs,...
    ~, ~, ~, ~, ~, ~, ~, FOI] = get_dbStruct(dataStruct, dbCount);
  
  phaseCohHalvesDB = {};
  for half = 1:2 % Loop through halves
    
    for sh = 1:numel(shankIDs) % Loop through shanks
      fprintf('%s shank %d -------------------\n', fnsData{dbCount}, sh);
      
      [~, ~, units, ~, ~, ~, ~,...
        ~, phaseCoh] = get_shankStruct(dbStruct, sh);
      if isempty(units)
        continue
      end
      if isempty(phaseCoh)
        continue
      end
      
      nUnits = numel(phaseCoh);
      for u = 1:nUnits % Loop through units
        fprintf('Processing %s shank %d unit %d\n', fnsData{dbCount}, sh, u);
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
        state.series = fnsData{dbCount};
        
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
  end
  
  
  % Save data
  dataString = ['dataStruct.seriesData.' entryName '.popData.phaseCohHalves = phaseCohHalvesDB;'];
  eval(dataString);
  dataStruct.popData.phaseCohHalves = phaseCohHalves;
  if intermediateSaving
    save(dataFile,'dataStruct','-v7.3'); %#ok<*UNRCH>
  end
end

%% SAVE DATA
if ~intermediateSaving
  save(dataFile,'dataStruct','-v7.3');
end
clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca dataStruct