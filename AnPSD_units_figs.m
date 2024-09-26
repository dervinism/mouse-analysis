% A part of the old AnPSD script adapted to display and save figures of
% coherence and phase analyses for unit and population spiking rate data.


%% LOAD PRE-PROCESSED DATA
if ~exist('dataStruct','var')
  load(dataFile);
end


%% INITIALISE PARAMETERS
params
lists
figsubdirname = '';
UPF = 5; % Units Per Figure
x_lim = [0.01 120];
visibility = 'on';
figType1 = 'ephysUnit';
figType2 = 'ephysGroups';
powerUnits = 'APs^2/Hz';


%% DRAW FIGURES
fnsData = fieldnames(dataStruct.seriesData);
if ~isempty(dbEntries) && dbEntries(1) == inf
  dbEntriesLocal = 1:numel(fnsData);
else
  dbEntriesLocal = dbEntries;
end
for dbCount = dbEntriesLocal % Loop through db entries
  % Load the contents of dbStruct
  [dbStruct, ~, ~, entryName, ~, ~, shankIDs,...
    ~, period, ~, ~, srData, ~, ~, ~,...
    MUAsAll] = get_dbStruct(dataStruct, dbCount);
  
  % Indices for truncating data
  [inds, MUAsAll] = determineInds(period, srData, MUAsAll);
  
  % Output folder
  strSep = strfind(entryName,'s');
  figsubdirname = entryName(strSep+1:end);
  if ~exist(figsubdirname,'dir')
    mkdir(figsubdirname)
  end

  for sh = 1:numel(shankIDs) % Loop through shanks
    fprintf('%s shank %d -------------------\n', fnsData{dbCount}, sh);
    
    % Load the contents of shankStruct
    [shankStruct, shank, units, ~, ~, ~, ~, MUAs, phaseCoh] = get_shankStruct(dbStruct, sh);
    if isempty(units)
      continue
    end
    if ~isfield(shankStruct,'phaseCoh')
      continue
    end
    if ~isempty(MUAs)
      MUAs = MUAs(:,inds);
    end
    
    nUnits = numel(units);
    for u = 1:nUnits % Loop through units
      fprintf('Processing %s shank %d unit %d\n', fnsData{dbCount}, sh, u);
      unitData = phaseCoh{u};
      if isfield(unitData,'coh')
          phaseCoherencePlots(sh, u, UPF, unitData, srData, units(u), nUnits, powerUnits, x_lim, figsubdirname, entryName, figType1,...
            visibility, MUAsAll, shank, MUAs);
      end
    end % loop over units
  end % loop over shanks
end % loop over db entries

clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca dataStruct