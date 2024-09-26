% A part of the old AnPSD script adapted to perform coherence and phase
% analyses for unit and population spiking rate data comparing cortical
% layers.


%% LOAD PRE-PROCESSED DATA
if ~exist('dataStruct','var')
  load(dataFile);
end


%% INITIALISE PARAMETERS
params
lists
intermediateSaving = false; % Save everytime a db series is finished being analysed (bad idea since saving is a lengthy process)


%% PERFORM COHERENCE ANALSYSES
fnsData = fieldnames(dataStruct.seriesData);
if ~isempty(dbEntries) && dbEntries(1) == inf % Here you can chose to execute only certain db entries
  dbEntriesLocal = 1:numel(fnsData);
else
  dbEntriesLocal = dbEntries;
end
for dbCount = dbEntriesLocal % Loop through db entries
  
  % Load the contents of dbStruct
  [dbStruct, repository, ~, entryName, ~, ~, shankIDs,...
    ~, period, ~, ~, srData, optCoh, exclRad, FOI,...
    MUAsAll, spkDB, spkDB_units] = get_dbStruct(dataStruct, dbCount);
  optCoh.maxFreq = maxFreq;
  optCoh.winfactor = winfactor;
  optCoh.freqfactor = freqfactor;
  optCoh.monotoneFreq = true;
  if numel(shankIDs) > 1
    continue
  end
  [~, ~, area] = determineArea(seriesFromEntry(entryName));
  if ~contains(area,'S1') && ~contains(area,'V1') && ~contains(area,'RSC')
    continue
  end
  if size(spkDB_units,2) > 1
    spkDB_units = spkDB_units(:,2);
  end
  
  % Indices for truncating data
  [inds, spkDB] = determineInds(period, srData, spkDB);
  
  % Load the contents of shankStruct
  [~, shank, units, unitMetadata, xcoords, ycoords, spk] = get_shankStruct(dbStruct, 1);
  if isempty(units)
    continue
  end
  
  % Truncate shank data
  if ~isempty(spk)
    if size(spk,2) < size(MUAsAll,2)
      spk = [full(spk) zeros(size(spk,1), size(MUAsAll,2)-size(spk,2))];
    end
    spk = spk(:,inds);
  end
  if strcmpi(repository, 'allensdk')
    divCh = mean([min(dbStruct.popData.muaMetadata(:,9)) max(dbStruct.popData.muaMetadata(:,9))]);
    MUAsAll1 = sum(spkDB(dbStruct.popData.muaMetadata(:,9) <= divCh,:), 1);
    MUAsAll2 = sum(spkDB(dbStruct.popData.muaMetadata(:,9) > divCh,:), 1);
  else
    divCh = min(dbStruct.dbSeries.chOI) + round(numel(dbStruct.dbSeries.chOI)/2) - 1;
    MUAsAll1 = sum(spkDB(dbStruct.popData.muaMetadata(:,3) <= divCh,:), 1);
    MUAsAll2 = sum(spkDB(dbStruct.popData.muaMetadata(:,3) > divCh,:), 1);
  end
  
  parfor u = 1:numel(units) % Loop through units
  %for u = 1:numel(units)
    fprintf('Started processing unit %i\n',units(u));
    
    unitData = struct();
    unitData.unit = units(u);
    if strcmpi(repository, 'allensdk')
      unitCh = unitMetadata(u,9);
    else
      unitCh = unitMetadata(u,3); %#ok<*PFBNS>
    end
    spkOI = full(spk(u,:));
    
    if unitCh <= divCh
      PR = full(MUAsAll2);
    else
      PR = full(MUAsAll1);
    end
    
    % Summary data
    [unitData.mfr, unitData.mfr_1sthalf, unitData.mfr_2ndhalf, unitData.lfr1, unitData.lfr5] = rateCalc(spkOI, srData);
    [unitData.psd_halves_freq, unitData.psd_halves, unitData.freq, unitData.psd, unitData.psd_conf,...
      unitData.psd_numelSignal, unitData.beta_halves, unitData.beta] = psdCalc(spkOI, srData, optCoh);
    freqPSD = unitData.freq;
    
    if sum(PR)
      [unitData.freq, unitData.coh, unitData.phase, unitData.coh_conf, unitData.phase_confU, unitData.phase_confL,...
        unitData.coh_halves_freq, unitData.coh_halves, unitData.coh_conf_halves, unitData.phase_halves,...
        unitData.phase_conf_halves] = phaseCohCalc(PR, spkOI, srData, optCoh);
      [unitData.mfr_pop, unitData.mfr_1sthalf_pop, unitData.mfr_2ndhalf_pop, unitData.lfr1_pop,...
        unitData.lfr5_pop] = rateCalc(PR, srData);
      
      % Compute the rate adjustment factors for coherence report, per Aoi et al.
      % we adjust as if the neuron's rate is 1 spk/s
      [unitData.rateadjust_kappa, unitData.rateadjust_kappa_halves] = kappaCalc(unitData.mfr, unitData.mfr_1sthalf,...
        unitData.mfr_2ndhalf, unitData.psd, unitData.psd_halves, unitData.coh, unitData.coh_halves);
      [unitData.psd_halves_freq, unitData.psd_halves, unitData.freq, unitData.psd, unitData.psd_conf,...
        unitData.coh, unitData.phase, unitData.coh_conf, unitData.phase_confU, unitData.phase_confL,...
        unitData.coh_halves_freq, unitData.coh_halves, unitData.coh_conf_halves, unitData.phase_halves,...
        unitData.phase_conf_halves, unitData.rateadjust_kappa, unitData.rateadjust_kappa_halves] = correctFreq(...
        unitData.psd_halves_freq, unitData.psd_halves, freqPSD, unitData.psd, unitData.psd_conf,...
        unitData.freq, unitData.coh, unitData.phase, unitData.coh_conf, unitData.phase_confU,...
        unitData.phase_confL, unitData.coh_halves_freq, unitData.coh_halves, unitData.coh_conf_halves,...
        unitData.phase_halves, unitData.phase_conf_halves, unitData.rateadjust_kappa, unitData.rateadjust_kappa_halves);
      
      % stPR with the local (same shank/tetrode) population rate
      window = 1e2;
      unitData.stPRsh = stprCalc(spkOI, PR, window);
      unitData.stPR = stprHalfCalc(spkOI, PR, window);
      
      % Mean vs var using different bin sizes
      [unitData.m, unitData.mv1, unitData.mv2, unitData.mean_var, unitData.mean_var_timeBin] = meanVarCalc(spkOI, srData);
      
      % Extract phase and coherence for every FOI
      [unitData.phaseFOI, unitData.cohFOI, unitData.coh_confFOI, unitData.actualFOI, unitData.fInds] = phaseCohFOI(FOI,...
        unitData.freq, unitData.phase, unitData.coh, unitData.coh_conf, unitData.rateadjust_kappa);
    else
      unitData.phaseFOI = nan(size(FOI));
      unitData.cohFOI = nan(size(FOI));
      unitData.coh_confFOI = nan(size(FOI));
      unitData.actualFOI = nan(size(FOI));
      unitData.fInds = nan(size(FOI));
    end
    
    % Save unit summary data
    saveParfor([entryName '_' unitNamer(u) '_phaseCoh.mat'], unitData);
    fprintf('Finished processing unit %i\n',unitData.unit);
    
  end % loop over units
  
  % Save shank summary data
  fileList = dir([entryName '_*_phaseCoh.mat']);
  if ~isempty(fileList)
    phaseCohDeep = struct([]);
    phaseCohSuperficial = struct([]);
    for i = 1:size(fileList,1)
      qPhaseCoh = load(fileList(i).name);
      unitInd = unitMetadata(:,1) == qPhaseCoh.q.unit;
      if strcmpi(repository, 'allensdk')
        unitCh2 = unitMetadata(unitInd,9);
      else
        unitCh2 = unitMetadata(unitInd,3);
      end
      if unitCh2 <= divCh
        phaseCohDeep{end+1} = qPhaseCoh.q; %#ok<SAGROW>
      else
        phaseCohSuperficial{end+1} = qPhaseCoh.q; %#ok<SAGROW>
      end
    end
    dataString = ['dataStruct.seriesDataLaminar.' entryName '.shankData.' shankIDs{1} '.phaseCohDeep = phaseCohDeep;'];
    eval(dataString);
    dataString = ['dataStruct.seriesDataLaminar.' entryName '.shankData.' shankIDs{1} '.phaseCohSuperficial = phaseCohSuperficial;'];
    eval(dataString);
    if intermediateSaving
      save(dataFile,'dataStruct','-v7.3'); %#ok<*UNRCH>
    end
    delete([entryName '_*.mat']);
  end
  
  fprintf('Finished processing db entry %i\n',dbCount);
end % loop over db entries


%% SAVE DATA
if ~intermediateSaving
  save(dataFile,'dataStruct','-v7.3');
end
clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca dataStruct