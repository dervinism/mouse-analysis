% A part of the old AnPSD script adapted to perform coherence and phase
% analyses for unit and population spiking rate data.


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
  [dbStruct, ~, ~, entryName, ~, ~, shankIDs,...
    ~, period, ~, ~, srData, optCoh, exclRad, FOI,...
    MUAsAll, spkDB] = get_dbStruct(dataStruct, dbCount);
  optCoh.maxFreq = maxFreq;
  optCoh.winfactor = winfactor;
  optCoh.freqfactor = freqfactor;
  optCoh.monotoneFreq = true;
  if ~isempty(dbStruct.popData.muaMetadata)
    xcoordsMUA = dbStruct.popData.muaMetadata(:,4);
    ycoordsMUA = dbStruct.popData.muaMetadata(:,5);
  end
  
  % Indices for truncating data
  [seriesName, animal] = seriesFromEntry(entryName);
  runningSpeedDataEntry = [animal '_s' seriesName(1:min([14 numel(seriesName)]))];
  if excludeRunning && isfield(dataStruct, 'runningSpeedData') && ~isempty(dataStruct.runningSpeedData.(runningSpeedDataEntry).maxQuietPeriod)
    period = combinePeriods(period, dataStruct.runningSpeedData.(runningSpeedDataEntry).maxQuietPeriod, srData);
  end
  [inds, spkDB] = determineInds(period, srData, spkDB);
  
  % Initialise storage variables
  phaseCohDB.mfr_1sthalf = [];
  phaseCohDB.mfr_2ndhalf = [];
  phaseCohDB.beta = [];
  phaseCohDB.phaseFOI = []; 
  phaseCohDB.cohFOI = [];
  phaseCohDB.coh_confFOI = [];

  for sh = 1:numel(shankIDs) % Loop through shanks
    fprintf('%s shank %d -------------------\n', fnsData{dbCount}, sh);
    
    % Load the contents of shankStruct
    [~, shank, units, ~, xcoords, ycoords, spk, MUAs] = get_shankStruct(dbStruct, sh);
    if isempty(units)
      dataString = ['dataStruct.seriesData.' entryName '.popData.phaseCoh = phaseCohDB;'];
      eval(dataString);
      continue
    end
    
    % Truncate shank data
    if ~isempty(spk)
      if size(spk,2) < size(MUAsAll,2)
        spk = [full(spk) zeros(size(spk,1), size(MUAsAll,2)-size(spk,2))];
      end
      spk = spk(:,inds);
    end
    if ~isempty(MUAs)
      if size(MUAs,2) < size(MUAsAll,2)
        MUAs = [MUAs zeros(size(MUAs,1), size(MUAsAll,2)-size(MUAs,2))]; %#ok<*AGROW>
      end
      MUAs = MUAs(:,inds);
    end

    parfor u = 1:numel(units) % Loop through units
    %for u = 1:numel(units)
      fprintf('Started processing unit %i\n',units(u));
      
      unitData = struct();
      unitData.unit = units(u);
      spkOI = full(spk(u,:));
      
      % Spikes for the vicinity of unit of interest (in this case the whole shank)
      if numel(shankIDs) > 1
        localMUA = []; % this is not a single-shank probe
      else
        localMUA = MUAs;
      end
      assert(all(spkOI <= MUAs))
      
      % Population rate
      unit = units(u);
      iCoords = [xcoords(u) ycoords(u)];
      distances = sqrt((xcoordsMUA-iCoords(1)).^2 + (ycoordsMUA-iCoords(2)).^2); %#ok<*PFBNS>
      unitsOI = distances;
      unitsOI(isnan(unitsOI)) = 0;
      unitsOI(unitsOI <= exclRad) = 0;
      unitsOI(unitsOI > exclRad) = 1;
      unitsOI = logical(unitsOI);
      PR = full(sum(spkDB(unitsOI,:),1));
      %PR = MUAs-spkOI;

      % Summary data
      [unitData.mfr, unitData.mfr_1sthalf, unitData.mfr_2ndhalf, unitData.lfr1, unitData.lfr5] = rateCalc(spkOI, srData);
      [unitData.psd_halves_freq, unitData.psd_halves, unitData.freq, unitData.psd, unitData.psd_conf,...
        unitData.psd_numelSignal, unitData.beta_halves, unitData.beta] = psdCalc(spkOI, srData, optCoh);
      freqPSD = unitData.freq;
      
      if sum(PR)
        if numel(shank) == 1 && ~isempty(localMUA) % all we've got is one shank
          [unitData.freq, unitData.coh, unitData.phase, unitData.coh_conf, unitData.phase_confU, unitData.phase_confL,...
            unitData.coh_halves_freq, unitData.coh_halves, unitData.coh_conf_halves, unitData.phase_halves,...
            unitData.phase_conf_halves] = phaseCohCalc(PR, spkOI, srData, optCoh);
          [unitData.mfr_pop, unitData.mfr_1sthalf_pop, unitData.mfr_2ndhalf_pop, unitData.lfr1_pop,...
            unitData.lfr5_pop] = rateCalc(PR, srData);
          [unitData.mfr_mua, unitData.mfr_1sthalf_mua, unitData.mfr_2ndhalf_mua, unitData.lfr1_mua,...
            unitData.lfr5_mua] = rateCalc(localMUA, srData);
          
        else % we have more than 1 shank
          [unitData.freq_local, unitData.coh_local, unitData.phase_local, unitData.coh_conf_local, unitData.phase_confU_local,...
            unitData.phase_confL_local, unitData.coh_halves_freq_local, unitData.coh_halves_local,...
            unitData.coh_conf_halves_local, unitData.phase_halves_local, unitData.phase_conf_halves_local] = phaseCohCalc(PR,...
            spkOI, srData, optCoh);
          [unitData.rateadjust_kappa_local, unitData.rateadjust_kappa_halves_local] = kappaCalc(unitData.mfr,...
            unitData.mfr_1sthalf, unitData.mfr_2ndhalf, unitData.psd, unitData.psd_halves, unitData.coh_local,...
            unitData.coh_halves_local);
          [unitData.mfr_pop_local, unitData.mfr_1sthalf_pop_local, unitData.mfr_2ndhalf_pop_local, unitData.lfr1_pop_local,...
            unitData.lfr5_pop_local] = rateCalc(PR, srData);
          
          PR_global = sum(MUAsAll,1)-MUAs;
          [unitData.freq, unitData.coh, unitData.phase, unitData.coh_conf, unitData.phase_confU, unitData.phase_confL,...
            unitData.coh_halves_freq, unitData.coh_halves, unitData.coh_conf_halves, unitData.phase_halves,...
            unitData.phase_conf_halves] = phaseCohCalc(PR_global, spkOI, srData, optCoh);
          [unitData.mfr_pop, unitData.mfr_1sthalf_pop, unitData.mfr_2ndhalf_pop, unitData.lfr1_pop,...
            unitData.lfr5_pop] = rateCalc(PR_global, srData);
        end
        
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
        
        % Spike-trigerred population rate (stPR)
        window = 1e2;
        if ~(numel(shank) == 1 && ~isempty(localMUA))
          unitData.stPRglobal = stprCalc(spkOI, PR_global, window);
        end
        % stPR with the local (same shank/tetrode) population rate
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
      fprintf('Finished processing unit %i\n',unit);
      
    end % loop over units
    
    % Save shank summary data
    fileList = dir([entryName '_*_phaseCoh.mat']);
    if ~isempty(fileList)
      phaseCoh = struct([]);
      for i = 1:size(fileList,1)
        qPhaseCoh = load(fileList(i).name);
        phaseCoh{end+1} = qPhaseCoh.q; %#ok<SAGROW>
        phaseCohDB.mfr_1sthalf = [phaseCohDB.mfr_1sthalf; phaseCoh{i}.mfr_1sthalf];
        phaseCohDB.mfr_2ndhalf = [phaseCohDB.mfr_2ndhalf; phaseCoh{i}.mfr_2ndhalf];
        phaseCohDB.beta = [phaseCohDB.beta; phaseCoh{i}.beta];
        if numel(phaseCoh{i}.phaseFOI) == 1 && isnan(phaseCoh{i}.phaseFOI)
          phaseCohDB.phaseFOI = [phaseCohDB.phaseFOI; nan(1,size(phaseCohDB.phaseFOI,2))];
          phaseCohDB.cohFOI = [phaseCohDB.cohFOI; nan(1,size(phaseCohDB.cohFOI,2))];
          phaseCohDB.coh_confFOI = [phaseCohDB.coh_confFOI; nan(1,size(phaseCohDB.coh_confFOI,2))];
        else
          phaseCohDB.phaseFOI = [phaseCohDB.phaseFOI; phaseCoh{i}.phaseFOI];
          phaseCohDB.cohFOI = [phaseCohDB.cohFOI; phaseCoh{i}.cohFOI];
          phaseCohDB.coh_confFOI = [phaseCohDB.coh_confFOI; phaseCoh{i}.coh_confFOI];
        end
      end
      dataString = ['dataStruct.seriesData.' entryName '.shankData.' shankIDs{sh} '.phaseCoh = phaseCoh;'];
      eval(dataString);
      dataString = ['dataStruct.seriesData.' entryName '.popData.phaseCoh = phaseCohDB;'];
      eval(dataString);
      if intermediateSaving
        save(dataFile,'dataStruct','-v7.3'); %#ok<*UNRCH>
      end
      delete([entryName '_*.mat']);
    end
    
    fprintf('Finished processing shank %i\n',sh);
  end % loop over shanks
  
  fprintf('Finished processing db entry %i\n',dbCount);
end % loop over db entries


%% SAVE DATA
if ~intermediateSaving
  save(dataFile,'dataStruct','-v7.3');
end
clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca dataStruct