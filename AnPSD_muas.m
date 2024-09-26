% A part of the old AnPSD script adapted to perform coherence and phase
% analyses for mua data.


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
    ~, period, ~, ~, srData, optCoh, exclRad, ~,...
    MUAsAll, spkDB, spkDB_units] = get_dbStruct(dataStruct, dbCount);
  optCoh.maxFreq = maxFreq;
  optCoh.winfactor = winfactor;
  optCoh.freqfactor = freqfactor;
  optCoh.monotoneFreq = true;
  MUAsAll = sum(MUAsAll,1);
  if isempty(MUAsAll) || ~sum(sum(MUAsAll))
    disp(['No spiking data for ' fnsData{dbCount} '. Skippig to the next db entry...']);
    continue
  end
  xcoordsMUA = dbStruct.popData.muaMetadata(:,4);
  ycoordsMUA = dbStruct.popData.muaMetadata(:,5);
  
  % Indices for truncating data
  [seriesName, animal] = seriesFromEntry(entryName);
  runningSpeedDataEntry = [animal '_s' seriesName(1:min([14 numel(seriesName)]))];
  if excludeRunning && isfield(dataStruct, 'runningSpeedData') && ~isempty(dataStruct.runningSpeedData.(runningSpeedDataEntry).maxQuietPeriod)
    period = combinePeriods(period, dataStruct.runningSpeedData.(runningSpeedDataEntry).maxQuietPeriod, srData);
  end
  [~, spkDB] = determineInds(period, srData, spkDB);
  
  % Get unit IDs
  units = [];
  phaseCoh = {};
  for sh = 1:numel(shankIDs) % Loop through shanks
    [~, ~, unitsSh, ~, ~, ~, ~, ~, phaseCohSh] = get_shankStruct(dbStruct, sh);
    units = [units; unitsSh]; %#ok<*AGROW>
    if ~isempty(unitsSh)
      phaseCoh = [phaseCoh phaseCohSh];
    end
  end
  
  if ~isempty(spkDB_units)
    
    parfor u = 1:numel(spkDB_units) % Loop through units
    %for u = 1:numel(spkDB_units)
      fprintf('Started processing unit %i\n',spkDB_units(u));
      
      unit = spkDB_units(u);
      if ismember(spkDB_units(u), units)
        unitInd = units == spkDB_units(u);
        unitData = phaseCoh{unitInd};
      else
        unitData = struct();
        unitData.unit = spkDB_units(u);
        spkOI = full(spkDB(u,:));
        
        % Population rate
        iCoords = [xcoordsMUA(u) ycoordsMUA(u)];
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
          [unitData.freq, unitData.coh, unitData.phase, unitData.coh_conf, unitData.phase_confU, unitData.phase_confL,...
            unitData.coh_halves_freq, unitData.coh_halves, unitData.coh_conf_halves, unitData.phase_halves,...
            unitData.phase_conf_halves] = phaseCohCalc(PR, spkOI, srData, optCoh);
          [unitData.mfr_pop, unitData.mfr_1sthalf_pop, unitData.mfr_2ndhalf_pop, unitData.lfr1_pop,...
            unitData.lfr5_pop] = rateCalc(PR, srData);
          [unitData.mfr_mua, unitData.mfr_1sthalf_mua, unitData.mfr_2ndhalf_mua, unitData.lfr1_mua,...
            unitData.lfr5_mua] = rateCalc(MUAsAll, srData);
          
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
        end
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
      end
      dataString = ['dataStruct.seriesData.' entryName '.popData.phaseCoh.phaseCoh = phaseCoh;'];
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