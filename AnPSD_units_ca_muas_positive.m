% A part of the old AnPSD script adapted to perform coherence and phase
% analyses for MUA vs population spiking rate of a different brain area
% (within positive subpopulation only).
% MUAsAll refers to the population rate of area 1
% PR refers to the population rate of area 2


%% LOAD PRE-PROCESSED DATA
if ~exist('dataStruct','var')
  load(dataFile);
end


%% INITILIASE PARAMETERS
params
lists
overwrite = true; % Overwrites existig data. Otherwise if false, keeps the old data.
intermediateSaving = false; % Save everytime a db series is finished being analysed (bad idea since saving is a lengthy process)


%% ADJUST THE NUMBER OF PARALLEL PROCESSORS (applies only when running on NSG Portal)
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    poolsize = 0;
else
    poolsize = p.NumWorkers;
end
if poolsize == 20
  delete(gcp('nocreate'))
  parpool(24);
end


%% DETERMINE SERIES TO COMPARE
series2compare


%% PERFORM COHERENCE ANALSYSES
animal = dataStruct.seriesData.(fnsData{1}).db(1).animal;
if ~isempty(dbEntries_c) && dbEntries_c(1) == inf % Here you can chose to execute only certain db entries
  dbEntries_cLocal = 1:numel(series_c);
else
  dbEntries_cLocal = dbEntries_c;
end
for dbCount = dbEntries_cLocal % Loop through db entries
  seriesName1 = series_c{dbCount};
  if ~isfield(dataStruct, 'seriesData_positive') || ~isfield(dataStruct.seriesData_positive, [animal '_s' seriesName1])
    continue
  else
    fnsData_positive = fieldnames(dataStruct.seriesData_positive);
  end
  
  % Load the contents of dbStruct
  seriesInd = find(endsWith(fnsData_positive,seriesName1));
  [dbStruct, ~, ~, entryName, ~, ~, shankIDs,...
    ~, period, ~, ~, srData, optCoh, ~, FOI,...
    MUAsAll, spk, muas, ~, ~, muaMetadata] = get_dbStruct(dataStruct, seriesInd, 'positive');
  optCoh.maxFreq = maxFreq_ca;
  optCoh.winfactor = winfactor;
  optCoh.freqfactor = freqfactor;
  optCoh.monotoneFreq = true;
  
  % Indices for truncating data and the population rate for the first series
  if ~sum(MUAsAll)
    continue
  end
  [seriesName, animal] = seriesFromEntry(entryName);
  runningSpeedDataEntry = [animal '_s' seriesName(1:min([14 numel(seriesName)]))];
  if excludeRunning && isfield(dataStruct, 'runningSpeedData') && ~isempty(dataStruct.runningSpeedData.(runningSpeedDataEntry).maxQuietPeriod)
    period = combinePeriods(period, dataStruct.runningSpeedData.(runningSpeedDataEntry).maxQuietPeriod, srData);
  end
  [inds1, MUAsAll] = determineInds(period, srData, MUAsAll);
  spk = spk(:,inds1);
  
  for sca = 1:numel(series_ca{dbCount}) % Loop through comparison areas
    sca_name = series_ca{dbCount}{sca};
    phaseCohDB.phaseFOI = []; 
    phaseCohDB.cohFOI = [];
    phaseCohDB.coh_confFOI = [];
    
    % Test if the analysis data already exists
    if ~overwrite && isfield(dataStruct, 'seriesData_ca_muas_positive')
      dataString = [animal '_s' num2str(series_c{dbCount}) '__' animal '_s' num2str(sca_name)];
      if isfield(dataStruct.seriesData_ca_muas_positive, dataString)
        continue
      end
    end
    
    dbStruct2 = dataStruct.seriesData_positive.([animal '_s' num2str(sca_name)]);
    if isempty(dbStruct2)
      continue
    end
    if iscell(dbStruct.dbSeries.period)
      assert(iscell(dbStruct2.dbSeries.period))
      for iPeriod = 1:numel(dbStruct.dbSeries.period)
        assert(dbStruct.dbSeries.period{iPeriod}(1) == dbStruct2.dbSeries.period{iPeriod}(1));
        assert(dbStruct.dbSeries.period{iPeriod}(2) == dbStruct2.dbSeries.period{iPeriod}(2));
      end
    else
      assert(dbStruct.dbSeries.period(1) == dbStruct2.dbSeries.period(1));
      assert(dbStruct.dbSeries.period(2) == dbStruct2.dbSeries.period(2));
    end
    PR = sum(dbStruct2.popData.MUAsAll,1); % PR for the second series
    if excludeRunning && isfield(dataStruct, 'runningSpeedData') && ~isempty(dataStruct.runningSpeedData.(runningSpeedDataEntry).maxQuietPeriod)
      period = combinePeriods(dbStruct2.dbSeries.period, dataStruct.runningSpeedData.(runningSpeedDataEntry).maxQuietPeriod, srData);
    else
      period = dbStruct2.dbSeries.period;
    end
    [inds2, PR] = determineInds(period, srData, PR);
    
    if isempty(spk) || isempty(PR)
      continue
    elseif size(spk,2) > numel(PR)
      PR = [PR zeros(1, size(spk,2) - numel(PR))]; %#ok<*AGROW>
    elseif size(spk,2) < numel(PR)
      PR = PR(1:size(spk,2));
    end
    
    dbStruct_ca = dataStruct.seriesData_ca_positive.([animal '_s' num2str(series_c{dbCount}) '__' animal '_s' num2str(sca_name)]);
    
    if ~isempty(muas)
      parfor u = 1:numel(muas) % Loop through muas
      %for u = 1:numel(muas)
        unit = muas(u);
        fprintf('Started processing unit %i\n',unit);
        
        unitData = [];
        unitData_ca = [];
        spkOI = full(spk(u,:));
        if isempty(spkOI)
          continue
        end
        
        if ~isempty(dbStruct_ca)
          for sh = 1:numel(shankIDs)
            units = dbStruct_ca.shankData.(['shank' num2str(sh)]).units;
            if ismember(unit, units)
              unitData_ca = dbStruct_ca.shankData.(['shank' num2str(sh)]).phaseCoh{unit == units};
              break
            end
          end
        end
        
        if isempty(unitData_ca)
          for sh = 1:numel(shankIDs)
            shankStruct = get_shankStruct(dbStruct, sh);
            if isempty(shankStruct)
              continue
            else
              units = shankStruct.units;
              if isfield(shankStruct, 'phaseCoh') && ismember(unit, units)
                unitData = shankStruct.phaseCoh{unit == units};
                break
              end
            end
          end
          if isempty(unitData)
            [unitData.mfr, unitData.mfr_1sthalf, unitData.mfr_2ndhalf] = rateCalc(full(spkOI), srData);
            [unitData.psd_halves_freq, unitData.psd_halves, freqPSD, unitData.psd, unitData.psd_conf] = psdCalc(full(spkOI), srData, optCoh);
          else
            freqPSD = unitData.freq;
          end
          
          % Sumarry data
          unitData_ca.unit = unit;
          [unitData_ca.freq, unitData_ca.coh, unitData_ca.phase, unitData_ca.coh_conf, unitData_ca.phase_confU,...
            unitData_ca.phase_confL, unitData_ca.coh_halves_freq, unitData_ca.coh_halves, unitData_ca.coh_conf_halves,...
            unitData_ca.phase_halves, unitData_ca.phase_conf_halves] = phaseCohCalc(PR, spkOI, srData, optCoh);
          
          % Compute the rate adjustment factors for coherence report, per Aoi et al.
          % we adjust as if the neuron's rate is 1 spk/s
          [unitData_ca.rateadjust_kappa, unitData_ca.rateadjust_kappa_halves] = kappaCalc(unitData.mfr, unitData.mfr_1sthalf,...
            unitData.mfr_2ndhalf, unitData.psd, unitData.psd_halves, unitData_ca.coh, unitData_ca.coh_halves);
          [~, ~, unitData_ca.freq, ~, ~, unitData_ca.coh, unitData_ca.phase, unitData_ca.coh_conf, unitData_ca.phase_confU,...
            unitData_ca.phase_confL, unitData_ca.coh_halves_freq, unitData_ca.coh_halves, unitData_ca.coh_conf_halves,...
            unitData_ca.phase_halves, unitData_ca.phase_conf_halves, unitData_ca.rateadjust_kappa, unitData_ca.rateadjust_kappa_halves...
            ] = correctFreq(unitData.psd_halves_freq, unitData.psd_halves, freqPSD, unitData.psd, unitData.psd_conf,...
            unitData_ca.freq, unitData_ca.coh, unitData_ca.phase, unitData_ca.coh_conf, unitData_ca.phase_confU,...
            unitData_ca.phase_confL, unitData_ca.coh_halves_freq, unitData_ca.coh_halves, unitData_ca.coh_conf_halves,...
            unitData_ca.phase_halves, unitData_ca.phase_conf_halves, unitData_ca.rateadjust_kappa, unitData_ca.rateadjust_kappa_halves);
          
          % Spike-trigerred population rate (stPR)
          window = 1e2;
          unitData_ca.stPR = stprCalc(spkOI, PR, window);
          unitData_ca.stPR_halves = stprHalfCalc(spkOI, PR, window);
          
          % Extract phase and coherence for every FOI
          [unitData_ca.phaseFOI, unitData_ca.cohFOI, unitData_ca.coh_confFOI, unitData_ca.actualFOI, unitData_ca.fInds] = phaseCohFOI(FOI,...
            unitData_ca.freq, unitData_ca.phase, unitData_ca.coh, unitData_ca.coh_conf, unitData_ca.rateadjust_kappa);
        end
        
        % Save unit summary data
        saveParfor([entryName '_' unitNamer(u) '_phaseCoh.mat'], unitData_ca);
        fprintf('Finished processing unit %i\n',unit);
        
      end % loop over muas
    end
    
    % Save shank summary data
    fileList = dir([entryName '_*_phaseCoh.mat']);
    if ~isempty(fileList)
      phaseCoh = struct([]);
      for i = 1:size(fileList,1)
        qPhaseCoh = load(fileList(i).name);
        phaseCoh{end+1} = qPhaseCoh.q; %#ok<*SAGROW>
        phaseCohDB.phaseFOI = [phaseCohDB.phaseFOI; phaseCoh{i}.phaseFOI];
        phaseCohDB.cohFOI = [phaseCohDB.cohFOI; phaseCoh{i}.cohFOI];
        phaseCohDB.coh_confFOI = [phaseCohDB.coh_confFOI; phaseCoh{i}.coh_confFOI];
      end
      dataString = ['dataStruct.seriesData_ca_muas_positive.' animal '_s' num2str(series_c{dbCount}) '__' animal '_s' num2str(sca_name) '.phaseCoh = phaseCoh;'];
      eval(dataString);
      dataString = ['dataStruct.seriesData_ca_muas_positive.' animal '_s' num2str(series_c{dbCount}) '__' animal '_s' num2str(sca_name) '.units = muas;'];
      eval(dataString);
      dataString = ['dataStruct.seriesData_ca_muas_positive.' animal '_s' num2str(series_c{dbCount}) '__' animal '_s' num2str(sca_name) '.unitMetadata = muaMetadata;'];
      eval(dataString);
      dataString = ['dataStruct.seriesData_ca_muas_positive.' animal '_s' num2str(series_c{dbCount}) '__' animal '_s' num2str(sca_name) '.popData.phaseCoh = phaseCohDB;'];
      eval(dataString);
      delete([entryName '_*.mat']);
    else
      dataString = ['dataStruct.seriesData_ca_muas_positive.' animal '_s' num2str(series_c{dbCount}) '__' animal '_s' num2str(sca_name) '.phaseCoh = [];'];
      eval(dataString);
      dataString = ['dataStruct.seriesData_ca_muas_positive.' animal '_s' num2str(series_c{dbCount}) '__' animal '_s' num2str(sca_name) '.units = [];'];
      eval(dataString);
      dataString = ['dataStruct.seriesData_ca_muas_positive.' animal '_s' num2str(series_c{dbCount}) '__' animal '_s' num2str(sca_name) '.unitMetadata = [];'];
      eval(dataString);
      dataString = ['dataStruct.seriesData_ca_muas_positive.' animal '_s' num2str(series_c{dbCount}) '__' animal '_s' num2str(sca_name) '.popData.phaseCoh = [];'];
      eval(dataString);
    end
    if intermediateSaving
      save(dataFile,'dataStruct','-v7.3'); %#ok<*UNRCH>
    end
    
    fprintf('Finished comparing against series %i\n',sca);
  end % loop over compared against areas
  fprintf('Finished processing db entry %i\n',dbCount);
end % loop over db entries

% Save data if it hasn't been already saved
if ~intermediateSaving
  save(dataFile,'dataStruct','-v7.3');
end
clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca dataStruct


function popData_ca = comparePops(PR, spkDB, FOI, srData, opt)

if numel(spkDB) > numel(PR)
  PR = [PR zeros(1, numel(spkDB) - numel(PR))];
elseif numel(spkDB) < numel(PR)
  spkDB = [spkDB zeros(1, numel(PR) - numel(spkDB))];
end

[popData_ca.mfr, popData_ca.mfr_1sthalf, popData_ca.mfr_2ndhalf, popData_ca.lfr1, popData_ca.lfr5] = rateCalc(spkDB, srData);
[popData_ca.psd_halves_freq, popData_ca.psd_halves, freqPSD, popData_ca.psd, popData_ca.psd_conf,...
  popData_ca.psd_numelSignal, popData_ca.beta_halves, popData_ca.beta] = psdCalc(spkDB, srData, opt);

[popData_ca.freq, popData_ca.coh, popData_ca.phase, popData_ca.coh_conf, popData_ca.phase_confU,...
  popData_ca.phase_confL, popData_ca.coh_halves_freq, popData_ca.coh_halves, popData_ca.coh_conf_halves,...
  popData_ca.phase_halves, popData_ca.phase_conf_halves] = phaseCohCalc(PR, spkDB, srData, opt);

% Compute the rate adjustment factors for coherence report, per Aoi et al.
% we adjust as if the neuron's rate is 1 spk/s
[popData_ca.rateadjust_kappa, popData_ca.rateadjust_kappa_halves] = kappaCalc(popData_ca.mfr, popData_ca.mfr_1sthalf,...
  popData_ca.mfr_2ndhalf, popData_ca.psd, popData_ca.psd_halves, popData_ca.coh, popData_ca.coh_halves);
[popData_ca.psd_halves_freq, popData_ca.psd_halves, popData_ca.freq, popData_ca.psd, popData_ca.psd_conf,...
  popData_ca.coh, popData_ca.phase, popData_ca.coh_conf, popData_ca.phase_confU, popData_ca.phase_confL,...
  popData_ca.coh_halves_freq, popData_ca.coh_halves, popData_ca.coh_conf_halves, popData_ca.phase_halves,...
  popData_ca.phase_conf_halves, popData_ca.rateadjust_kappa, popData_ca.rateadjust_kappa_halves] = correctFreq(...
  popData_ca.psd_halves_freq, popData_ca.psd_halves, freqPSD, popData_ca.psd, popData_ca.psd_conf,...
  popData_ca.freq, popData_ca.coh, popData_ca.phase, popData_ca.coh_conf, popData_ca.phase_confU,...
  popData_ca.phase_confL, popData_ca.coh_halves_freq, popData_ca.coh_halves, popData_ca.coh_conf_halves,...
  popData_ca.phase_halves, popData_ca.phase_conf_halves, popData_ca.rateadjust_kappa, popData_ca.rateadjust_kappa_halves);

% Spike-trigerred population rate (stPR)
window = 1e2;
popData_ca.stPR = stprCalc(spkDB, PR, window);
popData_ca.stPR_halves = stprHalfCalc(spkDB, PR, window);

% Extract phase and coherence for every FOI
[popData_ca.phaseFOI, popData_ca.cohFOI, popData_ca.coh_confFOI, popData_ca.actualFOI, popData_ca.fInds] = phaseCohFOI(FOI,...
  popData_ca.freq, popData_ca.phase, popData_ca.coh, popData_ca.coh_conf, popData_ca.rateadjust_kappa);
end