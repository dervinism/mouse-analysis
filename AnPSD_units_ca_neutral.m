% A part of the old AnPSD script adapted to perform coherence and phase
% analyses for unit vs population spiking rate of a different brain area
% (within uncorrelated subpopulation only).


%% LOAD PRE-PROCESSED DATA
if ~exist('dataStruct','var')
  load(dataFile);
end


%% INITILIASE PARAMETERS
params
lists
overwrite = true;
intermediateSaving = false;


%% ADJUST THE NUMBER OF PARALLEL PROCESSORS
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
if ~isempty(dbEntries_c) && dbEntries_c(1) == inf
  dbEntries_cLocal = 1:numel(series_c);
else
  dbEntries_cLocal = dbEntries_c;
end
for dbCount = dbEntries_cLocal % Loop through db entries
  seriesName1 = series_c{dbCount};
  if ~isfield(dataStruct, 'seriesData_neutral') || ~isfield(dataStruct.seriesData_neutral, [animal '_s' seriesName1])
    continue
  else
    fnsData_neutral = fieldnames(dataStruct.seriesData_neutral);
  end
  
  % Load the contents of dbStruct
  seriesInd = find(endsWith(fnsData_neutral,seriesName1));
  [dbStruct, ~, ~, entryName, ~, ~, shankIDs,...
    ~, period, ~, ~, srData, optCoh, ~, FOI,...
    MUAsAll] = get_dbStruct(dataStruct, seriesInd, 'neutral');
  optCoh.winfactor = winfactor;
  optCoh.freqfactor = freqfactor;
  optCoh.monotoneFreq = true;
  if isempty(dbStruct) || ~isfield(dbStruct, 'shankData') || isempty(dbStruct.shankData) || ~sum(MUAsAll)
    continue
  end
  
  % Indices for truncating data and the population rate for the first series
  [inds1, spkDB] = determineInds(period, srData, MUAsAll);
  
  % Obtain population rate constructed from pure units only for the first series
    spkDB2 = [];
    for sh = 1:numel(shankIDs)
      [shankStruct, shank, units, unitMetadata, xcoords, ycoords, spk,...
        MUAs, phaseCoh] = get_shankStruct(dbStruct, sh);
      if ~isempty(spk)
        if size(spk,2) < size(MUAsAll,2)
          spk = [full(spk) zeros(size(spk,1), size(MUAsAll,2)-size(spk,2))];
        end
        spk = spk(:,inds1);
      end
      for u = 1:numel(units)
        spkOI = full(spk(u,:));
        if numel(spkDB2) > numel(spkOI)
          spkDB2 = spkDB2 + [spkOI zeros(1,numel(spkDB2)-numel(spkOI))]; % pure units PR for the first series
        elseif numel(spkDB2) < numel(spkOI)
          spkDB2 = [spkDB2 zeros(1,numel(spkOI)-numel(spkDB2))] + spkOI;
        else
          spkDB2 = spkDB2 + spkOI;
        end
      end
    end
  
  for sca = 1:numel(series_ca{dbCount}) % Loop through comparison areas
    sca_name = series_ca{dbCount}{sca};
    phaseCohDB.phaseFOI = []; 
    phaseCohDB.cohFOI = [];
    phaseCohDB.coh_confFOI = [];
    
  % Test if the analysis data already exists
    if ~overwrite && isfield(dataStruct, 'seriesData_ca_neutral')
      dataString = [animal '_s' num2str(series_c{dbCount}) '__' animal '_s' num2str(sca_name)];
      if isfield(dataStruct.seriesData_ca_neutral, dataString)
        if isfield(dataStruct.seriesData_ca_neutral.(dataString), 'shankData')
          continue
        end
      end
    end
    
    if ~isfield(dataStruct.seriesData_neutral, [animal '_s' num2str(sca_name)])
      continue
    end
    dbStruct_ca = dataStruct.seriesData_neutral.([animal '_s' num2str(sca_name)]);
    if isempty(dbStruct_ca) || ~isfield(dbStruct_ca, 'shankData') || isempty(dbStruct_ca.shankData)
      continue
    end
    shankIDs_ca = fieldnames(dbStruct_ca.shankData);
    if iscell(dbStruct.dbSeries.period)
      assert(iscell(dbStruct_ca.dbSeries.period))
      for iPeriod = 1:numel(dbStruct.dbSeries.period)
        assert(dbStruct.dbSeries.period{iPeriod}(1) == dbStruct_ca.dbSeries.period{iPeriod}(1));
        assert(dbStruct.dbSeries.period{iPeriod}(2) == dbStruct_ca.dbSeries.period{iPeriod}(2));
      end
    else
      assert(dbStruct.dbSeries.period(1) == dbStruct_ca.dbSeries.period(1));
      assert(dbStruct.dbSeries.period(2) == dbStruct_ca.dbSeries.period(2));
    end
    PR = sum(dbStruct_ca.popData.MUAsAll,1); % PR for the second series
    [inds2, PR] = determineInds(dbStruct_ca.dbSeries.period, srData, PR);
    
    % Obtain population rate constructed from pure units only for the second series
    PR2 = [];
    for sh = 1:numel(shankIDs_ca)
      shankStruct = eval(['dbStruct_ca.shankData.' shankIDs_ca{sh}]);
      if isempty(shankStruct)
        continue
      end
      units = shankStruct.units;
      if ~isempty(shankStruct.spk)
        if size(shankStruct.spk,2) < numel(PR)
          spk_ca = [full(shankStruct.spk) zeros(size(shankStruct.spk,1), numel(PR)-size(shankStruct.spk,2))];
          spk_ca = spk_ca(:,inds2);
        else
          spk_ca = shankStruct.spk(:,inds2);
        end
      else
        spk_ca = shankStruct.spk;
      end
      for u = 1:numel(units)
        spkOI_ca = full(spk_ca(u,:));
        if numel(PR2) > numel(spkOI_ca)
          PR2 = PR2 + [spkOI_ca zeros(1,numel(PR2)-numel(spkOI_ca))]; % pure units PR for the sacond series
        elseif numel(PR2) < numel(spkOI_ca)
          PR2 = [PR2 zeros(1,numel(spkOI_ca)-numel(PR2))] + spkOI_ca;
        else
          PR2 = PR2 + spkOI_ca;
        end
      end
    end
    
    if ~(isempty(PR2) || isempty(spkDB2))
      if numel(PR2) > numel(spkDB2)
        PR2 = PR2(1:numel(spkDB2));
      elseif numel(PR2) < numel(spkDB2)
        PR2 = [PR2 zeros(1,numel(spkDB2)-numel(PR2))];
      end
    end
    
    for sh = 1:numel(shankIDs) % Loop through shanks
      fprintf('db %d %s v %s shank %d -------------------\n', dbCount, ['s' num2str(series_c{dbCount})], ['s' num2str(sca_name)], sh);
      
      % Load the contents of shankStruct
      [shankStruct, shank, units, unitMetadata, xcoords, ycoords, spk,...
        ~, phaseCoh] = get_shankStruct(dbStruct, sh);
      if ~isempty(spk)
        if size(spk,2) < size(MUAsAll,2)
          spk = [full(spk) zeros(size(spk,1), size(MUAsAll,2)-size(spk,2))];
        end
        spk = spk(:,inds1);
      end
      
      if isempty(spkDB) || isempty(PR)
        continue
      elseif size(spk,2) > numel(PR)
      	PR = [PR zeros(1, size(spk,2) - numel(PR))]; %#ok<*AGROW>
      elseif size(spk,2) < numel(PR)
        PR = PR(1:size(spk,2));
      end
      
      if ~isempty(units)
        parfor u = 1:numel(units) % Loop through units
        %for u = 1:numel(units)
          unit = units(u);
          fprintf('Started processing unit %i\n',unit);
          
          unitData_ca = struct();
          unitData_ca.unit = unit;
          spkOI = full(spk(u,:));
          if isempty(spkOI)
            continue
          end
          
          if isfield(shankStruct, 'phaseCoh')
            unitData = shankStruct.phaseCoh{u};
          else
            [unitData.mfr, unitData.mfr_1sthalf, unitData.mfr_2ndhalf] = rateCalc(full(spkOI), srData);
            [~, unitData.psd_halves, ~, unitData.psd] = psdCalc(full(spkOI), srData, optCoh);
          end

          % Sumarry data
          [unitData_ca.freq, unitData_ca.coh, unitData_ca.phase, unitData_ca.coh_conf, unitData_ca.phase_confU,...
            unitData_ca.phase_confL, unitData_ca.coh_halves_freq, unitData_ca.coh_halves, unitData_ca.coh_conf_halves,...
            unitData_ca.phase_halves, unitData_ca.phase_conf_halves] = phaseCohCalc(PR, spkOI, srData, optCoh);

          % Compute the rate adjustment factors for coherence report, per Aoi et al.
          % we adjust as if the neuron's rate is 1 spk/s
          [unitData_ca.rateadjust_kappa, unitData_ca.rateadjust_kappa_halves] = kappaCalc(unitData.mfr, unitData.mfr_1sthalf,...
            unitData.mfr_2ndhalf, unitData.psd, unitData.psd_halves, unitData_ca.coh, unitData_ca.coh_halves);
          if numel(unitData_ca.coh) == 1 && isnan(unitData_ca.coh)
            unitData_ca.rateadjust_kappa = 1;
            unitData_ca.rateadjust_kappa_halves = [1;1];
          else
            assert(size(unitData_ca.coh,2) == size(unitData_ca.rateadjust_kappa,2));
            if size(unitData_ca.coh_halves,2) > 1
              assert(size(unitData_ca.coh_halves,2) == size(unitData_ca.rateadjust_kappa_halves,2));
            elseif sum(isnan(unitData_ca.coh_halves)) == 2
              unitData_ca.rateadjust_kappa_halves = [NaN; NaN];
            end
          end

          % Spike-trigerred population rate (stPR)
          window = 1e2;
          unitData_ca.stPR = stprCalc(spkOI, PR, window);
          unitData_ca.stPR_halves = stprHalfCalc(spkOI, PR, window);

          % Extract phase and coherence for every FOI
          [unitData_ca.phaseFOI, unitData_ca.cohFOI, unitData_ca.coh_confFOI, unitData_ca.actualFOI, unitData_ca.fInds] = phaseCohFOI(FOI,...
            unitData_ca.freq, unitData_ca.phase, unitData_ca.coh, unitData_ca.coh_conf, unitData_ca.rateadjust_kappa);

          % Save unit summary data
          saveParfor([entryName '_' unitNamer(u) '_phaseCoh.mat'], unitData_ca);
          fprintf('Finished processing unit %i\n',unit);

        end % loop over units
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
        dataString = ['dataStruct.seriesData_ca_neutral.' animal '_s' num2str(series_c{dbCount}) '__' animal '_s' num2str(sca_name) '.shankData.' shankIDs{sh} '.phaseCoh = phaseCoh;'];
        eval(dataString);
        dataString = ['dataStruct.seriesData_ca_neutral.' animal '_s' num2str(series_c{dbCount}) '__' animal '_s' num2str(sca_name) '.shankData.' shankIDs{sh} '.shankID = sh;'];
        eval(dataString);
        dataString = ['dataStruct.seriesData_ca_neutral.' animal '_s' num2str(series_c{dbCount}) '__' animal '_s' num2str(sca_name) '.shankData.' shankIDs{sh} '.units = units;'];
        eval(dataString);
        dataString = ['dataStruct.seriesData_ca_neutral.' animal '_s' num2str(series_c{dbCount}) '__' animal '_s' num2str(sca_name) '.shankData.' shankIDs{sh} '.unitMetadata = unitMetadata;'];
        eval(dataString);
        dataString = ['dataStruct.seriesData_ca_neutral.' animal '_s' num2str(series_c{dbCount}) '__' animal '_s' num2str(sca_name) '.popData.phaseCoh = phaseCohDB;'];
        eval(dataString);
        if intermediateSaving
          save(dataFile,'dataStruct','-v7.3'); %#ok<*UNRCH>
        end
        delete([entryName '_*.mat']);
      else
        dataString = ['dataStruct.seriesData_ca_neutral.' animal '_s' num2str(series_c{dbCount}) '__' animal '_s' num2str(sca_name) '.shankData.' shankIDs{sh} '.phaseCoh = {};'];
        eval(dataString);
        dataString = ['dataStruct.seriesData_ca_neutral.' animal '_s' num2str(series_c{dbCount}) '__' animal '_s' num2str(sca_name) '.shankData.' shankIDs{sh} '.shankID = sh;'];
        eval(dataString);
        dataString = ['dataStruct.seriesData_ca_neutral.' animal '_s' num2str(series_c{dbCount}) '__' animal '_s' num2str(sca_name) '.shankData.' shankIDs{sh} '.units = [];'];
        eval(dataString);
        dataString = ['dataStruct.seriesData_ca_neutral.' animal '_s' num2str(series_c{dbCount}) '__' animal '_s' num2str(sca_name) '.shankData.' shankIDs{sh} '.unitMetadata = [];'];
        eval(dataString);
        dataString = ['dataStruct.seriesData_ca_neutral.' animal '_s' num2str(series_c{dbCount}) '__' animal '_s' num2str(sca_name) '.popData.phaseCoh = phaseCohDB;'];
        eval(dataString);
        if intermediateSaving
          save(dataFile,'dataStruct','-v7.3'); %#ok<*UNRCH>
        end
      end
      
      fprintf('Finished processing shank %i\n',sh);
    end % loop over shanks
    
    % Compare populations
    if ~(isempty(PR) || isempty(spkDB))
      popData_ca = comparePops(PR, spkDB, FOI, srData, optCoh);
    else
      popData_ca = [];
    end
    if ~(isempty(PR2) || isempty(spkDB2))
      popData_ca2 = comparePops(PR2, spkDB2, FOI, srData, optCoh);
    else
      popData_ca2 = [];
    end
    
    % Save
    dataString = ['dataStruct.seriesData_ca_neutral.' animal '_s' num2str(series_c{dbCount}) '__' animal '_s' num2str(sca_name) '.popData.phaseCohPop = popData_ca;'];
    eval(dataString);
    dataString = ['dataStruct.seriesData_ca_neutral.' animal '_s' num2str(series_c{dbCount}) '__' animal '_s' num2str(sca_name) '.popData.phaseCohPopUnitOnly = popData_ca2;'];
    eval(dataString);
    if intermediateSaving
      save(dataFile,'dataStruct','-v7.3');
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
[popData_ca.psd_halves_freq, popData_ca.psd_halves, popData_ca.freq, popData_ca.psd, popData_ca.psd_conf,...
  popData_ca.psd_numelSignal, popData_ca.beta_halves, popData_ca.beta] = psdCalc(spkDB, srData, opt);

[popData_ca.freq, popData_ca.coh, popData_ca.phase, popData_ca.coh_conf, popData_ca.phase_confU,...
  popData_ca.phase_confL, popData_ca.coh_halves_freq, popData_ca.coh_halves, popData_ca.coh_conf_halves,...
  popData_ca.phase_halves, popData_ca.phase_conf_halves] = phaseCohCalc(PR, spkDB, srData, opt);

% Compute the rate adjustment factors for coherence report, per Aoi et al.
% we adjust as if the neuron's rate is 1 spk/s
[popData_ca.rateadjust_kappa, popData_ca.rateadjust_kappa_halves] = kappaCalc(popData_ca.mfr, popData_ca.mfr_1sthalf,...
  popData_ca.mfr_2ndhalf, popData_ca.psd, popData_ca.psd_halves, popData_ca.coh, popData_ca.coh_halves);

% Spike-trigerred population rate (stPR)
window = 1e2;
popData_ca.stPR = stprCalc(spkDB, PR, window);
popData_ca.stPR_halves = stprHalfCalc(spkDB, PR, window);

% Extract phase and coherence for every FOI
[popData_ca.phaseFOI, popData_ca.cohFOI, popData_ca.coh_confFOI, popData_ca.actualFOI, popData_ca.fInds] = phaseCohFOI(FOI,...
  popData_ca.freq, popData_ca.phase, popData_ca.coh, popData_ca.coh_conf, popData_ca.rateadjust_kappa);
end