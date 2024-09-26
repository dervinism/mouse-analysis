% Run this script to perform analyses comparing LFP and MUA of different brain areas.


%% LOAD PRE-PROCESSED DATA
load(dataFile);


%% INITIALISE PARAMETERS
visibility = 'on';
intermediateSaving = false;


%% DETERMINE SERIES TO COMPARE
series2compare


%% CORRELATE LFP WITH CROSS-AREA LFP
figs = {};
animal = dataStruct.seriesData.(fnsData{1}).db(1).animal;
if ~isempty(dbEntries_c) && dbEntries_c(1) == inf
  dbEntries_cLocal = 1:numel(series_c);
else
  dbEntries_cLocal = dbEntries_c;
end
for dbCount = dbEntries_cLocal % Loop through db entries
  
  % Load the contents of dbStruct
  seriesName1 = series_c{dbCount};
  [dbStruct, ~, ~, entryName, ~, ~, ~,...
    ~, ~, ~, ~, srData] = get_dbStruct(dataStruct, find(endsWith(fnsData,seriesName1)));
  
  for sca = 1:numel(series_ca{dbCount})
    sca_name = [entryName '__' animal '_s' series_ca{dbCount}{sca}];
    ind = strfind(sca_name, '__');
    if ~strcmpi(entryName, sca_name(1:ind-1))
      continue
    end
    dbStruct_ca = dataStruct.seriesData.([animal '_s' series_ca{dbCount}{sca}]);
    
    % Loop through LFP channels
    if isfield(dbStruct, 'lfpPowerData')
      chOIDB = dbStruct.lfpPowerData.chOIDB;
    else
      continue
    end
    for iCh = 1:numel(chOIDB)
      
      if isfield(dbStruct_ca, 'lfpPowerData')
        chOIDB_ca = dbStruct_ca.lfpPowerData.chOIDB;
      else
        continue
      end
      lfpTimes = dbStruct_ca.lfpPowerData.lfpTimes;
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
      
      % Load LFP measures for series1
      fastPower1 = dbStruct.lfpPowerData.fastPower{iCh};
      if isempty(fastPower1)
        continue
      end
      
      % Clean NaNs
      fastPower1(isnan(fastPower1)) = 0;
      
      for iCh_ca = 1:numel(chOIDB_ca)
        
        % Load LFP measures
        fastPower2 = dbStruct_ca.lfpPowerData.fastPower{iCh_ca};
        if isempty(fastPower2)
          continue
        end
        
        % Clean NaNs
        fastPower2(isnan(fastPower2)) = 0;
        
        % Calculate phase and coherence measures
        optLFP = dbStruct.conf.optCoh;
        optLFP.winfactor = winfactor;
        optLFP.monotoneFreq = true;
        optLFP.freqfactor = freqfactor;
        optLFP.typespk1 = 'c';
        optLFP.typespk2 = 'c';
        % Calculations
        [lfpPhaseCoh.fastPower{iCh}{iCh_ca}.freq, lfpPhaseCoh.fastPower{iCh}{iCh_ca}.coh, lfpPhaseCoh.fastPower{iCh}{iCh_ca}.phase,...
          lfpPhaseCoh.fastPower{iCh}{iCh_ca}.coh_conf, lfpPhaseCoh.fastPower{iCh}{iCh_ca}.phase_confU, lfpPhaseCoh.fastPower{iCh}{iCh_ca}.phase_confL,...
          lfpPhaseCoh.fastPower{iCh}{iCh_ca}.coh_halves_freq, lfpPhaseCoh.fastPower{iCh}{iCh_ca}.coh_halves, lfpPhaseCoh.fastPower{iCh}{iCh_ca}.coh_conf_halves,...
          lfpPhaseCoh.fastPower{iCh}{iCh_ca}.phase_halves, lfpPhaseCoh.fastPower{iCh}{iCh_ca}.phase_conf_halves] = phaseCohLFP(fastPower2,...
          fastPower1, srData, optLFP);
        lfpPhaseCoh.chans = [chOIDB(iCh) chOIDB_ca(iCh_ca)];
      end
    end
    
    % SAVE DATA
    if exist('lfpPhaseCoh','var')
      dataString = ['dataStruct.seriesData_ca.' animal '_s' num2str(series_c{dbCount}) '__' animal '_s' series_ca{dbCount}{sca} '.lfpphaseCohDataLFP = lfpPhaseCoh;'];
    else
      dataString = ['dataStruct.seriesData_ca.' animal '_s' num2str(series_c{dbCount}) '__' animal '_s' series_ca{dbCount}{sca} '.lfpphaseCohDataLFP = [];'];
    end
    eval(dataString);
    if intermediateSaving
      save(dataFile,'dataStruct','-v7.3'); %#ok<*UNRCH>
    end
  end
end

if ~intermediateSaving
  save(dataFile,'dataStruct','-v7.3');
end
clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca



function [freq, coh, phase, coh_conf, phase_confU, phase_confL, coh_halves_freq, coh_halves, coh_conf_halves, phase_halves,...
  phase_conf_halves] = phaseCohLFP(LFPsignal2, LFPsignal1, srData, optLFP)

[freq, coh, phase, coh_conf, phase_confU, phase_confL, coh_halves_freq, coh_halves, coh_conf_halves, phase_halves,...
  phase_conf_halves] = phaseCohCalc(LFPsignal2-mean(LFPsignal2), LFPsignal1-mean(LFPsignal1), srData, optLFP);
end