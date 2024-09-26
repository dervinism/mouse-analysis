% Run this script to perform phase and coherence comparisons between series.


% LOAD PRE-PROCESSED DATA
if ~exist('dataStruct','var')
  load(dataFile);
end


% INITIALISE PARAMETERS
crossArea = true;
series = [7 13]+5;   % Cx: [9 17] [25 33] [41 49] OR %[10 19] [28 37] [46 55] +8
                    % CxVsTh: [7 13] [19 25] [31 37] OR [55 59] [63 67] [71 75]
                    % CxVsHp: [7 13] [19 25] [31 37]+1 OR [55 59] [63 67] [71 75]+1
                    % Th: [9 17] [25 33] [41 49] +4 OR %[10 19] [28 37] [46 55] +4
                    % ThVsCx: [7 13] [19 25] [31 37]+2 OR [55 59] [63 67] [71 75]+2
                    % ThVsHp: [7 13] [19 25] [31 37]+3
                    % Hp: [9 17] [25 33] [41 49] +7 OR %[10 19] [28 37] [46 55] +7
                    % HpVsCx: [7 13] [19 25] [31 37]+4 OR [55 59] [63 67] [71 75]+4
                    % HpVsTh: [7 13] [19 25] [31 37]+5
seriesStr = {'awake'; 'anaesthesia'};


% LOOP THROUGH DB ENTRIES
phaseCohSeries = {};
if crossArea
  fnsData = fieldnames(dataStruct.seriesData_ca);
else
  fnsData = fieldnames(dataStruct.seriesData); %#ok<*UNRCH>
end
for dbCount = series
  if crossArea
    dbStruct = dataStruct.seriesData_ca.(fnsData{dbCount});
  else
    dbStruct = dataStruct.seriesData.(fnsData{dbCount});
  end
  shankIDs = fieldnames(dbStruct.shankData);
  
% LOAD THE CONTENTS OF THE DB STRUCTURE VARIABLE
if crossArea
  fnsData2 = fieldnames(dataStruct.seriesData);
  FOI = dataStruct.seriesData.(fnsData2{1}).FOI;
else
  FOI = dbStruct.FOI;
end
  
% LOOP THROUGH SHANKS
  for sh = 1:numel(shankIDs)
    shankStruct = dbStruct.shankData.(shankIDs{sh});
    fprintf('%s shank %d -------------------\n', fnsData{dbCount}, sh);
    
% LOAD THE CONTENTS OF THE SHANK STRUCTURE VARIABLE
    units = shankStruct.units;
    uCentres = shankStruct.uCentres;
    if isempty(units)
      continue
    end
    phaseCoh = shankStruct.phaseCoh;
    
% LOOP THROUGH UNITS
    nUnits = numel(units);
    for u = 1:nUnits
      fprintf('Processing %s shank %d unit %d\n', fnsData{dbCount}, sh, u);
      unitData = phaseCoh{u};
      
      [state.phase, state.coh, state.coh_conf, state.FOI, state.fInds] = phaseCohFOI(FOI, unitData.freq, unitData.phase,...
        unitData.coh, unitData.coh_conf, unitData.rateadjust_kappa);
      state.unit = unitData.unit;
      state.uCentre = uCentres(u);
      state.series = fnsData{dbCount};
      
      currentState = seriesStr{find(dbCount==series)}; %#ok<FNDSB>
      if strcmp(currentState,'awake')
        phaseCohSeries{unitData.unit}.awake = state; %#ok<*SAGROW>
      elseif strcmp(currentState,'anaesthesia')
        phaseCohSeries{unitData.unit}.anaesthesia = state;
      end
    end
  end
end


% SAVE DATA
if crossArea
  animal = dataStruct.seriesData.(fnsData2{1}).db(1).animal;
else
  animal = dbStruct.db(series(1)).animal;
end
series1 = fnsData{series(1)}; %dbStruct.db(series(1)).series;
series2 = fnsData{series(2)}; %dbStruct.db(series(2)).series;
if ~crossArea
  i = find(series1 == 's');
  series1 = series1(i+1:end);
  series2 = series2(i+1:end);
end
%dataFileSeries = [dataFile(1:end-4) '_s' num2str(series1) '_s' num2str(series2) '.mat'];
dataFileSeries = [dataFile(1:end-4) '_s' series1 '_s' series2 '.mat'];
save(dataFileSeries,'phaseCohSeries','animal','series1','series2','FOI','-v7.3');

clearvars -except dataFile dataFileSeries dbEntries dbEntries_c dbEntries_ca dataStruct