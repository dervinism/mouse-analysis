% This script renames series data subfields.


% LOAD PRE-PROCESSED DATA
dataFile = 'M190128_A_MD.mat';
load(dataFile);


% LOOP THROUGH DB ENTRIES
dataStruct2 = {};
fnsData = fieldnames(dataStruct.seriesData);
for dbCount = 1:numel(fnsData)
  if isfield(dataStruct.seriesData.(fnsData{dbCount}), 'lfpPowerData')
    if isfield(dataStruct.seriesData.(fnsData{dbCount}).lfpPowerData, 'areaInterpFilt')
      dataStruct2.seriesData.(fnsData{dbCount}).lfpPowerData.areaInterpFilt = dataStruct.seriesData.(fnsData{dbCount}).lfpPowerData.areaInterpFilt;
      dataStruct2.seriesData.(fnsData{dbCount}).lfpPowerData.areaInterpTimes = dataStruct.seriesData.(fnsData{dbCount}).lfpPowerData.areaInterpTimes;
    end
    if isfield(dataStruct.seriesData.(fnsData{dbCount}).lfpPowerData, 'motionInterpFilt')
      dataStruct2.seriesData.(fnsData{dbCount}).lfpPowerData.motionInterpFilt = dataStruct.seriesData.(fnsData{dbCount}).lfpPowerData.motionInterpFilt;
      dataStruct2.seriesData.(fnsData{dbCount}).lfpPowerData.motionInterpTimes = dataStruct.seriesData.(fnsData{dbCount}).lfpPowerData.motionInterpTimes;
    end
    if isfield(dataStruct.seriesData.(fnsData{dbCount}).lfpPowerData, 'areaInterpFilt')
      dataStruct.seriesData.(fnsData{dbCount}).lfpPowerData = rmfield(dataStruct.seriesData.(fnsData{dbCount}).lfpPowerData, 'areaInterpFilt');
      dataStruct.seriesData.(fnsData{dbCount}).lfpPowerData = rmfield(dataStruct.seriesData.(fnsData{dbCount}).lfpPowerData, 'areaInterpTimes');
    end
    if isfield(dataStruct.seriesData.(fnsData{dbCount}).lfpPowerData, 'motionInterpFilt')
      dataStruct.seriesData.(fnsData{dbCount}).lfpPowerData = rmfield(dataStruct.seriesData.(fnsData{dbCount}).lfpPowerData, 'motionInterpFilt');
      dataStruct.seriesData.(fnsData{dbCount}).lfpPowerData = rmfield(dataStruct.seriesData.(fnsData{dbCount}).lfpPowerData, 'motionInterpTimes');
    end
  else
    dataStruct2.seriesData.(fnsData{dbCount}).lfpPowerData = [];
  end
  if isfield(dataStruct.seriesData.(fnsData{dbCount}), 'lfpphaseCohData')
    dataStruct2.seriesData.(fnsData{dbCount}).lfpphaseCohData = dataStruct.seriesData.(fnsData{dbCount}).lfpphaseCohData;
    dataStruct.seriesData.(fnsData{dbCount}) = rmfield(dataStruct.seriesData.(fnsData{dbCount}), 'lfpphaseCohData');
  else
    dataStruct2.seriesData.(fnsData{dbCount}).lfpphaseCohData = [];
  end
  if isfield(dataStruct.seriesData.(fnsData{dbCount}), 'lfpphaseCohDataMotion')
    dataStruct2.seriesData.(fnsData{dbCount}).lfpphaseCohDataMotion = dataStruct.seriesData.(fnsData{dbCount}).lfpphaseCohDataMotion;
    dataStruct.seriesData.(fnsData{dbCount}) = rmfield(dataStruct.seriesData.(fnsData{dbCount}), 'lfpphaseCohDataMotion');
  else
    dataStruct2.seriesData.(fnsData{dbCount}).lfpphaseCohDataMotion = [];
  end
end % loop over db entries

save(dataFile,'dataStruct','-v7.3');
save([dataFile(1:end-4) '_lfp.mat'],'dataStruct2','-v7.3');