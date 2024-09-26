% This script renames series data subfields.


% LOAD PRE-PROCESSED DATA
dataFile = 'M190128_C_MD.mat';
load(dataFile);


% LOOP THROUGH DB ENTRIES
dataStruct2 = {};
fnsData = fieldnames(dataStruct.seriesData);
for dbCount = 1:numel(fnsData)
  if isfield(dataStruct.seriesData.(fnsData{dbCount}), 'lfpPowerData')
    dataStruct2.seriesData.(fnsData{dbCount}).lfpPowerData = dataStruct.seriesData.(fnsData{dbCount}).lfpPowerData;
    dataStruct.seriesData.(fnsData{dbCount}) = rmfield(dataStruct.seriesData.(fnsData{dbCount}), 'lfpPowerData');
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
save([dataFile(1:end-4) '_lfpFull.mat'],'dataStruct2','-v7.3');