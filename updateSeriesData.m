% This script renames series data subfields.


% LOAD PRE-PROCESSED DATA
dataFile = 'M190523_B_MD.mat';
load(dataFile);


% LOOP THROUGH DB ENTRIES
fnsData = fieldnames(dataStruct.seriesData);
for dbCount = 1:numel(fnsData)
  if isfield(dataStruct.seriesData.(fnsData{dbCount}), 'lfpphaseCohData')
    dataStruct.seriesData.(fnsData{dbCount}).lfpphaseCohDataMotion = dataStruct.seriesData.(fnsData{dbCount}).lfpphaseCohData;
    dataStruct.seriesData.(fnsData{dbCount}) = rmfield(dataStruct.seriesData.(fnsData{dbCount}), 'lfpphaseCohData');
  end
%   if isfield(dataStruct.seriesData.(fnsData{dbCount}), 'lfpPowerDataMotion')
%     dataStruct.seriesData.(fnsData{dbCount}) = rmfield(dataStruct.seriesData.(fnsData{dbCount}), 'lfpPowerDataMotion');
%   end
end % loop over db entries

save(dataFile,'dataStruct','-v7.3');