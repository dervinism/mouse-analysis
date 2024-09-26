% Copy NSG dependencies

load(dataFile)

fnsData = fieldnames(dataStruct.seriesData);
if ~isempty(dbEntries) && dbEntries(1) == inf
  dbEntries = 1:numel(fnsData);
end
for dbCount = dbEntries
  dbStruct = dataStruct.seriesData.(fnsData{dbCount});
  entryName = dbStruct.db(dbCount).entryName;
  strSep = strfind(entryName, 's');
  seriesName = entryName(strSep+1:end);
  baseFilename = dbStruct.io.baseFilename;
  sourceDir = fileparts(baseFilename);
  
  if ~exist(seriesName,'dir')
    mkdir(seriesName)
  end
  copyfile([baseFilename '.qua.1.mat'], seriesName);
  copyfile([sourceDir filesep 'forPRB*'], seriesName);
  if exist([sourceDir filesep 'waveforms.mat'], 'file')
    copyfile([sourceDir filesep 'waveforms.mat'], seriesName);
  else
    disp(['Non-existent waveform file: ' sourceDir filesep 'waveforms.mat']);
  end
end