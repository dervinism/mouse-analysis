% A part of the old AnPSD script adapted to perform coherence and phase
% analyses for unit and population spiking rate data.


% LOAD PRE-PROCESSED DATA
dataFile = 'M190128_A_MD.mat';
load(dataFile);


% LOOP THROUGH DB ENTRIES
fnsData = fieldnames(dataStruct.seriesData);
for dbCount = 1:numel(fnsData)
  dbStruct = dataStruct.seriesData.(fnsData{dbCount});
  shankIDs = fieldnames(dbStruct.shankData);
  
% LOAD THE CONTENTS OF THE DB STRUCTURE VARIABLE
  dataDir = dbStruct.io.dataDir;
  baseFilename = dbStruct.io.baseFilename;
  opt = dbStruct.conf.opt;
  srRecording = dbStruct.conf.params.srRecording;
  srData = dbStruct.conf.params.srData;
  probe = dbStruct.conf.probe;
  MUAsGrps = dbStruct.popData.MUAsGrps;
  if ~iscell(MUAsGrps)
    MUAsGrps = {MUAsGrps};
  end
  spkDB = dbStruct.popData.spkDB;
  spkDB_units = dbStruct.popData.spkDB_units;
  entryName = dbStruct.db(dbCount).entryName;
  chOI = dbStruct.db(dbCount).chOI;
  FOI = dbStruct.FOI;
  exclRad = dbStruct.exclRad;
  
% CORRECT ENTRIES
  strSep = strfind(baseFilename, 'No');
  correctedFilename = [baseFilename(1:strSep-1) baseFilename(strSep+2:end)];
  dbStruct.io.baseFilename = correctedFilename;
  dbStruct.db(dbCount).baseFilename = correctedFilename;
  
  dataStruct.seriesData.(fnsData{dbCount}) = dbStruct;
  save(dataFile,'dataStruct','-v7.3');
  
  fprintf('Finished processing db entry %i\n',dbCount);
end % loop over db entries
clearvars
