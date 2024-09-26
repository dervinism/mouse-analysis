function [dbStruct, repository, dataDir, entryName, baseFilename, probe, shankIDs,...
  chOI, period, periodLFP, srRecording, srData, optCoh, exclRad, FOI,...
  MUAsAll, spkDB, spkDB_units, phaseCoh, phaseCohHalves, muaMetadata] = get_dbStruct(dataStruct, dbCount, subpop)

if nargin < 3
  fnsData = fieldnames(dataStruct.seriesData);
  dbStruct = dataStruct.seriesData.(fnsData{dbCount});
elseif strcmp(subpop, 'positive')
  fnsData = fieldnames(dataStruct.seriesData_positive);
  dbStruct = dataStruct.seriesData_positive.(fnsData{dbCount});
elseif strcmp(subpop, 'negative')
  fnsData = fieldnames(dataStruct.seriesData_negative);
  dbStruct = dataStruct.seriesData_negative.(fnsData{dbCount});
elseif strcmp(subpop, 'neutral')
  fnsData = fieldnames(dataStruct.seriesData_neutral);
  dbStruct = dataStruct.seriesData_neutral.(fnsData{dbCount});
end

if ~isfield(dbStruct,'db')
  dbStruct = [];
  repository = [];
  dataDir = [];
  entryName = [];
  baseFilename = [];
  probe = [];
  shankIDs = [];
  chOI = [];
  period = [];
  periodLFP = [];
  srRecording = [];
  srData = [];
  optCoh = [];
  exclRad = [];
  FOI = [];
  MUAsAll = [];
  spkDB = [];
  spkDB_units = [];
  phaseCoh = [];
  phaseCohHalves = [];
  muaMetadata = [];
  return
end
repository = dbStruct.db(dbCount).repository;
dataDir = dbStruct.io.dataDir;
entryName = dbStruct.db(dbCount).entryName;
baseFilename = dbStruct.io.baseFilename;
probe = dbStruct.conf.probe;
shankIDs = fieldnames(dbStruct.shankData);
chOI = dbStruct.db(dbCount).chOI;
period = dbStruct.db(dbCount).period;
periodLFP = dbStruct.db(dbCount).periodLFP;
srRecording = dbStruct.conf.samplingParams.srRecording;
srData = dbStruct.conf.samplingParams.srData;

optCoh = dbStruct.conf.optCoh;
exclRad = dbStruct.conf.exclRad;
FOI = dbStruct.conf.FOI;

MUAsAll = dbStruct.popData.MUAsAll;
spkDB = dbStruct.popData.spkDB;
spkDB_units = dbStruct.popData.spkDB_units;
if isfield(dbStruct.popData, 'muaMetadata')
  muaMetadata = dbStruct.popData.muaMetadata;
else
  muaMetadata = [];
end

if isfield(dbStruct.popData, 'phaseCoh')
  phaseCoh = dbStruct.popData.phaseCoh;
else
  phaseCoh = [];
end
if isfield(dbStruct.popData, 'phaseCohHalves')
  phaseCohHalves = dbStruct.popData.phaseCohHalves;
else
  phaseCohHalves = [];
end