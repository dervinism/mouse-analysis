% A script for loading motion data and detecting frames

% LOAD PRE-PROCESSED DATA
if ~exist('dataStruct','var')
  load(dataFile);
end

motionDataFiles % motionData, period, and croppingInstructions variables are defined inside
fnsData = fieldnames(dataStruct.seriesData);
for dbCount = 1:numel(fnsData)
  
  % Load the contents of dbStruct
  [dbStruct, repository, ~, entryName, baseFilename, probe, ~,...
    ~, ~, ~, srRecording] = get_dbStruct(dataStruct, dbCount);
  
  % Load motion data
  if strcmp(repository, 'uol')
    if ~isempty(motionData{dbCount}) && ~isempty(dbStruct.db(dbCount).eyeCameraCh)
      load(motionData{dbCount});
      if ~exist('frameTimes','var') || ~exist('frameInd','var')
        fprintf(1, 'Detecting frame times in %s\n', motionData{dbCount});
        if strcmpi(probe, 'Neuropixels')
          dataFilename = [baseFilename '.bin'];
        else
          dataFilename = [baseFilename '.dat'];
        end
        opt.updateMotionData = 'smaller';
        opt.motionDataFile = motionData{dbCount};
        opt.dataCrop = croppingInstructions{dbCount};
        [~, ~, results] = syncFuncDualNeuropix(dataFilename, [], dbStruct.db(dbCount).chN, [], srRecording, opt);
        s = results.s;
        sa = results.sa;
        frameTimes = results.frameTimes;
      end
      s = torow(s);
      sa = torow(sa);
      frameTimes = torow(frameTimes);
      periodDB = period{dbCount}; %#ok<*NASGU>
    else
      s = []; sa = []; frameTimes = []; periodDB = [];
    end
  else
    s = []; sa = []; frameTimes = []; periodDB = [];
  end
  
  % Update dataStruct
  [seriesName, animal] = seriesFromEntry(entryName);
  entryName = [animal '_s' seriesName(1:14)];
  dataString = ['dataStruct.motionData.' entryName '.s = s;'];
  eval(dataString);
  dataString = ['dataStruct.motionData.' entryName '.sa = sa;'];
  eval(dataString);
  dataString = ['dataStruct.motionData.' entryName '.frameTimes = frameTimes;'];
  eval(dataString);
  if ~isempty(s)
    dataString = ['dataStruct.motionData.' entryName '.period = periodDB;'];
  else
    dataString = ['dataStruct.motionData.' entryName '.period = [];'];
  end
  eval(dataString);
  clear s sa frameTimes periodDB
end

% Save data
save(dataFile,'dataStruct','-v7.3');
clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca dataStruct