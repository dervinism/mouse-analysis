% Synchronisation function wrap script

% Initialise probe setup variables and arousal measure files
if iscell(db(1).series)
  dataDir = [topDir filesep db(1).animal filesep num2str(db(1).series{1})];
else
  dataDir = [topDir filesep db(1).animal filesep num2str(db(1).series)];
end
if strcmp(repository, 'uol')
  fileList = dir([dataDir filesep 'forPRB*.mat']);
  if isempty(fileList)
    error('The data processing folder is missing the forPRB*.mat file containing probe configuration.');
  elseif numel(fileList) > 1
    error('There is more than one forPRB*.mat file in the data processing folder. Please remove conflicting files.');
  else
    probeSetup = load([dataDir filesep fileList.name]);
    probe = fileList.name(8:end-4);
    if strcmpi(probe,'mecP3opt3')
      probe = 'imecP3opt3';
    end
  end
  
  eyeDataFiles
  motionDataFiles
end

% Initialise synchronisation function using the existing one or else create a blank one.
fnsData = {db.entryName};
syncFunc = {};
for dbCount = 1:length(db)
  if exist('dataStruct','var') && isfield(dataStruct, 'seriesData') && dbCount <= numel(dataStruct.seriesData)
    dbStruct = dataStruct.seriesData.(fnsData{dbCount});
    if isfield(dbStruct.conf,'syncFunc')
      syncFunc{dbCount} = dbStruct.conf.syncFunc; %#ok<*SAGROW>
    else
      syncFunc{dbCount}.a = 0; syncFunc{dbCount}.b = 1;
    end
  else
    syncFunc{dbCount}.a = 0; syncFunc{dbCount}.b = 1;
  end
end

% Calculate the synchronisation function (only applies to uol dual neuropixels recordings)
if overwrite && syncFuncCalc
  if ~isempty(dbEntries) && dbEntries(1) == inf
    dbEntries = 1:numel(fnsData);
  end
  for dbCount = dbEntries
    if strcmp(repository, 'uol') && strcmpi(probe,'Neuropixels')
      [series, animal] = seriesFromEntry(fnsData{dbCount});
      if endsWith([animal '_s' series(1:14) '1'], fnsData) && endsWith([animal '_s' series(1:14) '2'], fnsData)
        if iscell(db(dbCount).series)
          dataDir = [topDir filesep db(dbCount).animal filesep num2str(db(dbCount).series{1})];
        else
          dataDir = [topDir filesep db(dbCount).animal filesep num2str(db(dbCount).series)];
        end
        if ~isempty(db(dbCount).eyeCameraCh) && numel(seriesFromEntry(fnsData{dbCount})) == 15 && strcmpi(fnsData{dbCount}(end), '1')...
            && dbCount ~= numel(fnsData) && numel(seriesFromEntry(fnsData{dbCount+1})) == 15 && strcmpi(fnsData{dbCount+1}(end), '2')
          baseFilename = [dataDir filesep db(dbCount).basefilename];
          dataFilename1 = [baseFilename '.bin'];
          if iscell(db(dbCount+1).series)
            dataDir = [topDir filesep db(dbCount+1).animal filesep num2str(db(dbCount+1).series{1})];
          else
            dataDir = [topDir filesep db(dbCount+1).animal filesep num2str(db(dbCount+1).series)];
          end
          baseFilename = [dataDir filesep db(dbCount+1).basefilename];
          dataFilename2 = [baseFilename '.bin'];
          
          %         if ~isempty(eyeData) && ~isempty(eyeData{dbCount})
          %           opt2.updatePupilData = 'smaller';
          %           opt2.pupilDataFile = eyeData{dbCount};
          %         else
          opt2.updatePupilData = 'no';
          %         end
          %         if ~isempty(motionData) && ~isempty(motionData{dbCount})
          %           opt2.updateMovementData = 'smaller';
          %           opt2.movementDataFile = motionData{dbCount};
          %         else
          opt2.updateMovementData = 'no';
          %         end
          opt2.dataCrop = croppingInstructions{dbCount};
          [syncFunc_1to2, syncFunc_2to1] = syncFuncDualNeuropix(dataFilename1, dataFilename2, db(dbCount).chN, db(dbCount+1).chN,...
            srRecording, opt2);
          if syncFunc_1to2.b < syncFunc_2to1.b
            syncFunc{dbCount} = syncFunc_1to2;
          else
            syncFunc{dbCount}.a = 0; syncFunc{dbCount}.b = 1;
          end
        elseif ~isempty(db(dbCount).eyeCameraCh)...
            && (strcmpi(fnsData{dbCount}(end), '1') || strcmpi(fnsData{dbCount}(end), '8') ||...
            strcmpi(fnsData{dbCount}(end), '9') || strcmpi(fnsData{dbCount}(end-1:end), '10') ||...
            strcmpi(fnsData{dbCount}(end-1:end), '11') || strcmpi(fnsData{dbCount}(end-2:end), '810'))
          if exist('syncFunc_1to2','var')
            if syncFunc_1to2.b < syncFunc_2to1.b
              syncFunc{dbCount} = syncFunc_1to2;
            else
              syncFunc{dbCount}.a = 0; syncFunc{dbCount}.b = 1;
            end
          else
            syncFunc{dbCount}.a = 0; syncFunc{dbCount}.b = 1;
          end
        elseif ~isempty(db(dbCount).eyeCameraCh)...
            && (strcmpi(fnsData{dbCount}(end), '2') || strcmpi(fnsData{dbCount}(end), '3') ||...
            strcmpi(fnsData{dbCount}(end), '4') || strcmpi(fnsData{dbCount}(end), '5') ||...
            strcmpi(fnsData{dbCount}(end), '6') || strcmpi(fnsData{dbCount}(end), '7') ||...
            strcmpi(fnsData{dbCount}(end-1:end), '24') || strcmpi(fnsData{dbCount}(end-1:end), '56'))
          if exist('syncFunc_1to2','var')
            if syncFunc_1to2.b < syncFunc_2to1.b
              syncFunc{dbCount}.a = 0; syncFunc{dbCount}.b = 1;
            else
              syncFunc{dbCount} = syncFunc_2to1;
            end
          else
            syncFunc{dbCount}.a = 0; syncFunc{dbCount}.b = 1;
          end
        end
      else
        syncFunc{dbCount}.a = 0; syncFunc{dbCount}.b = 1;
      end
    else
      syncFunc{dbCount}.a = 0; syncFunc{dbCount}.b = 1;
    end
  end
end