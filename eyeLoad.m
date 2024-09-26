% A script for loading eye data and detecting frames

% LOAD PRE-PROCESSED DATA
if ~exist('dataStruct','var')
  load(dataFile);
end

eyeDataFiles % eyeData, period, and croppingInstructions variables are defined inside
fnsData = fieldnames(dataStruct.seriesData);
for dbCount = 1:numel(fnsData)
  
  % Load the contents of dbStruct
  [dbStruct, repository, ~, entryName, baseFilename, probe, ~,...
    ~, ~, ~, srRecording] = get_dbStruct(dataStruct, dbCount);
  
  % Load eye data
  if strcmp(repository, 'uol')
    if ~isempty(eyeData{dbCount}) && ~isempty(dbStruct.db(dbCount).eyeCameraCh)
      if strcmp(eyeData{dbCount}(end-2:end), 'mat') % Eliptic fit data (most is no longer used except for a few instances)
        load(eyeData{dbCount});
        if ~isfield(results,'frameTimes') || ~isfield(results,'frameInd')
          fprintf(1, 'Detecting frame times in %s\n', eyeData{dbCount});
          if strcmpi(probe, 'Neuropixels')
            dataFilename = [baseFilename '.bin'];
          else
            dataFilename = [baseFilename '.dat'];
          end
          opt.updatePupilData = 'smaller'; % options parameters are explained in syncFuncDualNeuropix
          opt.pupilDataFile = eyeData{dbCount};
          opt.dataCrop = croppingInstructions{dbCount};
          [~, ~, results] = syncFuncDualNeuropix(dataFilename, [], dbStruct.db(dbCount).chN, [], srRecording, opt); % load pupil data
        end
        pupilArea = torow(results.area); %#ok<*NASGU>
        if isfield(results,'blink')
          blink = torow(results.blink);
        end
        frameTimes = torow(results.frameTimes);
        periodDB = period{dbCount};
      elseif strcmp(eyeData{dbCount}(end-2:end), 'csv') % deeplabcut data (most of the eye-tracking data now comes in this format)
        [seriesName, animal] = seriesFromEntry(entryName);
        if exist('prevEntryName', 'var')
          prevSeriesName = seriesFromEntry(prevEntryName);
        else
          prevSeriesName = '00000000000000';
        end
        if (numel(seriesName) >= 14 && strcmpi(seriesName(1:14), prevSeriesName(1:14))) ||...
            (numel(seriesName) < 14 && strcmpi(seriesName, prevSeriesName)) % Do not reload data if the previous series has the same basename meaning they share pupil data
          pupilArea = prevPupilArea;
          frameTimes = prevFrameTimes;
        else
          if isfield(dataStruct,'eyeData') &&...
              ((numel(seriesName) >= 14 && isfield(dataStruct.eyeData, [animal '_s' seriesName(1:14)]) && isfield(dataStruct.eyeData.([animal '_s' seriesName(1:14)]),'frameTimes')) ||...
              (numel(seriesName) < 14 && isfield(dataStruct.eyeData, [animal '_s' seriesName]) && isfield(dataStruct.eyeData.([animal '_s' seriesName]),'frameTimes'))) % if eye data was loaded already some time ago
            pupilArea = torow(calcPupilArea(eyeData{dbCount}));
            try
              frameTimes = dataStruct.eyeData.([animal '_s' seriesName(1:14)]).frameTimes;
            catch
              frameTimes = dataStruct.eyeData.([animal '_s' seriesName]).frameTimes;
            end
            if numel(pupilArea) ~= numel(frameTimes) % check if the existing and newly loaded data match in size. If not, detect frames again.
              fprintf(1, 'Detecting frame times in %s\n', eyeData{dbCount});
              if strcmpi(probe, 'Neuropixels')
                dataFilename = [baseFilename '.bin'];
              else
                dataFilename = [baseFilename '.dat'];
              end
              opt.updatePupilData = 'smaller';
              opt.pupilDataFile = eyeData{dbCount};
              opt.dataCrop = croppingInstructions{dbCount};
              [~, ~, results] = syncFuncDualNeuropix(dataFilename, [], dbStruct.db(dbCount).chN, [], srRecording, opt);
              pupilArea = torow(results.area);
              frameTimes = torow(results.frameTimes);
            end
          else % If no previously loaded data exist
            fprintf(1, 'Detecting frame times in %s\n', eyeData{dbCount});
            if strcmpi(probe, 'Neuropixels')
              dataFilename = [baseFilename '.bin'];
            else
              dataFilename = [baseFilename '.dat'];
            end
            opt.updatePupilData = 'smaller';
            opt.pupilDataFile = eyeData{dbCount};
            opt.dataCrop = croppingInstructions{dbCount};
            [~, ~, results] = syncFuncDualNeuropix(dataFilename, [], dbStruct.db(dbCount).chN, [], srRecording, opt);
            pupilArea = torow(results.area);
            frameTimes = torow(results.frameTimes);
          end
        end
        periodDB = period{dbCount};
      end
    else
      pupilArea = []; frameTimes = []; periodDB = [];
    end
  else
    eyeData = [baseFilename '.mat'];
    load(eyeData)
    if exist('gazeData','var')
      [pupilArea, frameTimes] = eyeLoad_allensdk(eyeData); %#ok<*ASGLU>
      periodDB = period;
    else
      pupilArea = []; frameTimes = []; periodDB = [];
    end
  end
  
  % Crop period
  if ~isempty(periodDB)
    minTime = min(frameTimes);
    maxTime = max(frameTimes);
    if iscell(periodDB)
      periodDBnew = [];
      for iCell = 1:numel(periodDB)
        if periodDB{iCell}(1) >= maxTime
          break
        elseif iCell == 1
          periodDBnew = [max([periodDB{iCell}(1) minTime]) min([periodDB{iCell}(2) maxTime])];
        elseif iCell == 2
          periodDBnew = {periodDBnew; [max([periodDB{iCell}(1) minTime]) min([periodDB{iCell}(2) maxTime])]}; %#ok<*AGROW>
        else
          periodDBnew{iCell} = [max([periodDB{iCell}(1) minTime]) min([periodDB{iCell}(2) maxTime])]; %#ok<*SAGROW>
        end
      end
      periodDB = periodDBnew;
    else
      periodDB = [max([periodDB(1) minTime]) min([periodDB(2) maxTime])];
    end
  end
  
  % Update dataStruct
  [seriesName, animal] = seriesFromEntry(entryName);
  entryName = [animal '_s' seriesName(1:14)];
  dataString = ['dataStruct.eyeData.' entryName '.pupilArea = pupilArea;'];
  eval(dataString);
  if exist('blink','var')
    dataString = ['dataStruct.eyeData.' entryName '.blink = blink;'];
    eval(dataString);
  elseif isfield(dataStruct.eyeData.(entryName), 'blink')
    dataStruct.eyeData.(entryName) = rmfield(dataStruct.eyeData.(entryName), 'blink');
  end
  dataString = ['dataStruct.eyeData.' entryName '.frameTimes = frameTimes;'];
  eval(dataString);
  if ~isempty(pupilArea)
    dataString = ['dataStruct.eyeData.' entryName '.period = periodDB;'];
  else
    dataString = ['dataStruct.eyeData.' entryName '.period = [];'];
  end
  eval(dataString);
  
  prevEntryName = entryName;
  prevPupilArea = pupilArea;
  prevFrameTimes = frameTimes;
  clear pupilArea blink frameTimes periodDB
end

% Save data
save(dataFile,'dataStruct','-v7.3');
clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca dataStruct