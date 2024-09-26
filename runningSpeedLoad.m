% A script for loading motion data and detecting frames

% LOAD PRE-PROCESSED DATA
if ~exist('dataStruct','var')
  load(dataFile);
end

fnsData = fieldnames(dataStruct.seriesData);
for dbCount = 1:numel(fnsData)
  % Load the contents of dbStruct
  [~, repository, ~, entryName, baseFilename] = get_dbStruct(dataStruct, dbCount);
  
  if strcmp(repository, 'allensdk')
    
    if ~isfield(dataStruct, 'runningSpeedData')
      % Load the file with running speed data
      runningSpeedData = [baseFilename '.mat'];
      load(runningSpeedData)
      
      % Extract relevant data
      measurementTimes = runningSpeedData(:,1)'; %#ok<*NASGU>
      runningSpeed = runningSpeedData(:,2)'; % m/s
    else
      [seriesName, animal] = seriesFromEntry(fnsData{dbCount});
      fnsData2 = fieldnames(dataStruct.runningSpeedData);
      ind = ismember([animal '_s' seriesName(1:min([14 numel(seriesName)]))], fnsData2);
      measurementTimes = dataStruct.runningSpeedData.(fnsData2{ind}).measurementTimes;
      runningSpeed = dataStruct.runningSpeedData.(fnsData2{ind}).runningSpeed;
    end
    
    % Find periods without running
    sdRunningSpeed = std(runningSpeed);
    meanRunningSpeed = mean(runningSpeed);
    runningThr = 7; %meanRunningSpeed + norminv(1-0.05/2)*sdRunningSpeed;
    runningLocs = runningSpeed > runningThr;
    [~, startTimesInds, startTimes, ~, endTimesInds, endTimes] = findRunningStartEndTimes(measurementTimes, runningLocs);
    for iTime = 2:numel(startTimes)
      if startTimes(iTime) - endTimes(iTime-1) < 1
        runningLocs(endTimesInds(iTime-1):startTimesInds(iTime)) = 1;
      end
    end
    [~, startTimesInds, startTimes, ~, endTimesInds, endTimes] = findRunningStartEndTimes(measurementTimes, runningLocs);
    for iTime = 1:numel(endTimes)
      if endTimes(iTime) - startTimes(iTime) < max([2*min([measurementTimes(2:end) - measurementTimes(1:end-1)]) 0.02])
        runningLocs(startTimesInds(iTime):endTimesInds(iTime)) = 0;
      end
    end
    runningTimesInit = measurementTimes(logical(runningLocs));
    runningLocs = false(size(measurementTimes));
    for iTime = 2:numel(runningTimesInit)
      if runningTimesInit(iTime) - runningTimesInit(iTime-1) < 10
        runningLocs(measurementTimes >= runningTimesInit(iTime-1) & measurementTimes <= runningTimesInit(iTime)) = true;
      end
    end
    runningLocs(1:end-2) = logical(runningLocs(1:end-2) + runningLocs(3:end));
    runningLocs(3:end) = logical(runningLocs(3:end) + runningLocs(1:end-2));
    [~, ~, startQuietTimes, ~, ~, endQuietTimes] = findRunningStartEndTimes(measurementTimes, ~runningLocs); %#ok<*ASGLU>
    [~, maxQuietPeriodInd] = max(endQuietTimes - startQuietTimes);
    maxQuietPeriod = [startQuietTimes(maxQuietPeriodInd) endQuietTimes(maxQuietPeriodInd)];
    
    % Plot data
    figure; plot(measurementTimes, runningSpeed); hold on
    plot(measurementTimes(runningLocs), meanRunningSpeed*ones(1,sum(runningLocs)), 'r.', 'MarkerSize',5); hold off
  else
    runningSpeed = []; measurementTimes = []; runningLocs = []; startQuietTimes = []; endQuietTimes = []; maxQuietPeriod = [];
  end
  
  % Update dataStruct
  [seriesName, animal] = seriesFromEntry(entryName);
  entryName = [animal '_s' seriesName(1:min([14 numel(seriesName)]))];
  dataString = ['dataStruct.runningSpeedData.' entryName '.runningSpeed = runningSpeed;'];
  eval(dataString);
  dataString = ['dataStruct.runningSpeedData.' entryName '.measurementTimes = measurementTimes;'];
  eval(dataString);
  dataString = ['dataStruct.runningSpeedData.' entryName '.runningLocs = runningLocs;'];
  eval(dataString);
  dataString = ['dataStruct.runningSpeedData.' entryName '.startQuietTimes = startQuietTimes;'];
  eval(dataString);
  dataString = ['dataStruct.runningSpeedData.' entryName '.endQuietTimes = endQuietTimes;'];
  eval(dataString);
  dataString = ['dataStruct.runningSpeedData.' entryName '.maxQuietPeriod = maxQuietPeriod;'];
  eval(dataString);
  clear runningSpeed measurementTimes runningLocs startQuietTimes endQuietTimes maxQuietPeriod
end

% Save data
save(dataFile,'dataStruct','-v7.3');
clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca dataStruct



%% Local functions
function [startLocs, startTimesInds, startTimes, endLocs, endTimesInds, endTimes] = findRunningStartEndTimes(measurementTimes, runningLocs)

runningTimesInds = find(logical(runningLocs));
startLocs = [runningLocs(1) runningLocs(2:end) - runningLocs(1:end-1)];
startLocs(startLocs < 1) = 0;
startTimesInds = find(startLocs);
startTimes = measurementTimes(logical(startLocs));
endLocs = [runningLocs(1:end-1) - runningLocs(2:end) runningLocs(end)];
endLocs(endLocs < 1) = 0;
endTimesInds = find(endLocs);
endTimes = measurementTimes(logical(endLocs));
assert(sum(startLocs) == sum(endLocs));
end