% Wrapping script for the detectFrames function as used with total movement data

motionDataFiles
fnsData = fieldnames(dataStruct.seriesData);
probe = dataStruct.seriesData.(fnsData{1}).conf.probe;
for iData = 1:numel(motionData)
  dbStruct = dataStruct.seriesData.(fnsData{iData});
  if ~isempty(motionData{iData}) && ~isempty(dbStruct.db(iData).eyeCameraCh)
    clear frameTimes frameInd
    load(motionData{iData});
    if ~exist('frameTimes','var') || ~exist('frameInd','var')
      fprintf(1, 'Detecting frame times in %s\n', motionData{iData});
      dbStruct = dataStruct.seriesData.(fnsData{iData});
      if strcmpi(probe, 'Neuropixels')
        dataFilename = [dbStruct.io.baseFilename '.bin'];
      else
        dataFilename = [dbStruct.io.baseFilename '.dat'];
      end
      opt.updateMovementData = 'smaller';
      opt.movementDataFile = motionData{iData};
      opt.dataCrop = croppingInstructions{iData};
      syncFuncDualNeuropix(dataFilename, [], dbStruct.db(iData).chN, [], dbStruct.conf.params.srRecording, opt)
    end
  end
end