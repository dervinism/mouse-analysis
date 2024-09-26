% Wrapping script for the detectFrames function as used with pupil area data

eyeDataFiles
fnsData = fieldnames(dataStruct.seriesData);
probe = dataStruct.seriesData.(fnsData{1}).conf.probe;
for iData = 1:numel(eyeData)
  dbStruct = dataStruct.seriesData.(fnsData{iData});
  if ~isempty(eyeData{iData}) && ~isempty(dbStruct.db(iData).eyeCameraCh)
    load(eyeData{iData});
    if ~isfield(results,'frameTimes') || ~isfield(results,'frameInd')
      fprintf(1, 'Detecting frame times in %s\n', eyeData{iData});
      dbStruct = dataStruct.seriesData.(fnsData{iData});
      if strcmpi(probe, 'Neuropixels')
        dataFilename = [dbStruct.io.baseFilename '.bin'];
      else
        dataFilename = [dbStruct.io.baseFilename '.dat'];
      end
      opt.updatePupilData = 'smaller';
      opt.pupilDataFile = eyeData{iData};
      opt.dataCrop = croppingInstructions{iData};
      syncFuncDualNeuropix(dataFilename, [], dbStruct.db(iData).chN, [], dbStruct.conf.params.srRecording, opt)
    end
  end
end