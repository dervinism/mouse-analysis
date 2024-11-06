% Downsample conversion script

clear all %#ok<*CLALL>

% Load data
params
D = dir('M*_old.mat');
dataFilename = D.name;
oldDataContainer = load(dataFilename);
dataStruct = oldDataContainer.dataStruct;
vectorWise = true;

% Loop over structure fields and update them
fields = fieldnames(dataStruct);
for iField = 1:numel(fields)
  if strcmpi(fields{iField}, 'seriesData')
    dataType = 'seriesData';
    dataFields = fieldnames(dataStruct.(dataType));
    for jField = 1:numel(dataFields)
      if ~isempty(oldDataContainer.dataStruct.(dataType).(dataFields{jField}))

        % Downsample unit data
        oldSR = oldDataContainer.dataStruct.(dataType).(dataFields{jField}).conf.samplingParams.srData;
        if oldSR > srData
          if ~isempty(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).shankData.shank1) && ...
            ~isempty(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).shankData.shank1.MUAs)
            if ~isempty(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).shankData.shank1.MUAs)
              dataStruct.(dataType).(dataFields{jField}).shankData.shank1.MUAs = ...
                resampleSpikeCounts(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).shankData.shank1.MUAs, ...
                stepsize=1/oldSR, newStepsize=1/srData);
            end
            if ~isempty(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk)
              try
                if vectorWise
                  error('Skipping to processing by vectors');
                end
                dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk = ...
                  sparse(resampleSpikeCounts(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk, ...
                  stepsize=1/oldSR, newStepsize=1/srData)); %#ok<*UNRCH>
              catch
                nUnits = size(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk,1);
                for iUnit = 1:nUnits
                  spkRow = resampleSpikeCounts(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk(iUnit,:), ...
                    stepsize=1/oldSR, newStepsize=1/srData);
                  if iUnit == 1
                    dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk = zeros(nUnits, numel(spkRow));
                  end
                  dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk(iUnit,:) = spkRow;
                  if iUnit == nUnits
                    dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk = ...
                      sparse(dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk);
                  end
                end
              end
            end
          end

          % Downsample MUA data
          if ~isempty(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).popData.MUAsAll)
            if ~isempty(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).popData.MUAsAll)
              dataStruct.(dataType).(dataFields{jField}).popData.MUAsAll = ...
                resampleSpikeCounts(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).popData.MUAsAll, ...
                stepsize=1/oldSR, newStepsize=1/srData);
            end
            if ~isempty(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).popData.spkDB)
              try
                if vectorWise
                  error('Skipping to processing by vectors');
                end
                dataStruct.(dataType).(dataFields{jField}).popData.spkDB = ...
                  sparse(resampleSpikeCounts(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).popData.spkDB, ...
                  stepsize=1/oldSR, newStepsize=1/srData));
              catch
                nUnits = size(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).popData.spkDB,1);
                for iUnit = 1:nUnits
                  spkRow = resampleSpikeCounts(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).popData.spkDB(iUnit,:), ...
                    stepsize=1/oldSR, newStepsize=1/srData);
                  if iUnit == 1
                    dataStruct.(dataType).(dataFields{jField}).popData.spkDB = zeros(nUnits, numel(spkRow));
                  end
                  dataStruct.(dataType).(dataFields{jField}).popData.spkDB(iUnit,:) = spkRow;
                  if iUnit == nUnits
                    dataStruct.(dataType).(dataFields{jField}).popData.spkDB = ...
                      sparse(dataStruct.(dataType).(dataFields{jField}).popData.spkDB);
                  end
                end
              end
            end
          end
        end

        % Update parameters
        dataStruct.(dataType).(dataFields{jField}).conf.optCoh.maxFreq = maxFreq;
        dataStruct.(dataType).(dataFields{jField}).conf.optCoh.winfactor = winfactor;
        dataStruct.(dataType).(dataFields{jField}).conf.optCoh.freqfactor = freqfactor;
        dataStruct.(dataType).(dataFields{jField}).conf.optCoh.tapers = tapers;
        dataStruct.(dataType).(dataFields{jField}).conf.exclRad = exclRad;
        dataStruct.(dataType).(dataFields{jField}).conf.FOI = FOI;
        dataStruct.(dataType).(dataFields{jField}).conf.samplingParams.srData = srData;
        for iEntry = 1:numel(dataStruct.(dataType).(dataFields{jField}).db)
          dataStruct.(dataType).(dataFields{jField}).db(iEntry).period = [0 99999];
          dataStruct.(dataType).(dataFields{jField}).db(iEntry).periodLFP = [0 99999];
        end
        dataStruct.(dataType).(dataFields{jField}).dbSeries.period = [0 99999];
        dataStruct.(dataType).(dataFields{jField}).dbSeries.periodLFP = [0 99999];
      end
    end
  elseif strcmpi(fields{iField}, 'seriesData_positive')
    dataType = 'seriesData_positive';
    dataFields = fieldnames(dataStruct.(dataType));
    for jField = 1:numel(dataFields)
      if ~isempty(oldDataContainer.dataStruct.(dataType).(dataFields{jField}))

        % Downsample unit data
        oldSR = oldDataContainer.dataStruct.(dataType).(dataFields{jField}).conf.samplingParams.srData;
        if oldSR > srData
          if ~isempty(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).shankData.shank1) && ...
            ~isempty(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).shankData.shank1.MUAs)
            if ~isempty(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).shankData.shank1.MUAs)
              dataStruct.(dataType).(dataFields{jField}).shankData.shank1.MUAs = ...
                resampleSpikeCounts(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).shankData.shank1.MUAs, ...
                stepsize=1/oldSR, newStepsize=1/srData);
            end
            if ~isempty(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk)
              try
                if vectorWise
                  error('Skipping to processing by vectors');
                end
                dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk = ...
                  sparse(resampleSpikeCounts(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk, ...
                  stepsize=1/oldSR, newStepsize=1/srData));
              catch
                nUnits = size(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk,1);
                for iUnit = 1:nUnits
                  spkRow = resampleSpikeCounts(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk(iUnit,:), ...
                    stepsize=1/oldSR, newStepsize=1/srData);
                  if iUnit == 1
                    dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk = zeros(nUnits, numel(spkRow));
                  end
                  dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk(iUnit,:) = spkRow;
                  if iUnit == nUnits
                    dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk = ...
                      sparse(dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk);
                  end
                end
              end
            end
          end

          % Downsample MUA data
          if ~isempty(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).popData.MUAsAll)
            if ~isempty(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).popData.MUAsAll)
              dataStruct.(dataType).(dataFields{jField}).popData.MUAsAll = ...
                resampleSpikeCounts(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).popData.MUAsAll, ...
                stepsize=1/oldSR, newStepsize=1/srData);
            end
            if ~isempty(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).popData.spkDB)
              try
                if vectorWise
                  error('Skipping to processing by vectors');
                end
                dataStruct.(dataType).(dataFields{jField}).popData.spkDB = ...
                  sparse(resampleSpikeCounts(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).popData.spkDB, ...
                  stepsize=1/oldSR, newStepsize=1/srData));
              catch
                nUnits = size(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).popData.spkDB,1);
                for iUnit = 1:nUnits
                  spkRow = resampleSpikeCounts(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).popData.spkDB(iUnit,:), ...
                    stepsize=1/oldSR, newStepsize=1/srData);
                  if iUnit == 1
                    dataStruct.(dataType).(dataFields{jField}).popData.spkDB = zeros(nUnits, numel(spkRow));
                  end
                  dataStruct.(dataType).(dataFields{jField}).popData.spkDB(iUnit,:) = spkRow;
                  if iUnit == nUnits
                    dataStruct.(dataType).(dataFields{jField}).popData.spkDB = ...
                      sparse(dataStruct.(dataType).(dataFields{jField}).popData.spkDB);
                  end
                end
              end
            end
          end
        end

        % Update parameters
        dataStruct.(dataType).(dataFields{jField}).conf.optCoh.maxFreq = maxFreq_ca;
        dataStruct.(dataType).(dataFields{jField}).conf.optCoh.winfactor = winfactor;
        dataStruct.(dataType).(dataFields{jField}).conf.optCoh.freqfactor = freqfactor;
        dataStruct.(dataType).(dataFields{jField}).conf.optCoh.tapers = tapers;
        dataStruct.(dataType).(dataFields{jField}).conf.exclRad = exclRad;
        dataStruct.(dataType).(dataFields{jField}).conf.FOI = FOI;
        dataStruct.(dataType).(dataFields{jField}).conf.samplingParams.srData = srData;
        for iEntry = 1:numel(dataStruct.(dataType).(dataFields{jField}).db)
          dataStruct.(dataType).(dataFields{jField}).db(iEntry).period = [0 99999];
          dataStruct.(dataType).(dataFields{jField}).db(iEntry).periodLFP = [0 99999];
        end
        dataStruct.(dataType).(dataFields{jField}).dbSeries.period = [0 99999];
        dataStruct.(dataType).(dataFields{jField}).dbSeries.periodLFP = [0 99999];
      end
    end
  elseif strcmpi(fields{iField}, 'seriesData_negative')
    dataType = 'seriesData_negative';
    dataFields = fieldnames(dataStruct.(dataType));
    for jField = 1:numel(dataFields)
      if ~isempty(oldDataContainer.dataStruct.(dataType).(dataFields{jField}))

        % Downsample unit data
        oldSR = oldDataContainer.dataStruct.(dataType).(dataFields{jField}).conf.samplingParams.srData;
        if oldSR > srData
          if ~isempty(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).shankData.shank1) && ...
            ~isempty(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).shankData.shank1.MUAs)
            if ~isempty(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).shankData.shank1.MUAs)
              dataStruct.(dataType).(dataFields{jField}).shankData.shank1.MUAs = ...
                resampleSpikeCounts(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).shankData.shank1.MUAs, ...
                stepsize=1/oldSR, newStepsize=1/srData);
            end
            if ~isempty(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk)
              try
                if vectorWise
                  error('Skipping to processing by vectors');
                end
                dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk = ...
                  sparse(resampleSpikeCounts(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk, ...
                  stepsize=1/oldSR, newStepsize=1/srData));
              catch
                nUnits = size(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk,1);
                for iUnit = 1:nUnits
                  spkRow = resampleSpikeCounts(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk(iUnit,:), ...
                    stepsize=1/oldSR, newStepsize=1/srData);
                  if iUnit == 1
                    dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk = zeros(nUnits, numel(spkRow));
                  end
                  dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk(iUnit,:) = spkRow;
                  if iUnit == nUnits
                    dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk = ...
                      sparse(dataStruct.(dataType).(dataFields{jField}).shankData.shank1.spk);
                  end
                end
              end
            end
          end

          % Downsample MUA data
          if ~isempty(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).popData.MUAsAll)
            if ~isempty(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).popData.MUAsAll)
              dataStruct.(dataType).(dataFields{jField}).popData.MUAsAll = ...
                resampleSpikeCounts(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).popData.MUAsAll, ...
                stepsize=1/oldSR, newStepsize=1/srData);
            end
            if ~isempty(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).popData.spkDB)
              try
                if vectorWise
                  error('Skipping to processing by vectors');
                end
                dataStruct.(dataType).(dataFields{jField}).popData.spkDB = ...
                  sparse(resampleSpikeCounts(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).popData.spkDB, ...
                  stepsize=1/oldSR, newStepsize=1/srData));
              catch
                nUnits = size(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).popData.spkDB,1);
                for iUnit = 1:nUnits
                  spkRow = resampleSpikeCounts(oldDataContainer.dataStruct.(dataType).(dataFields{jField}).popData.spkDB(iUnit,:), ...
                    stepsize=1/oldSR, newStepsize=1/srData);
                  if iUnit == 1
                    dataStruct.(dataType).(dataFields{jField}).popData.spkDB = zeros(nUnits, numel(spkRow));
                  end
                  dataStruct.(dataType).(dataFields{jField}).popData.spkDB(iUnit,:) = spkRow;
                  if iUnit == nUnits
                    dataStruct.(dataType).(dataFields{jField}).popData.spkDB = ...
                      sparse(dataStruct.(dataType).(dataFields{jField}).popData.spkDB);
                  end
                end
              end
            end
          end
        end

        % Update parameters
        dataStruct.(dataType).(dataFields{jField}).conf.optCoh.maxFreq = maxFreq_ca;
        dataStruct.(dataType).(dataFields{jField}).conf.optCoh.winfactor = winfactor;
        dataStruct.(dataType).(dataFields{jField}).conf.optCoh.freqfactor = freqfactor;
        dataStruct.(dataType).(dataFields{jField}).conf.optCoh.tapers = tapers;
        dataStruct.(dataType).(dataFields{jField}).conf.exclRad = exclRad;
        dataStruct.(dataType).(dataFields{jField}).conf.FOI = FOI;
        dataStruct.(dataType).(dataFields{jField}).conf.samplingParams.srData = srData;
        for iEntry = 1:numel(dataStruct.(dataType).(dataFields{jField}).db)
          dataStruct.(dataType).(dataFields{jField}).db(iEntry).period = [0 99999];
          dataStruct.(dataType).(dataFields{jField}).db(iEntry).periodLFP = [0 99999];
        end
        dataStruct.(dataType).(dataFields{jField}).dbSeries.period = [0 99999];
        dataStruct.(dataType).(dataFields{jField}).dbSeries.periodLFP = [0 99999];
      end
    end
  end
end

% Save the new data
newDataFilename = [dataFilename(1:end-8) '.mat'];
save(newDataFilename, 'dataStruct', '-v7.3');