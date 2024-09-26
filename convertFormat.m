function [convertedFilename, dataStruct] = convertFormat(path2file)
% convertedFilename = convertFormat(path2file)
%
% Function converts Martynas's data format to the one that Michael prefers.
% Input: path2file - full file path.
% Output: convertedFilename.
%         dataStruct - data structure after conversion.


load(path2file) %#ok<*LOAD>


%% Allocate new data fields
fnsData = fieldnames(dataStruct.seriesData); %#ok<*NODEF>
dataStructNew = struct();
prevSeries = '';
recCount = 0;
for dbCount = 1:numel(fnsData)
  [seriesName, animal] = seriesFromEntry(fnsData{dbCount});
  [~, ~, area] = determineArea(seriesName);
  currentSeries = seriesName(1:min([numel(seriesName) 14]));
  if ~strcmp(prevSeries, currentSeries)
    if dbCount == 1
      dataStructNew.animalID = animal;
    end
    recCount = recCount + 1;
    dataStructNew.experiments{recCount}.recordingID = currentSeries;
    dataStructNew.experiments{recCount}.areas = {};
    dataStructNew.experiments{recCount}.areaComparisons = {};
    dataStructNew.experiments{recCount}.eyeData = struct();
    dataStructNew.experiments{recCount}.motionData = struct();
    dataStructNew.experiments{recCount}.runningSpeedData = struct();
    areaCount = 0;
    prevSeries = currentSeries;
  end
  areaCount = areaCount + 1;
  dataStructNew.experiments{recCount}.areas{areaCount}.areaAcronym = area;
end

if isfield(dataStruct, 'seriesData_ca')
  fnsData_ca = fieldnames(dataStruct.seriesData_ca);
  for dbCount = 1:numel(fnsData_ca)
    [seriesName1, seriesName2] = seriesNames(fnsData_ca{dbCount});
    assert(strcmp(seriesName1(1:min([numel(seriesName1) 14])), seriesName2(1:min([numel(seriesName2) 14]))));
    [~, ~, area1] = determineArea(seriesName1);
    [~, ~, area2] = determineArea(seriesName2);
    for iExp = 1:numel(dataStructNew.experiments)
      if strcmp(dataStructNew.experiments{iExp}.recordingID, seriesName1(1:min([numel(seriesName1) 14])))
        dataStructNew.experiments{iExp}.areaComparisons{numel(dataStructNew.experiments{iExp}.areaComparisons)+1}.areaAcronym1 = area1;
        dataStructNew.experiments{iExp}.areaComparisons{numel(dataStructNew.experiments{iExp}.areaComparisons)}.areaAcronym2 = area2;
        dataStructNew.experiments{iExp}.areaComparisons{numel(dataStructNew.experiments{iExp}.areaComparisons)}.areaComparisonAcronym = [area1 '_wrt_' area2];
        break
      end
    end
  end
end


%% Allocate data to new fields
topFields = fieldnames(dataStruct);
for iField = 1:numel(topFields)
  if strcmp(topFields{iField}, 'popData') || strcmp(topFields{iField}, 'popData_ca') ||...
      strcmp(topFields{iField}, 'popData_ca_positive') || strcmp(topFields{iField}, 'popData_ca_negative') ||...
      strcmp(topFields{iField}, 'popData_ca_neutral')
    continue
  elseif strcmp(topFields{iField}, 'eyeData') || strcmp(topFields{iField}, 'motionData') || strcmp(topFields{iField}, 'runningSpeedData')
    fnsData = fieldnames(dataStruct.(topFields{iField}));
    for dbCount = 1:numel(fnsData)
      seriesName = seriesFromEntry(fnsData{dbCount});
      for iExp = 1:numel(dataStructNew.experiments)
        if strcmp(dataStructNew.experiments{iExp}.recordingID, seriesName(1:min([numel(seriesName) 14])))
          dataStructNew.experiments{iExp}.(topFields{iField}) = dataStruct.(topFields{iField}).(fnsData{dbCount});
          if isfield(dataStructNew.experiments{iExp}.(topFields{iField}), 'shankData')
            fnsShankData = fieldnames(dataStructNew.experiments{iExp}.(topFields{iField}).shankData);
            for iShank = 1:numel(fnsShankData)
              dataStructNew.experiments{iExp}.(topFields{iField}).shankData.shank{iShank} =...
                dataStructNew.experiments{iExp}.(topFields{iField}).shankData.(['shank' num2str(iShank)]);
              dataStructNew.experiments{iExp}.(topFields{iField}).shankData =...
                rmfield(dataStructNew.experiments{iExp}.(topFields{iField}).shankData, ['shank' num2str(iShank)]);
            end
          end
          break
        end
      end
    end
  else
    fnsData = fieldnames(dataStruct.(topFields{iField}));
    for dbCount = 1:numel(fnsData)
      try
        [seriesName1, seriesName2] = seriesNames(fnsData{dbCount});
        assert(strcmp(seriesName1(1:min([numel(seriesName1) 14])), seriesName2(1:min([numel(seriesName2) 14]))));
        [~, ~, area1] = determineArea(seriesName1);
        [~, ~, area2] = determineArea(seriesName2);
        for iExp = 1:numel(dataStructNew.experiments)
          if strcmp(dataStructNew.experiments{iExp}.recordingID, seriesName1(1:min([numel(seriesName1) 14])))
            for iAreaComp = 1:numel(dataStructNew.experiments{iExp}.areaComparisons)
              if strcmp(dataStructNew.experiments{iExp}.areaComparisons{iAreaComp}.areaComparisonAcronym, [area1 '_wrt_' area2])
                dataStructNew.experiments{iExp}.areaComparisons{iAreaComp}.(topFields{iField}) = dataStruct.(topFields{iField}).(fnsData{dbCount});
                if isfield(dataStructNew.experiments{iExp}.areaComparisons{iAreaComp}.(topFields{iField}), 'shankData')
                  fnsShankData = fieldnames(dataStructNew.experiments{iExp}.areaComparisons{iAreaComp}.(topFields{iField}).shankData);
                  for iShank = 1:numel(fnsShankData)
                    dataStructNew.experiments{iExp}.areaComparisons{iAreaComp}.(topFields{iField}).shankData.shank{iShank} =...
                      dataStructNew.experiments{iExp}.areaComparisons{iAreaComp}.(topFields{iField}).shankData.(['shank' num2str(iShank)]);
                    dataStructNew.experiments{iExp}.areaComparisons{iAreaComp}.(topFields{iField}).shankData =...
                      rmfield(dataStructNew.experiments{iExp}.areaComparisons{iAreaComp}.(topFields{iField}).shankData, ['shank' num2str(iShank)]);
                  end
                end
                break
              end
            end
            break
          end
        end
      catch
        try
          seriesName = seriesFromEntry(fnsData{dbCount});
          [~, ~, area] = determineArea(seriesName);
          for iExp = 1:numel(dataStructNew.experiments)
            if strcmp(dataStructNew.experiments{iExp}.recordingID, seriesName(1:min([numel(seriesName) 14])))
              for iArea = 1:numel(dataStructNew.experiments{iExp}.areas)
                if strcmp(dataStructNew.experiments{iExp}.areas{iArea}.areaAcronym, area)
                  dataStructNew.experiments{iExp}.areas{iArea}.(topFields{iField}) = dataStruct.(topFields{iField}).(fnsData{dbCount});
                  if isfield(dataStructNew.experiments{iExp}.areas{iArea}.(topFields{iField}), 'shankData')
                    fnsShankData = fieldnames(dataStructNew.experiments{iExp}.areas{iArea}.(topFields{iField}).shankData);
                    for iShank = 1:numel(fnsShankData)
                      dataStructNew.experiments{iExp}.areas{iArea}.(topFields{iField}).shankData.shank{iShank} =...
                        dataStructNew.experiments{iExp}.areas{iArea}.(topFields{iField}).shankData.(['shank' num2str(iShank)]);
                      dataStructNew.experiments{iExp}.areas{iArea}.(topFields{iField}).shankData =...
                        rmfield(dataStructNew.experiments{iExp}.areas{iArea}.(topFields{iField}).shankData, ['shank' num2str(iShank)]);
                    end
                  end
                  break
                end
              end
              break
            end
          end
        catch me
          disp(me);
          error('Potential error in the original data structure');
        end
      end
    end
  end
end


%% Save data
path = fileparts(path2file);
convertedFilename = [path filesep animal '_converted.mat'];
dataStruct = dataStructNew;
save(convertedFilename, 'dataStruct', '-v7.3');