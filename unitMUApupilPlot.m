clearvars -except allData

repository = 'uol';
opt.srDataNew = 10;
opt.rasterType = 'compressed';
opt.outputDir = 'R:\CSN\Shared\Dynamics\Data\mouse_analysis_data\rasters\uol\Cx+Hp+pupil\spiking+pupil\M191119_A_MD';
opt.calcRanksOnly = false;
opt.period = [4000 4600];

lists
if strcmp(repository,'all')
  animals = animalsOI;
elseif strcmp(repository,'uol')
  animals = animalsNeuropixelsUOL; %animalsUOLOI;
elseif strcmp(repository,'allensdk')
  animals = animalsAllensdk;
end

%seriesOI = {{'1','56','7','11'},{'14','1921'}};
seriesOI = {{'1','24','56','7','810','11'}};


%% Load data
if ~exist('allData', 'var')
  params
  if strcmp(repository,'all')
    load([dataDir filesep 'allData.mat']);
  elseif strcmp(repository,'uol')
    load([dataDir filesep 'allData_uol.mat']);
  elseif strcmp(repository,'allensdk')
    load([dataDir filesep 'allData_allensdk.mat']);
  end
end


%% Draw raster images
animalList = {};
seriesList = {};
percentMUAsList = [];
percentUnitsList = [];
for animal = 1:numel(animals)
  if isfield(allData, animals{animal})
    if isfield(allData.(animals{animal}).dataStruct, 'eyeData')
      fnsData = fieldnames(allData.(animals{animal}).dataStruct.eyeData);
    else
      continue
    end
    for iRec = 1:numel(fnsData)
      recName = fnsData{iRec};
      for iSeriesGroup = 1:numel(seriesOI)
        seriesGroup = seriesOI{iSeriesGroup};
        for iSeries = 1:numel(seriesGroup)
          seriesGroup{iSeries} = [recName seriesGroup{iSeries}];
        end
        [~, percentages] = unitMUApupilPlotOI(allData.(animals{animal}).dataStruct, seriesGroup, opt);
        if ~isempty(percentages) && opt.calcRanksOnly
          animalList = [animalList; animals{animal}]; %#ok<*AGROW>
          seriesList = [seriesList; recName];
          percentMUAsList = [percentMUAsList; percentages(1)];
          percentUnitsList = [percentUnitsList; percentages(2)];
        end
        close all
      end
    end
  end
end


%% Sort pupil-correlated unit percentages
if opt.calcRanksOnly
  [~, MUAOrder] = sort(percentMUAsList, 'descend');
  [~, unitOrder] = sort(percentUnitsList, 'descend');
  positivityRanks = struct();
  for iSeries = 1:numel(seriesList)
    positivityRanks.(animalList{iSeries}).(seriesList{iSeries}).prctPositiveMUAs = percentMUAsList(iSeries);
    positivityRanks.(animalList{iSeries}).(seriesList{iSeries}).prctPositiveUnits = percentUnitsList(iSeries);
    positivityRanks.(animalList{iSeries}).(seriesList{iSeries}).rankMUAs = find(MUAOrder == iSeries);
    positivityRanks.(animalList{iSeries}).(seriesList{iSeries}).rankUnits = find(unitOrder == iSeries);
  end
  if ~exist(opt.outputDir, 'file')
    mkdir(opt.outputDir);
  end
  save([opt.outputDir filesep 'positivityRanks.mat'], 'positivityRanks', '-v7.3');
end