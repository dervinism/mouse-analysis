%load('D:\infraslow-dynamics\04_data_analysis\001_uol\noRun\allData_uol.mat');
clearvars -except allData

repository = 'uol';
opt.normType = 'regular';
opt.srDataNew = 10;
opt.rasterType = 'compressed';
outputDir = 'D:\infraslow-dynamics\04_data_analysis\006_analysis_results\noRun\rasters';
if ~exist(outputDir, 'dir')
  mkdir(outputDir)
end
if strcmp(repository, 'uol')
  if strcmp(opt.normType, 'normalised')
    opt.outputDir = fullfile(outputDir, 'uol\Th+Cx+Hp_normalised');
  else
    opt.outputDir = fullfile(outputDir, 'uol\Th+Cx+Hp');
  end
elseif strcmp(repository, 'allensdk') %#ok<*UNRCH>
  if strcmp(opt.normType, 'normalised')
    opt.outputDir = fullfile(outputDir, 'allen\Th+Cx+Hp_normalised');
  else
    opt.outputDir = fullfile(outputDir, 'allen\Th+Cx+Hp');
  end
end

lists
if strcmp(repository,'all')
  animals = animalsOI;
elseif strcmp(repository,'uol')
  animals = animalsNeuropixelsUOL; %animalsUOLOI;
elseif strcmp(repository,'allensdk')
  animals = animalsAllensdk;
end

% seriesOI = {{'1'},{'7'},{'24','810'},{'56','11'},{'1','56','7','11'},{'1','24','56','7','810','11'},...
%   {'14'},{'15'},{'18'},{'1518'},{'19'},{'20'},{'21'},{'1921'},{'14','1921'}};
% seriesOI = {{'1','56','7','11'}};
seriesOI = {{'1','24','56','7','810','11'},{'14','1518','1921'}};


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
for animal = [9 11] %1:numel(animals) %[9 11] %4:numel(animals)
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
      rasterPlotOI(allData.(animals{animal}).dataStruct, seriesGroup, opt);
      close all
    end
  end
end