clearvars -except allData

repository = 'uol';
opt.normType = 'regular';
opt.srDataNew = 10;
opt.rasterType = 'compressed';
opt.outputDir = 'R:\CSN\Shared\Dynamics\Data\mouse_analysis_data\Fig1';
opt.unitIDsTh = [0152 0110]; %[951028663 951024131];
opt.unitIDsCx = [0646 0823 0861 0913]; %[950997930 951019949];
opt.unitIDsHp = [0750 0518]; %[951015782 951017309];
opt.unitIDsDG = [0602 0606]; %[0602 0570];
opt.period = [4000 4600];
animal = 'M191119_A_MD';
recName = '20191202210205';

lists
if strcmp(repository,'all')
  animals = animalsOI;
elseif strcmp(repository,'uol')
  animals = animalsNeuropixelsUOL; %animalsUOLOI;
elseif strcmp(repository,'allensdk')
  animals = animalsAllensdk;
end

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
for iSeriesGroup = 1:numel(seriesOI)
  seriesGroup = seriesOI{iSeriesGroup};
  for iSeries = 1:numel(seriesGroup)
    seriesGroup{iSeries} = [recName seriesGroup{iSeries}];
  end
  unitMUARasterPupilPlotOI(allData.(animal).dataStruct, seriesGroup, opt);
  %close all
end