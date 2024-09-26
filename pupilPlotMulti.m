clearvars -except allData
params
lists

repository = 'all';

if strcmp(repository,'all')
  animals = animalsOI;
elseif strcmp(repository,'uol')
  animals = animalsUOLOI;
elseif strcmp(repository,'allensdk')
  animals = animalsAllensdk;
end
mainFolder = [dataDir filesep pupilDir];


%% Load data
if ~exist('allData', 'var')
  params
  if strcmp(repository,'all') && exist([dataDir filesep 'allData.mat'], 'file')
    load([dataDir filesep 'allData.mat']);
  elseif strcmp(repository,'uol') && exist([dataDir filesep 'allData_uol.mat'], 'file')
    load([dataDir filesep 'allData_uol.mat']);
  elseif strcmp(repository,'allensdk') && exist([dataDir filesep 'allData_allensdk.mat'], 'file')
    load([dataDir filesep 'allData_allensdk.mat']);
  end
end


%% Plot pupil area
for animal = 1:numel(animals) % Loop through animals
  disp(['processing animal ' animals{animal}]);
  if exist('allData', 'var')
    try
      dataStruct = allData.(animals{animal}).dataStruct;
    catch
      load([dataDir filesep animals{animal} filesep animals{animal} '.mat'])
    end
  else
    load([dataDir filesep animals{animal} filesep animals{animal} '.mat'])
  end
  fnsData = fieldnames(dataStruct.seriesData);
  
  for dbCount = 1:numel(fnsData) % Loop through database entries
    seriesName = seriesFromEntry(fnsData{dbCount});
    if dbCount == 1 || ~strcmp(prevSeriesName, seriesName(1:min([numel(seriesName) 14])))
      prevSeriesName = seriesName(1:min([numel(seriesName) 14]));
      plotData = true;
    else
      plotData = false;
    end
    
    % Determine if series pupil data exist
    eyeDataEntry = [animals{animal} '_s' seriesName(1:min([numel(seriesName) 14]))];
    if ~isfield(dataStruct, 'eyeData') || ~isfield(dataStruct.eyeData, eyeDataEntry)
      continue
    end
    
    % Test for exceptions
%     if exceptionTest(except, seriesName)
%       prevSeriesName = 'exception occured';
%       continue
%     end
    
    % Plot the pupil area and save the figure
    if plotData
      plotEyeOrMotion(dataStruct, fnsData{dbCount}, [], 'pupil', [], mainFolder, true);
    end
  end
end