% Output location
outputFolder = 'D:\infraslow-dynamics\04_data_analysis\001_uol';

% Get recording names
recNames = fieldnames(allData);
nRecs = numel(recNames);

% Create recording data structures
for iRec = 1:nRecs
  dataStruct = allData.(recNames{iRec}).dataStruct;

  % Save data
  filename = fullfile(outputFolder, recNames{iRec}, [recNames{iRec} '.mat']);
  save(filename, 'dataStruct', '-v7.3');
end