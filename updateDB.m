dbCount = 6;
dbStruct = dataStruct.seriesData.(fnsData{dbCount});
dbStruct.db = db;
dbStruct.io.baseFilename = [dbStruct.io.dataDir filesep db(dbCount).basefilename];
dbEntry = db(dbCount).entryName;
seriesData = dataStruct.seriesData;
seriesData = setfield(seriesData,dbEntry,dbStruct); %#ok<SFLD>
dataStruct.seriesData = seriesData;
save(dataFile,'dataStruct','-v7.3');