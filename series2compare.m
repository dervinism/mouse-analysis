% Determine series to compare

areaIDs = {'1'; '2'; '3'; '24'; '5'; '6'; '56'; '7';...
  '8'; '28'; '10'; '810'; '210'; '11'; '511'; '611';...
  '13'; '14'; '15'; '17'; '18'; '1518'; '19'; '20'; '21'; '2021'; '1921'; '31'; '32'; '38'};

series_c = {};
series_ca = {};
fnsData = fieldnames(dataStruct.seriesData);
for dbCount = 1:numel(fnsData)
  dbStruct = dataStruct.seriesData.(fnsData{dbCount});
  entryName = dbStruct.db(dbCount).entryName;
  
  strSep = strfind(entryName,'s');
  seriesName = entryName(strSep+1:end);
  seriesBase = entryName(strSep+1:strSep+14);
  areaID = entryName(strSep+15:end);
  
  areaInd = startsWith(areaIDs, areaID) & endsWith(areaIDs, areaID);
  if sum(areaInd)
    series_c{numel(series_c)+1} = seriesName; %#ok<*SAGROW>
  else
    continue
  end
  
  series_ca_db = {};
  areaIDs2 = areaIDs(~areaInd);
  for dbCount2 = 1:numel(fnsData)
    dbStruct = dataStruct.seriesData.(fnsData{dbCount2});
    entryName = dbStruct.db(dbCount2).entryName;

    strSep = strfind(entryName,'s');
    seriesName = entryName(strSep+1:end);
    seriesBase2 = entryName(strSep+1:strSep+14);
    areaID2 = entryName(strSep+15:end);
    
    if strcmpi(seriesBase, seriesBase2)
      areaInd2 = startsWith(areaIDs2, areaID2) & endsWith(areaIDs2, areaID2);
      if sum(areaInd2)
        series_ca_db{numel(series_ca_db)+1} = seriesName;
      end
    end
    if dbCount2 == numel(fnsData)
      series_ca{numel(series_ca)+1} = series_ca_db;
    end
  end
end