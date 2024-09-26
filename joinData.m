fnsData = fieldnames(a.dataStruct.seriesData_ca);
for iField = 1:numel(fnsData)
  dataStruct.seriesData_ca.(fnsData{iField}) = a.dataStruct.seriesData_ca.(fnsData{iField})
end

numel(fieldnames(dataStruct.seriesData_ca))

seriesCount = 0;
for iSeries = 1:numel(series_ca)
  seriesCount = seriesCount + numel(series_ca{iSeries});
end
seriesCount