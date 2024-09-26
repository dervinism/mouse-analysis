SOI = c;
fieldNames = fieldnames(SOI.seriesData_ca);
for i = 1:numel(fieldNames)
  dataStruct.seriesData_ca.(fieldNames{i}) = SOI.seriesData_ca.(fieldNames{i});
end