fnsData = fieldnames(a.dataStruct.seriesData_ca);
for iField = 1:numel(fnsData)
  dataStruct.seriesData_ca.(fnsData{iField}) = a.dataStruct.seriesData_ca.(fnsData{iField});
end

fnsData = fieldnames(a.dataStruct.seriesData_ca_positive);
for iField = 1:numel(fnsData)
  dataStruct.seriesData_ca_positive.(fnsData{iField}) = a.dataStruct.seriesData_ca_positive.(fnsData{iField});
end

fnsData = fieldnames(a.dataStruct.seriesData_ca_negative);
for iField = 1:numel(fnsData)
  dataStruct.seriesData_ca_negative.(fnsData{iField}) = a.dataStruct.seriesData_ca_negative.(fnsData{iField});
end

fnsData = fieldnames(a.dataStruct.seriesData_ca_neutral);
for iField = 1:numel(fnsData)
  dataStruct.seriesData_ca_neutral.(fnsData{iField}) = a.dataStruct.seriesData_ca_neutral.(fnsData{iField});
end