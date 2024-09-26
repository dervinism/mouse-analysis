% dataFile = 'M190128_A_MD.mat';
% load(dataFile);

series_caOI = {'20190205085751';'20190312165540'};
fnsData_ca = fieldnames(dataStruct.seriesData_ca);
for iSeries = 1:numel(series_caOI)
  for iCell = 1:numel(fnsData_ca)
    i = strfind(fnsData_ca{iCell},series_caOI{iSeries});
    if ~isempty(i)
      disp(iCell);
    end
  end
end