function db = allensdkDB(baseSeries, series)

nDigits = numel(baseSeries);
if nDigits < 14
  baseSeriesLong = [repmat('0',1,14-nDigits) baseSeries];
end

db = [];
for iSeries = 1:numel(series)
  [~, suffixSeries] = determineArea(series{iSeries});
  db(iSeries).series = {str2double(baseSeries); [baseSeriesLong suffixSeries]}; %#ok<*AGROW>
  db(iSeries).basefilename = baseSeries;
end