function [commonPeriod, commonInds] = combinePeriods(period1, period2, sr)
% Given two time period arrays and the data sampling rate, determine the
% overlapping period.

if ~iscell(period1)
  period1 = {period1};
end
if ~iscell(period2)
  period2 = {period2};
end

if numel(period1) == 1 && numel(period2) == 1
  if max([period1{1}(1) period2{1}(1)]) >= min([period1{1}(2) period2{1}(2)])
    commonPeriod = [];
  else
    commonPeriod = [max([period1{1}(1) period2{1}(1)]) min([period1{1}(2) period2{1}(2)])];
  end
else
  period1vec = [];
  for iP = 1:numel(period1)
    period1vec = [period1vec round(period1{iP}(1)*sr):round(period1{iP}(2)*sr)];
  end
  period2vec = [];
  for iP = 1:numel(period2)
    period2vec = [period2vec round(period2{iP}(1)*sr):round(period2{iP}(2)*sr)];
  end
  commonPeriodVec = intersect(period1vec, period2vec);
  commonPeriodVecDiff = diff(commonPeriodVec);
  commonPeriodVecDiff = [commonPeriodVecDiff(1) commonPeriodVecDiff];
  [~, pksInds] = findpeaks(commonPeriodVecDiff);
  if isempty(pksInds)
    commonPeriod = [commonPeriodVec(1)/sr commonPeriodVec(end)/sr];
  else
    commonPeriod = {};
    for iPks = 1:numel(pksInds)
      if iPks == 1
        commonPeriod{1} = [commonPeriodVec(1)/sr commonPeriodVec(pksInds(iPks)-1)/sr];
      else
        commonPeriod{iPks} = [commonPeriodVec(pksInds(iPks-1))/sr commonPeriodVec(pksInds(iPks)-1)/sr]; %#ok<*AGROW>
      end
    end
    commonPeriod{iPks+1} = [commonPeriodVec(pksInds(iPks))/sr commonPeriodVec(end)/sr];
  end
end

if nargout > 1
  if iscell(commonPeriod)
    commonInds = {};
    for iCell = 1:numel(commonPeriod)
      commonInds{iCell} = round(commonPeriod{iCell}.*sr);
      commonInds{iCell}(1) = max([1 commonInds{iCell}(1)]);
    end
  else
    commonInds = round(commonPeriod.*sr);
    commonInds(1) = max([1 commonInds(1)]);
  end
end