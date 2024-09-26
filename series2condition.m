function [breakClause, condition] = series2condition(awakeList, anaesthesiaList, seriesName1, seriesName2)
% condition = series2condition(awake, anaesthesia, seriesName1, seriesName2)
%
% Function determines the condition of series (could be either one or two
% series). If two series are supplied, both need to have the same
% condition. Otherwise true breakClause is issued. The output condition
% variable is either 1 for awake and 2 for anaesthesia conditions.

seriesName1Length = numel(seriesName1);
if seriesName1Length > 14
  seriesName1 = seriesName1(1:14);
end

breakClause = false;
if nargin < 4
  for entry = 1:numel(awakeList)
    if strcmpi(awakeList{entry}, seriesName1)
      condition = 1;
      return
    end
  end
  for entry = 1:numel(anaesthesiaList)
    if strcmpi(anaesthesiaList{entry}, seriesName1)
      condition = 2;
      return
    end
  end
else
  seriesName2Length = numel(seriesName2);
  if seriesName2Length > 14
    seriesName2 = seriesName2(1:14);
  end
  for entry = 1:numel(awakeList)
    if strcmpi(awakeList{entry}, seriesName1) && strcmpi(awakeList{entry}, seriesName2)
      condition = 1;
      return
    end
  end
  for entry = 1:numel(anaesthesiaList)
    if strcmpi(anaesthesiaList{entry}, seriesName1) && strcmpi(anaesthesiaList{entry}, seriesName2)
      condition = 2;
      return
    end
  end
end

condition = [];
breakClause = true;