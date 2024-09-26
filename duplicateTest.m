function breakClause = duplicateTest(compName, seriesData_ca, animal, recording)
% breakClause = duplicateTest(compName, seriesData_ca, animal, recording)
%
% Function determines if a given area comparison is duplicate and, if true,
% issues a true breakClause.
% Input: compName - comparison name.
%        seriesData_ca - series data for cross-area comparisons.
%        animal
%        recording - recording ID.
% Output: breakClause.

duplicateClause = duplicateComparison(compName);

breakClause = false;
if duplicateClause
  if strcmpi(compName, 'S1VsCA1')
    dataString = [animal '_s' recording '1' '__' animal '_s' recording '11'];
  elseif strcmpi(compName, 'S1VsVB1')
    dataString = [animal '_s' recording '1' '__' animal '_s' recording '8'];
  elseif strcmpi(compName, 'S1VsTh1')
    dataString = [animal '_s' recording '1' '__' animal '_s' recording '810'];
  else
    breakClause = true;
    return
  end
  if isfield(seriesData_ca, dataString) && isfield(seriesData_ca.(dataString), 'shankData')
    breakClause = true;
  end
end