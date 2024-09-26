function breakClause = duplicateComparison(comparison)
% breakClause = duplicateComparison(comparison)
%
% Function checks if area comparison is potentially duplicate. A comparison
% is potentially duplicate if area1 on probe1 is being compared to area2 on
% probe2 when area2 potentially exists on probe1. When comparison is
% potentially duplicate a true breakClause is issued. S1 comparisons can
% only be potential duplicates as more information is needed to establish
% it with certainty (whether probe 1 extends beyond S1). Comparisons not
% involving S1 are real duplicates.

unwantedComparisonList = {
  'VB1VsCA3';
  'CA3VsVB1';
  'VB2VsCA1';
  'CA1VsVB2';
  'Th1VsCA3';
  'CA3VsTh1';
  'Th2VsCA1';
  'CA1VsTh2';
  'S1VsCA1';
  'S1VsVB1';
  'S1VsTh1';
  'RSCVsCA3';
  'RSCVsVB2';
  'RSCVsTh2'};

for iComp = 1:numel(unwantedComparisonList)
  if strcmpi(comparison, unwantedComparisonList{iComp})
    breakClause = true;
    return
  else
    breakClause = false;
  end
end