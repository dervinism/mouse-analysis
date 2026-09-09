function [breakClause, comp, compName, comparisonsList] = series2Comparison(comparisonsList, seriesName1, seriesName2, reverse, augmented)
% [breakClause, comp, compName, comparisonsList] = series2Comparison(comparisonsList, seriesName1, seriesName2, reverse, repository)
% Function tests whether two series are comparable and, if so, identifies
% which comparison it is.
% Input: comparisonsList - a cell array containing a list of comparisons
%                          with area names in the form of
%                          areaName1VsareaName2 or areaName1vsareaName2 or
%                          areaName1vareaName2
%        seriesName1, seriesName2
%        reverse - true if you are interested in reverse comparisons,
%                  default is false.
%        augmented - if true, outputs an augmented comparisons lists
%                    suitable for single probe recordings. Default is
%                    false.
% Output: breakClause - true if two areas (that series correspond to) being
%                       compared are not in the area comparisons list,
%                       false otherwise.
%         comp is the list of comparisons that given series are part of.
%         compName is the list of the corresponding comparison names.
%         comparisonsList - useful to output when using the reverse option.

if nargin < 5
  augmented = false;
end
if nargin < 4
  reverse = false;
end

% Reverse the comparison list if needed
if reverse
  comparisonsList = reverseList(comparisonsList);
end

% Check if comparison is in the list
[~, ~, areaName1] = determineArea(seriesName1);
[~, ~, areaName2] = determineArea(seriesName2);
if strcmp(areaName1(1),'r') || strcmp(areaName1(1),'l')
  areaName1 = areaName1(2:end);
end
if strcmp(areaName1,'VB1') || strcmp(areaName1,'VB2')
  areaName1 = 'VB';
end
if strcmp(areaName1,'Th1') || strcmp(areaName1,'Th2')
  areaName1 = 'Th';
end
if strcmp(areaName2(1),'r') || strcmp(areaName2(1),'l')
  areaName2 = areaName2(2:end);
end
if strcmp(areaName2,'VB1') || strcmp(areaName2,'VB2')
  areaName2 = 'VB';
end
if strcmp(areaName2,'Th1') || strcmp(areaName2,'Th2')
  areaName2 = 'Th';
end
comparison = [areaName1 'Vs' areaName2];
for iComp = 1:numel(comparisonsList)
  if strcmpi(comparison, comparisonsList{iComp})
    comp(1) = iComp;
    compName{1} = comparison;
    breakClause = false;
    break
  elseif iComp == numel(comparisonsList)
    comp = [];
    compName = [];
    breakClause = true;
    return
  end
end

% Check if comparison belongs to a comparison group
comparison = [areaName1(2:end) 'Vs' areaName2];
[comp, compName] = findComp(comparisonsList, comparison, comp, compName);

comparison = [areaName1 'Vs' areaName2(2:end)];
[comp, compName] = findComp(comparisonsList, comparison, comp, compName);

comparison = [areaName1(2:end) 'Vs' areaName2(2:end)];
[comp, compName] = findComp(comparisonsList, comparison, comp, compName);

if augmented
  % 1st quadruple
  comparison = [areaName1(1:end-1) 'Vs' areaName2];
  [comp, compName] = findComp(comparisonsList, comparison, comp, compName);
  
  comparison = [areaName1(2:end-1) 'Vs' areaName2(1:end)];
  [comp, compName] = findComp(comparisonsList, comparison, comp, compName);
  
  comparison = [areaName1(1:end-1) 'Vs' areaName2(2:end)];
  [comp, compName] = findComp(comparisonsList, comparison, comp, compName);
  
  comparison = [areaName1(2:end-1) 'Vs' areaName2(2:end)];
  [comp, compName] = findComp(comparisonsList, comparison, comp, compName);
  
  % 2nd quadruple
  comparison = [areaName1 'Vs' areaName2(1:end-1)];
  [comp, compName] = findComp(comparisonsList, comparison, comp, compName);
  
  comparison = [areaName1(2:end) 'Vs' areaName2(1:end-1)];
  [comp, compName] = findComp(comparisonsList, comparison, comp, compName);
  
  comparison = [areaName1(1:end) 'Vs' areaName2(2:end-1)];
  [comp, compName] = findComp(comparisonsList, comparison, comp, compName);
  
  comparison = [areaName1(2:end) 'Vs' areaName2(2:end-1)];
  [comp, compName] = findComp(comparisonsList, comparison, comp, compName);
  
  % 3rd quadruple
  comparison = [areaName1(1:end-1) 'Vs' areaName2(1:end-1)];
  [comp, compName] = findComp(comparisonsList, comparison, comp, compName);
  
  comparison = [areaName1(2:end-1) 'Vs' areaName2(1:end-1)];
  [comp, compName] = findComp(comparisonsList, comparison, comp, compName);
  
  comparison = [areaName1(1:end-1) 'Vs' areaName2(2:end-1)];
  [comp, compName] = findComp(comparisonsList, comparison, comp, compName);
  
  comparison = [areaName1(2:end-1) 'Vs' areaName2(2:end-1)];
  [comp, compName] = findComp(comparisonsList, comparison, comp, compName);
  
  % Special case: Cx
  % Sub-areas that pool into Cx: S1, RSC (somatosensory/retrosplenial).
  % VIS (visual cortex, relevant for the Allen dataset) is included on its
  % own, NOT alongside its own V1/V2 sub-components: for Allen, VIS is
  % constructed upstream (determineSeries.m) as a combined V1+V2 reference
  % signal built from the SAME source-side units as the individual V1/V2
  % comparisons, not a disjoint unit set - adding V1 and V2 on top of VIS
  % would repeat the same source units up to three times (confirmed via
  % identical total row counts across ThVsV1/ThVsVIS/ThVsV2-style entries
  % in globalUnits_ca.mat).
  if strcmpi(areaName1,'S1') || strcmpi(areaName1,'RSC') || strcmpi(areaName1,'VIS') ||...
      strcmpi(areaName2,'S1') || strcmpi(areaName2,'RSC') || strcmpi(areaName2,'VIS')
    if strcmpi(areaName1,'S1') || strcmpi(areaName1,'RSC') || strcmpi(areaName1,'VIS')
      areaName1 = 'Cx';
    end
    if strcmpi(areaName2,'S1') || strcmpi(areaName2,'RSC') || strcmpi(areaName2,'VIS')
      areaName2 = 'Cx';
    end
    comparison = [areaName1 'Vs' areaName2];
    [comp, compName] = findComp(comparisonsList, comparison, comp, compName);
  end

  % Special case: Hp
  % Sub-areas that pool into Hp: CA (merged CA1+CA3), DG, and the individual
  % CA1/CA3 sub-fields, since determineArea/determineAreaFromSeries never
  % actually returns the merged 'CA' label for a real per-channel UoL
  % recording (only 'CA1' or 'CA3' individually) - without these, units
  % recorded in CA1 or CA3 were silently excluded from the Hp pool entirely.
  % NOTE: this branch is UoL-specific in its actual effect (it is currently
  % never exercised for Allen, since the allensdk call site does not pass
  % augmented=true). If that ever changes, do NOT enable this branch as-is
  % for Allen: Allen's 'Hp'/'CA' labels are themselves pre-combined
  % population-level reference signals built upstream (determineSeries.m),
  % analogous to VIS above, so pooling CA1/CA3/DG alongside Hp/CA would
  % repeat the same source-side units rather than combining disjoint ones -
  % the same class of bug the VIS fix above addresses.
  if strcmpi(areaName1,'CA') || strcmpi(areaName1,'DG') || strcmpi(areaName1,'CA1') || strcmpi(areaName1,'CA3') ||...
      strcmpi(areaName2,'CA') || strcmpi(areaName2,'DG') || strcmpi(areaName2,'CA1') || strcmpi(areaName2,'CA3')
    if strcmpi(areaName1,'CA') || strcmpi(areaName1,'DG') || strcmpi(areaName1,'CA1') || strcmpi(areaName1,'CA3')
      areaName1 = 'Hp';
    end
    if strcmpi(areaName2,'CA') || strcmpi(areaName2,'DG') || strcmpi(areaName2,'CA1') || strcmpi(areaName2,'CA3')
      areaName2 = 'Hp';
    end
    comparison = [areaName1 'Vs' areaName2];
    [comp, compName] = findComp(comparisonsList, comparison, comp, compName);
  end
end



function reverseComparisonsList = reverseList(comparisonsList)

reverseComparisonsList = {};
for iComp = 1:numel(comparisonsList)
  reverseComparisonsList{iComp} = reverseComparison(comparisonsList{iComp}); %#ok<*AGROW>
end


function reverseComparison = reverseComparison(comparison)

strInd = strfind(comparison, 'Vs');
if isempty(strInd)
  strInd = strfind(comparison, 'vs');
end
if isempty(strInd)
  strInd = strfind(comparison, 'v');
  areaName1 = comparison(1:strInd-1);
  areaName2 = comparison(strInd+1:end);
else
  areaName1 = comparison(1:strInd-1);
  areaName2 = comparison(strInd+2:end);
end
reverseComparison = [areaName2 'Vs' areaName1];


function [comp, compName] = findComp(comparisonsList, comparison, comp, compName)

for iComp = 1:numel(comparisonsList)
  if strcmpi(comparison, comparisonsList{iComp})
    comp(numel(comp)+1) = iComp;
    compName{numel(compName)+1} = comparison;
    break
  end
end