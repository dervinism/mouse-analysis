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
  if strcmpi(areaName1,'S1') || strcmpi(areaName1,'RSC') ||...
      strcmpi(areaName2,'S1') || strcmpi(areaName2,'RSC')
    if strcmpi(areaName1,'S1') || strcmpi(areaName1,'RSC')
      areaName1 = 'Cx';
    end
    if strcmpi(areaName2,'S1') || strcmpi(areaName2,'RSC')
      areaName2 = 'Cx';
    end
    comparison = [areaName1 'Vs' areaName2];
    [comp, compName] = findComp(comparisonsList, comparison, comp, compName);
  end
  
  % Special case: Hp
  if strcmpi(areaName1,'CA') || strcmpi(areaName1,'DG') ||...
      strcmpi(areaName2,'CA') || strcmpi(areaName2,'DG')
    if strcmpi(areaName1,'CA') || strcmpi(areaName1,'DG')
      areaName1 = 'Hp';
    end
    if strcmpi(areaName2,'CA') || strcmpi(areaName2,'DG')
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