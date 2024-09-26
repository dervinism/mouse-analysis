function [pvalMWW, pvalWEst, pvalWObs] = dispPhaseHistoStats(fileName, area, frequency, display)
% dispPhaseHistoStats(fileName, area, frequency)
%
% Function displays non-parametric histogram mean comparison statistics
% (Mardia-Watson-Wheeler Uniform-Scores Test and Watson Test).
% Input: fileName - the name of the file containing histogram data.
%        area - area comparison of interest.
%        frequency - frequency of interest
%        display - true if you want to print out the stats in the command
%          window, default is false.
% Output: pvalMWW - p-value of the Mardia-Watson-Wheeler Uniform-Scores
%                   Test.
%         pvalWEst - estimated p-value of Watson Test.
%         pvalWObs - observed average p-value of Watson Test after running
%                    multiple permutation tests.

if nargin < 4
  display = false;
end

if numel(area) > 6 && strcmp(area(end-6:end), 'VsPupil')
  area = area(1:end-7);
elseif numel(area) > 7 && strcmp(area(end-7:end), 'VsMotion')
  area = area(1:end-8);
end

load(fileName) %#ok<*LOAD>

% Determine the area index
for iArea = 1:numel(areas) %#ok<*USENS>
  if strcmpi(area, areas{iArea})
    break
  end
  if iArea == numel(areas)
    error('The area comparison you supplied does not exist')
  end
end

% Determine the frequency index
for iF = 1:numel(FOI)
  if frequency == FOI(iF)
    break
  end
  if iF == numel(FOI)
    error('The frequency value you supplied does not exist')
  end
end

% Display stats
if isempty(fPEst)
  pvalMWW = NaN;
else
  pvalMWW = fPEst{iArea}(iF);
end
if isempty(pEst)
  pvalWEst = NaN;
else
  pvalWEst = pEst{iArea}(iF);
end
if isempty(pObs)
  pvalWObs = NaN;
else
  pvalWObs = pObs{iArea}(iF);
end
if display
  disp(['p-value of Mardia-Watson-Wheeler Uniform-Scores Test: ' num2str(pvalMWW)]);
  disp(['Estimated p-value of Watson Test: ' num2str(pvalWEst)]);
  disp(['Observed p-value of Watson Test: ' num2str(pvalWObs)]);
end