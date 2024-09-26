function [adjPhase, phaseROI, phaseCI95ROI] = adjustPhase(phase, freq, phaseCI95, cutoffFreq)
% [adjPhase, phaseROI, phaseCI95ROI] = adjustPhase(phase, freq, phaseCI95, cutoffFreq)
%
% Function adjusts phase values in the best possible way so far.
% Input: phase, freq, phaseCI95, cutoffFreq (default is 30 Hz).
% Output: adjPhase - adjusted phase values in radians.
%         phaseROI - adjusted phase for the frequency region of interest.
%         phaseCI95ROI - adjusted phase 95% confidence interval for the
%                        frequency region of interest. Can be left empty if
%                        not needed.

if nargin < 4
  cutoffFreq = 30;
end

adjPhase = bestUnwrap(phase);
if ~isempty(adjPhase)
  phaseROI = adjPhase(:,freq <= cutoffFreq);
  if (max(phaseROI(~isnan(phaseROI)))+min(phaseROI(~isnan(phaseROI))))/2 < -pi %mean(options.phaseLim)-pi
    adjPhase = adjPhase + 2*pi;
  elseif (max(phaseROI(~isnan(phaseROI)))+min(phase(~isnan(phaseROI))))/2 > pi %mean(options.phaseLim)+pi
    adjPhase = adjPhase - 2*pi;
  end
else
  return
end
phaseROI = adjPhase(freq <= cutoffFreq);
if median(phaseROI(~isnan(phaseROI))) < -pi %mean(options.phaseLim)-pi
  adjPhase = adjPhase + 2*pi;
elseif median(phaseROI(~isnan(phaseROI))) > pi %mean(options.phaseLim)+pi
  adjPhase = adjPhase - 2*pi;
end
phaseROI = adjPhase(freq <= cutoffFreq);
if isempty(phaseCI95)
  phaseCI95ROI = [];
else
  phaseCI95ROI = phaseCI95(freq <= cutoffFreq);
end