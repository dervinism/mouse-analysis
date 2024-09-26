function phasePeak = circPeak(phase, edges, range)
% phasePeak = circPeak(phase, edges, range)
%
% Function finds a peak phase or the most common phase in the phase sample.
% Inputs: phase - phase sample vector.
%         edges - histcounts bin edges for producing a phase histogram.
%         range - a vector defining the ends of the phase range of
%                 interest. For example, [-pi pi].
% Output:phasePeak - the most common phase value.

% Produce phase histogram
phase = phase(~isnan(phase));
phaseHist = histcounts(phase, edges);
locs = edges(1:end-1) + (edges(2)-edges(1))/2;

% Identify range start index
rangeStart = locs - range(1);
rangeStart(rangeStart >= 0) = 1;
rangeStart(rangeStart < 0) = 0;
iRangeStart = find(rangeStart,1,'first');
if isempty(iRangeStart)
  iRangeStart = 1;
end

% Identify range end index
rangeEnd = locs - range(2);
rangeEnd(rangeEnd >= 0) = 1;
rangeEnd(rangeEnd < 0) = 0;
iRangeEnd = find(rangeEnd,1,'first');
if isempty(iRangeEnd)
  iRangeEnd = numel(locs);
end

% Find the peak
[~, iPeak] = max(phaseHist(iRangeStart:iRangeEnd));
phasePeak = locs(iRangeStart + iPeak - 1);