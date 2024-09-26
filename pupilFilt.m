function [areaFilt, interpTimes, areaInd, area] = pupilFilt(eyeDataDB, sr, ephysSignal, filtArtefact, period, srRecording)

plotFigs = false;

area = abs(eyeDataDB.pupilArea);
frameTimes = eyeDataDB.frameTimes;

% Descriptive measures
if plotFigs
  %figure; plot(frameTimes, area); hold on %#ok<*UNRCH>
end
%meanArea = mean(area,'omitnan');
medianArea = median(area,'omitnan');
%stdArea = std(area,'omitnan');
%stdAreaTop = medianArea+stdArea;
%stdAreaBottom = medianArea-stdArea;

% Interpolate blinks and noise and resample the data
if isfield(eyeDataDB,'blink')
  area(eyeDataDB.blink) = NaN;
end
area(area > 20*medianArea) = NaN;
area(area < 0.05*medianArea) = NaN;
%   plot(frameTimes,area, 'LineWidth',5); hold on
%   plot(frameTimes,ones(1,numel(area))*medianArea)
%   plot(frameTimes,ones(1,numel(area))*5*medianArea)
%   plot(frameTimes,ones(1,numel(area))*0.2*medianArea)
interpTimes = eyeDataDB.frameTimes(1):1/sr:eyeDataDB.frameTimes(end);
area = interp1(frameTimes(~isnan(area)), area(~isnan(area)), interpTimes);
interpTimesPrev = interpTimes(~isnan(area));
area = area(~isnan(area));
interpTimes = interpTimesPrev(1):1/sr:interpTimesPrev(end);
area = interp1(interpTimesPrev, area, interpTimes);
unfilteredTimes = interpTimes;

% Filter data
d = designFilterLP(1.5, 2, 0.5, 65, sr);
areaFilt = filtfilt(d,area);
%areaFilt = medfilt1(area);

% Remove the filtering artefacts and trim the data to match the ephys data
areaFilt = areaFilt(max([1 filtArtefact]):end-filtArtefact);
areaInd = round(interpTimes*sr);
areaInd = unique([areaInd(areaInd > 0) numel(interpTimes)]);
areaInd = areaInd(max([1 filtArtefact]):end-filtArtefact);
interpTimes = interpTimes(max([1 filtArtefact]):end-filtArtefact);
if interpTimes(end) > numel(ephysSignal)/srRecording
  validInd = floor((numel(ephysSignal)/srRecording)*sr);
  validInd = find(areaInd == validInd);
  areaInd = areaInd(1:validInd);
  areaFilt = areaFilt(1:validInd);
  interpTimes = interpTimes(1:validInd);
end
%plot(interpTimes,areaFilt)

iTime = [];
if ~iscell(period)
  period = {period};
end
for iPeriod = 1:numel(period)
  interpTimesInd1 = abs(interpTimes - period{iPeriod}(1));
  [~, iTime1] = min(interpTimesInd1);
  if period{iPeriod}(2) > interpTimes(end)
    period{iPeriod}(2) = interpTimes(end);
  end
  interpTimesInd2 = abs(interpTimes - period{iPeriod}(2));
  [~, iTime2] = min(interpTimesInd2);
  iTime = [iTime iTime1:iTime2]; %#ok<AGROW>
end
areaFilt = areaFilt(iTime);
interpTimes = interpTimes(iTime);
areaInd = areaInd(iTime);
if plotFigs
  %plot(interpTimes,areaFilt, '.', 'MarkerSize',5)
end

area = interp1(unfilteredTimes, area, interpTimes, 'linear','extrap');