function [motionFilt, interpTimes, motionInd] = motionFilt(motion, frameTimes, sr, ephysSignal, filtArtefact, period, srRecording)

% figure; plot(frameTimes, motion); hold on

% Resample the data
interpTimes = frameTimes(1):1/sr:frameTimes(end);
motion = interp1(frameTimes, motion, interpTimes);

% Filter data
d = designFilterLP(1.5, 2, 0.5, 65, sr);
motionFilt = filtfilt(d,motion);
%motionFilt = medfilt1(motion);

% Remove the filtering artefacts and trim the data to match the ephys data
motionFilt = motionFilt(max([1 filtArtefact]):end-filtArtefact);
motionInd = round(interpTimes*sr);
motionInd = motionInd(max([1 filtArtefact]):end-filtArtefact);
interpTimes = interpTimes(max([1 filtArtefact]):end-filtArtefact);
if interpTimes(end) > numel(ephysSignal)/srRecording
  validInd = floor((numel(ephysSignal)/srRecording)*sr);
  validInd = find(motionInd == validInd);
  motionInd = motionInd(1:validInd);
  motionFilt = motionFilt(1:validInd);
  interpTimes = interpTimes(1:validInd);
end
% plot(interpTimes,motionFilt)

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
motionFilt = motionFilt(iTime);
interpTimes = interpTimes(iTime);
motionInd = motionInd(iTime);
% plot(interpTimes,motionFilt, '.', 'MarkerSize',5)