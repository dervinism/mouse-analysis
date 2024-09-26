function [frameTimes, frameInd] = detectFrames(dataFilename, nChans, sr)
% A helper function for eyeAnalysis script. It finds the times and indices
% of every consecutive peak in the voltage data recorded in an extra Open
% Ephys channel. These peaks correspond to the eye-tracker video frame
% refresh instances.
% Input: dataFilename - a string with the name of the data file (in dat or
%                       bin formats).
%        nChans - a total number of channels in the recording. The voltage
%                 channel is supposed to be the last one.
%        sr - a sampling rate.
% Output: frameTimes and frameInd.

chunkSize = 1000000;

try
  fid = fopen(dataFilename, 'r');
  d = dir(dataFilename);
  nSampsTotal = d.bytes/nChans/2;
catch
  fid = fopen([dataFilename(1:end-3) 'dat'], 'r');
  d = dir([dataFilename(1:end-3) 'dat']);
  nSampsTotal = d.bytes/nChans/2;
end

nChunksTotal = ceil(nSampsTotal/chunkSize);

% Load voltage data
chunkInd = 1;
eyeData = [];
while 1
  fprintf(1, 'chunk %d/%d\n', chunkInd, nChunksTotal);
  dat = fread(fid, [nChans chunkSize], '*int16');
  if ~isempty(dat)
    eyeData = [eyeData dat(end,:)]; %#ok<AGROW>
  else
    break
  end
  
%   plot(eyeData,'r')
%   hold on
%   displacement = 1;
%   threshold = (max(eyeData) - min(eyeData))/3;
%   eyeDataTh = eyeData;
%   eyeDataTh(eyeDataTh<=threshold) = 0;
%   eyeDataTh(eyeDataTh>threshold) = 1;
%   frameOnsets = eyeDataTh - [eyeDataTh(1:displacement) eyeDataTh(1:end-displacement)];
%   frameOnsets(frameOnsets<0.5) = 0;
%   frameOnsets(frameOnsets>0.5) = 1;
%   frameInd = find(logical(int8(round(frameOnsets))));
%   plot(frameInd(2:end),(frameInd(2:end)-frameInd(1:end-1))*10, '.g', 'MarkerSize',20)
%   hold off
  
  chunkInd = chunkInd+1;
end

% Detect frame times
displacement = 1;
peaks = findpeaks(single(eyeData));
if isstruct(peaks)
    peaks = peaks.loc;
end
threshold = double(mean(peaks)/2);
% threshold = maxk(eyeData,100);
% threshold = threshold(end)/2;
eyeDataTh = eyeData;
eyeDataTh(eyeDataTh<=threshold) = 0;
eyeDataTh(eyeDataTh>threshold) = 1;
frameOnsets = eyeDataTh - [eyeDataTh(1:displacement) eyeDataTh(1:end-displacement)];
frameOnsets(frameOnsets<=0.5) = 0;
frameOnsets(frameOnsets>0.5) = 1;

frameInd = find(logical(int8(round(frameOnsets))))';
frameTimes = frameInd*(1/sr)';

frameIntervals = frameInd(2:end)-frameInd(1:end-1);
frameIntervalsNorm = 2*threshold*(frameIntervals./mean(frameIntervals));
[largeFrameIntervals, largeFrameIntervalLocations] = maxk(frameIntervals,20);
plot(1/sr:1/sr:numel(eyeData)/sr, eyeData,'r')
hold on
plot(frameTimes(2:end), frameIntervalsNorm, '.g', 'MarkerSize',20)
% plot(frameTimes(2:end)-frameTimes(1:end-1))
hold off