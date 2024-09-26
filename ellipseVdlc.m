cleanUp


%% Data files
eyeDataEllipse = 'R:\Neuropix\Shared\Data\M200316_MD\20200316232145\videos_98028_98324_2949_processed.mat';
series = [];
eyeDataDLC = 'R:\Neuropix\Shared\Data\M200316_MD\20200316232145\videos_98028_98324_2949DLC_resnet50_pupilAreaNov23shuffle1_2400.csv';

%% Load and plot elliptic fit data
if ~isempty(eyeDataEllipse)
  if endsWith(eyeDataEllipse, 'processed.mat')
    load(eyeDataEllipse, 'results');
    pupilAreaEllipse = results.area;
    pupilAreaEllipse(isnan(pupilAreaEllipse)) = 0;
    if isfield(results, 'frameTimes')
      frameTimes = results.frameTimes;
      sr = 1/median(frameTimes(2:end)-frameTimes(1:end-1));
    else
      sr = 24;
      frameTimes = (1:numel(pupilAreaEllipse)).*(1/sr);
    end
    d = designFilterLP(1.5, 2, 0.5, 65, sr);
    pupilAreaEllipse_filt = filtfilt(d,pupilAreaEllipse);
    pupilAreaEllipse_filt(pupilAreaEllipse_filt < 0) = 0;
    figH = figure; plot(frameTimes, pupilAreaEllipse_filt)
  elseif endsWith(eyeDataEllipse, '.mat')
    [figH, frameTimes] = plotEyeOrMotion(eyeDataEllipse, series);
    sr = 1/median(frameTimes(2:end)-frameTimes(1:end-1));
  else
    figH = figure;
    obj = VideoReader(eyeDataEllipse);
    sr = obj.FrameRate;
  	frameTimes = (1:(obj.Duration*sr)).*(1/sr);
%     while(hasFrame(obj))
%       frame = readFrame(obj);
%       imshow(frame);
%       title(sprintf('Current Time = %.3f sec', obj.CurrentTime));
%     end
  end
else
  figH = figure;
  sr = 24;
end
hold on


%% Load and plot deeplabcut data
if ~isempty(eyeDataDLC)
  pupilAreaDLC = calcPupilArea(eyeDataDLC);
  d = designFilterLP(1.5, 2, 0.5, 65, sr);
  pupilAreaDLC_filt = filtfilt(d,pupilAreaDLC);
  pupilAreaDLC_filt(pupilAreaDLC_filt < 0) = 0;
  if ~exist('frameTimes','var')
    frameTimes = (1:numel(pupilAreaDLC_filt)).*(1/sr);
  end
  if numel(frameTimes) > numel(pupilAreaDLC_filt)
    warning(['deeplabcut is missing ' num2str(numel(frameTimes) - numel(pupilAreaDLC_filt)) ' video frames']);
%     dt = (frameTimes(end)-frameTimes(1))/(numel(pupilAreaDLC_filt)-1);
%     frameTimes = frameTimes(1):dt:frameTimes(end);
%     plot(frameTimes, pupilAreaDLC_filt)
    plot(frameTimes(1:numel(pupilAreaDLC_filt)), pupilAreaDLC_filt)
  else
    plot(frameTimes, pupilAreaDLC_filt(1:numel(frameTimes)))
  end
end
hold off


%% Tidy the figure
if ~isempty(eyeDataEllipse) && endsWith(eyeDataEllipse, '.mat') && ~isempty(eyeDataDLC)
  legend('Ellipse fit','deeplabcut');
elseif isempty(eyeDataEllipse) && ~isempty(eyeDataDLC)
  legend('deeplabcut');
elseif ~isempty(eyeDataEllipse) && isempty(eyeDataDLC)
  legend('Ellipse fit');
else
  legend('deeplabcut');
end
title('Pupil area');
xlabel('Time (s)');
ylabel('Pupil area (px^2)');