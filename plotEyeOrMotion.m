function [figH, frameTimes, pupilSignal] = plotEyeOrMotion(filename, series, period, dataType, ylimits, outputDir, margins, saveFig, produceFig)
% [figH, frameTimes] = plotEyeOrMotion(filename, series, period, dataType, ylimits, outputDir, margins, saveFig)
% This function plots the pupil area signal for a given time interval
% Input: filename
%        series
%        period
%        dataType - 'pupil' or 'motion'
%        ylimits
%        outputDir
%        margins - true or false
%        saveFig - true or false
%        produceFig - true or false
% Output: figH - figure handle.
%         frameTimes
%         pupilSignal - filtered pupil area

if nargin < 9
  produceFig = true;
end
if nargin < 8
  saveFig = true;
end
if nargin < 7
  margins = false;
end
if nargin < 6
  outputDir = pwd;
end
if nargin < 5
  ylimits = [];
end
if nargin < 4
  dataType = 'pupil';
end
if nargin < 3
  period = [];
end
if isstring(filename)
  load(filename); %#ok<*LOAD>
elseif isstruct(filename)
  dataStruct = filename;
else
  error('Unknown input variable filename type');
end

[seriesName, animal] = seriesFromEntry(series);
if strcmpi(dataType, 'pupil')
  pupilSignal = dataStruct.eyeData.([animal '_s' seriesName(1:min([numel(seriesName) 14]))]).pupilArea;
  frameTimes = dataStruct.eyeData.([animal '_s' seriesName(1:min([numel(seriesName) 14]))]).frameTimes;
  totalPeriod = dataStruct.eyeData.([animal '_s' seriesName(1:min([numel(seriesName) 14]))]).period;
elseif strcmpi(dataType, 'motion')
  pupilSignal = dataStruct.motionData.([animal '_s' seriesName(1:min([numel(seriesName) 14]))]).s;
  frameTimes = dataStruct.motionData.([animal '_s' seriesName(1:min([numel(seriesName) 14]))]).frameTimes;
  totalPeriod = dataStruct.eyeData.([animal '_s' seriesName(1:min([numel(seriesName) 14]))]).period;
end
pupilSignal(isnan(pupilSignal)) = 0;
pupilSignal(pupilSignal > 5*median(pupilSignal)) = 5*median(pupilSignal);
sr = 1/median(frameTimes(2:end)-frameTimes(1:end-1));
try
  d = designFilterLP(1.5, 2, 0.5, 65, sr);
catch
  figH = []; frameTimes = [];
  return
end
pupilSignal = filtfilt(d,pupilSignal);
pupilSignal(pupilSignal < 0) = 0;

if produceFig
  figH = figure;
  plot(frameTimes,pupilSignal, 'LineWidth',2);
  if ~isempty(period)
    xlim(period);
  else
    period = xlim;
  end
  if ~isempty(ylimits)
    ylim(ylimits);
  else
    ylimits = ylim;
  end
  
  if exist('totalPeriod', 'var')
    hold on
    plot(xlim, [0 0], 'r', 'LineWidth',2)
    if iscell(totalPeriod)
      for iPeriod = 1:numel(totalPeriod)
        if iPeriod == numel(totalPeriod)
          plot([totalPeriod{iPeriod}(1) min([totalPeriod{iPeriod}(2) frameTimes(end)])], [0 0], 'k', 'LineWidth',2)
        else
          plot(totalPeriod{iPeriod}, [0 0], 'k', 'LineWidth',2)
        end
      end
    else
      plot([totalPeriod(1) min([totalPeriod(2) frameTimes(end)])], [0 0], 'k', 'LineWidth',2)
    end
    hold off
  end
  
  if strcmpi(dataType, 'pupil')
    ax1 = axesProperties('Pupil area', 1, 'normal', 'off', 'w', 'Calibri', 25, 4/3, 2, [0.005 0], 'out',...
      'on', 'k', {'Time (s)'}, period, xticks,...
      'on', 'k', {'Pupil area (a.u.)'}, ylimits, yticks);
  elseif strcmpi(dataType, 'motion')
    ax1 = axesProperties('Total movement', 1, 'normal', 'off', 'w', 'Calibri', 25, 4/3, 2, [0.005 0], 'out',...
      'on', 'k', {'Time (s)'}, period, xticks,...
      'on', 'k', {'Total movement (a.u.)'}, ylimits, yticks);
  end
  
  if margins
    label = [4 3];
    margin = [1.25 1.25];
  else
    label = [0 0];
    margin = [0 0];
  end
  width = 2*15-label(1)-margin(1);
  height = 1*15-label(2)-margin(2);
  paperSize = resizeFig(figH, ax1, width, height, label, margin, 0);
  if saveFig
    [seriesName, animal] = seriesFromEntry(series);
    if strcmpi(dataType, 'pupil')
      figFileName = [animal '_s' seriesName(1:min([numel(seriesName) 14])) '_pupilSignal'];
    elseif strcmpi(dataType, 'motion')
      figFileName = [animal '_s' seriesName(1:min([numel(seriesName) 14])) '_motionSignal'];
    end
    set(figH, 'Name',figFileName);
    
    if ~exist(outputDir, 'file')
      mkdir(outputDir);
    end
    hgsave(figH, [outputDir filesep figFileName '.fig']);
    exportFig(figH, [outputDir filesep figFileName '.png'],'-dpng','-r300', paperSize);
    close(figH);
  end
else
  figH = [];
end