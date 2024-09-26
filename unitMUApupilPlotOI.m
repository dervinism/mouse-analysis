function [fH, percentages] = unitMUApupilPlotOI(dataStruct, series, ops)
% [fH, percentages] = unitMUApupilPlotOI(dataStruct, series, ops)
%
% Function produces two raster plots for a given data and series of
% interest. The first plot contains rasters of units and MUAs while the
% second one contains only unit rasters.
% Input: dataStruct
%        series - a cell with series names of interest. If the cell array
%                 has more than one element, rasters from multiple series
%                 are combined into a single raster plot.
%        ops - options structure variable with following fields:
%          outputDir - directory for saving figures. By default figures are
%                      saved in the current working directory.
%          srDataNew - raster sampling rate used for producing the figure.
%                      Default is 10 Hz.
%          rasterType - 'compressed' for drawing raster as a logical vector
%                       or 'regular' for drawing actual spike counts.
%                       Default is 'regular'.
%          xLim - range. Default is full range.
%          fraction - fraction of top-most and bottom-most units to be
%                     displayed For example, if ops.fraction = 0.1, 5% of
%                     top-most and 5% of bottom-most units will only be
%                     displayed. Default is 1.
%          calcRanksOnly - true or false. If true, calculates positively
%                          pupil area-correlated unit and MUA fractions
%                          only without producing figures. If false, the
%                          function produces figures and looks for fraction
%                          estimates file named positivityRanks.mat in the
%                          output folder (if provided) or in the present
%                          working directory. The default option is false.
%          period - restrict figure traces to a certain time period.
% Output: fH - a cell array of figure handles.
%         percentages - an array of positively pupil area size-correlated
%                       units and MUAs and units only.


%% Initialise input
if nargin < 3 || ~isfield(ops, 'outputDir')
  outputDir = pwd;
else
  outputDir = ops.outputDir;
end
if nargin < 3 || ~isfield(ops, 'srDataNew')
  srDataNew = 10;
else
  srDataNew = ops.srDataNew;
end
if nargin < 3 || ~isfield(ops, 'rasterType')
  type = 'regular';
else
  type = ops.rasterType;
end
if nargin < 3 || ~isfield(ops, 'xLim')
  xLim = [];
else
  xLim = ops.xLim;
end
if nargin < 3 || ~isfield(ops, 'fraction')
  fraction = 1;
else
  fraction = ops.fraction;
end
if nargin < 3 || ~isfield(ops, 'calcRanksOnly')
  calcRanksOnly = false;
else
  calcRanksOnly = ops.calcRanksOnly;
end
if nargin < 3 || ~isfield(ops, 'period')
  ops.period = [0 99999];
end
plotRasters = false;


fH = {};
percentages = [];
fnsData = fieldnames(dataStruct.seriesData);
animal = dataStruct.seriesData.(fnsData{1}).dbSeries.animal;
for dbCount = 1:numel(series)
  clear dsTimeNew
  if strcmpi(animal, series{dbCount}(1:numel(animal)))
    seriesName = series{dbCount};
  else
    seriesName = [animal '_s' series{dbCount}];
  end
  if ~isfield(dataStruct.seriesData, seriesName)
    disp(['Series ' seriesName ' does not exist. Skipping to the next series...']);
    continue
  end
  disp(['Processing series ' seriesName]);
  
  
  %% Get acceptable electrophysiological and eye data recording periods
  period = dataStruct.seriesData.(seriesName).dbSeries.period;
  seriesNameShort = seriesFromEntry(seriesName);
  eyePeriod = dataStruct.eyeData.([animal '_s' seriesNameShort(1:min([numel(seriesNameShort) 14]))]).period;
  
  
  %% Get MUA data
  dbStruct = dataStruct.seriesData.(seriesName);
  srData = dbStruct.conf.samplingParams.srData;
  
  if isfield(dataStruct, 'seriesData_positive') && isfield(dataStruct.seriesData_positive, seriesName)
    dbStruct_positive = dataStruct.seriesData_positive.(seriesName);
    if ~isempty(dbStruct_positive.popData.spkDB)
      if exist('prPositive', 'var')
        prPositive = sum(concatenateMat(prPositive, sum(dbStruct_positive.popData.spkDB,1)),1);
      else
        prPositive = sum(dbStruct_positive.popData.spkDB,1);
      end
      [spkDB_positive, dsTimeNew] = downsampleRasterMatrix(full(dbStruct_positive.popData.spkDB), srData, srDataNew);
      if min(size(spkDB_positive)) == 1
        spkDB_positive = torow(spkDB_positive);
      end
      spkDBUnits_positive = dbStruct_positive.popData.spkDB_units(:,end);
      rSpearman_positive = dbStruct_positive.popData.rSpearman;
    else
      spkDB_positive = []; spkDBUnits_positive = []; rSpearman_positive = [];
    end
  else
    spkDB_positive = []; spkDBUnits_positive = []; rSpearman_positive = [];
  end
  
  if isfield(dataStruct, 'seriesData_negative') && isfield(dataStruct.seriesData_negative, seriesName)
    dbStruct_negative = dataStruct.seriesData_negative.(seriesName);
    if ~isempty(dbStruct_negative.popData.spkDB)
      if exist('prNegative', 'var')
        prNegative = sum(concatenateMat(prNegative, sum(dbStruct_negative.popData.spkDB,1)),1);
      else
        prNegative = sum(dbStruct_negative.popData.spkDB,1);
      end
      [spkDB_negative, dsTimeNew] = downsampleRasterMatrix(full(dbStruct_negative.popData.spkDB), srData, srDataNew);
      if min(size(spkDB_negative)) == 1
        spkDB_negative = torow(spkDB_negative);
      end
      spkDBUnits_negative = dbStruct_negative.popData.spkDB_units(:,end);
      rSpearman_negative = dbStruct_negative.popData.rSpearman;
    else
      spkDB_negative = []; spkDBUnits_negative = []; rSpearman_negative = [];
    end
  else
    spkDB_negative = []; spkDBUnits_negative = []; rSpearman_negative = [];
  end
  
  if (isempty(spkDB_positive) && isempty(spkDB_negative)) || ~exist('dsTimeNew', 'var')
    continue
  end
  if exist('spkDB', 'var')
    spkDB_temp = concatenateMat(spkDB_positive, spkDB_negative);
    clear spkDB_positive spkDB_negative
    spkDB = concatenateMat(spkDB, spkDB_temp);
    clear spkDB_temp
    rSpearman = [rSpearman; rSpearman_positive'; rSpearman_negative']; %#ok<*AGROW>
    spkDBUnits = [spkDBUnits; spkDBUnits_positive; spkDBUnits_negative];
    spkDBUnitsPos = [spkDBUnitsPos; spkDBUnits_positive];
    spkDBUnitsNeg = [spkDBUnitsNeg; spkDBUnits_negative];
    if numel(dsTimeNew) > numel(dsTime)
      dsTime = dsTimeNew;
    end
  else
    spkDB = concatenateMat(spkDB_positive, spkDB_negative);
    clear spkDB_positive spkDB_negative
    rSpearman = [rSpearman_positive'; rSpearman_negative'];
    spkDBUnits = [spkDBUnits_positive; spkDBUnits_negative];
    spkDBUnitsPos = spkDBUnits_positive;
    spkDBUnitsNeg = spkDBUnits_negative;
    dsTime = dsTimeNew;
  end
  
  
  %% Get unit data
  %shankStruct = dbStruct.shankData.shank1;
  
  if isfield(dataStruct, 'seriesData_positive') && isfield(dataStruct.seriesData_positive, seriesName)
    shankStruct_positive = dbStruct_positive.shankData.shank1;
    if ~isempty(shankStruct_positive.spk)
      if exist('urPositive', 'var')
        urPositive = concatenateMat(urPositive, shankStruct_positive.spk);
      else
        urPositive = shankStruct_positive.spk;
      end
      [spk_positive, dsTimeNew] = downsampleRasterMatrix(full(shankStruct_positive.spk), srData, srDataNew);
      if min(size(spk_positive)) == 1
        spk_positive = torow(spk_positive);
      end
      units_positive = shankStruct_positive.units(:,end);
      rSpearmanUnits_positive = rSpearman_positive(ismember(spkDBUnits_positive, units_positive));
    else
      spk_positive = []; units_positive = []; rSpearmanUnits_positive = [];
    end
  else
    spk_positive = []; units_positive = []; rSpearmanUnits_positive = [];
  end
  
  if isfield(dataStruct, 'seriesData_negative') && isfield(dataStruct.seriesData_negative, seriesName)
    shankStruct_negative = dbStruct_negative.shankData.shank1;
    if ~isempty(shankStruct_negative.spk)
      if exist('urNegative', 'var')
        urNegative = concatenateMat(urNegative, shankStruct_negative.spk);
      else
        urNegative = shankStruct_negative.spk;
      end
      [spk_negative, dsTimeNew] = downsampleRasterMatrix(full(shankStruct_negative.spk), srData, srDataNew);
      if min(size(spk_negative)) == 1
        spk_negative = torow(spk_negative);
      end
      units_negative = shankStruct_negative.units(:,end);
      rSpearmanUnits_negative = rSpearman_negative(ismember(spkDBUnits_negative, units_negative));
    else
      spk_negative = []; units_negative = []; rSpearmanUnits_negative = [];
    end
  else
    spk_negative = []; units_negative = []; rSpearmanUnits_negative = [];
  end
  
  if ~exist('dsTimeNew', 'var')
    continue
  end
  if exist('spk', 'var')
    spk_temp = concatenateMat(spk_positive, spk_negative);
    clear spk_positive spk_negative
    spk = concatenateMat(spk, spk_temp);
    clear spk_temp
    rSpearmanUnits = [rSpearmanUnits; rSpearmanUnits_positive'; rSpearmanUnits_negative']; %#ok<*AGROW>
    rSpearmanUnitsPos = [rSpearmanUnitsPos; rSpearmanUnits_positive'];
    rSpearmanUnitsNeg = [rSpearmanUnitsNeg; rSpearmanUnits_negative'];
    units = [units; units_positive; units_negative];
    unitsPos = [unitsPos; units_positive];
    unitsNeg = [unitsNeg; units_negative];
    if numel(dsTimeNew) > numel(dsTime)
      dsTime = dsTimeNew;
    end
  else
    spk = concatenateMat(spk_positive, spk_negative);
    clear spk_positive spk_negative
    rSpearmanUnits = [rSpearmanUnits_positive'; rSpearmanUnits_negative'];
    rSpearmanUnitsPos = rSpearmanUnits_positive';
    rSpearmanUnitsNeg = rSpearmanUnits_negative';
    units = [units_positive; units_negative];
    unitsPos = units_positive;
    unitsNeg = units_negative;
    dsTime = dsTimeNew;
  end
end
if ~exist('dsTime', 'var')
  return
end
assert(numel(rSpearman) == size(spkDB,1));
assert(numel(rSpearmanUnits) == size(spk,1));
dsTime = dsTime - dsTime(1)/2; %#ok<*NASGU>


%% Calculate positively pupil area-correlated fractions
percentageMUAs = round(numel(spkDBUnitsPos)/(numel(spkDBUnitsPos) + numel(spkDBUnitsNeg))*100, 2);
percentageUnits = round(numel(unitsPos)/(numel(unitsPos) + numel(unitsNeg))*100, 2);
percentages = [percentageMUAs percentageUnits];
if calcRanksOnly
  return
else
  load([outputDir filesep 'positivityRanks.mat']); %#ok<*LOAD>
  seriesNameShort = seriesFromEntry(seriesName);
  rank = positivityRanks.(animal).([animal '_s' seriesNameShort(1:min([14 numel(seriesNameShort)]))]).rankMUAs;
  rank = positivityRanks.(animal).([animal '_s' seriesNameShort(1:min([14 numel(seriesNameShort)]))]).rankUnits;
  rankStr = num2str(rank);
  if numel(rankStr) == 1
    rankStr = ['0' rankStr];
  elseif numel(rankStr) == 0
    error('Empty rank variable');
  end
end

  
%% Plot pupil area
[fH1, frameTimes, pupilArea] = plotEyeOrMotion(dataStruct, series{dbCount}, [], 'pupil', [], pwd, true, false);
grid on
grid minor
xLim = xlim;
xAxisLength = xLim(2)-xLim(1);
rightLabelXPos = xLim(2) + xAxisLength*0.001;
text(rightLabelXPos, mean([min(pupilArea) max(pupilArea)]), 'pupil', 'Color',matlabColours(1));

interpTimes = frameTimes(1:3:end);
dtInterp = mean(interpTimes(2:end) - interpTimes(1:end-1));


%% Mark common excluded times
commonPeriod = combinePeriods(period, eyePeriod, srData);
if isempty(commonPeriod)
  close(fH1)
  return
end

hold on
plot(xlim, [0 0], 'r', 'LineWidth',2)
if iscell(commonPeriod)
  for iPeriod = 1:numel(commonPeriod)
    if iPeriod == numel(commonPeriod)
      plot([commonPeriod{iPeriod}(1) min([commonPeriod{iPeriod}(2) frameTimes(end)])], [0 0], 'k', 'LineWidth',2)
    else
      plot(commonPeriod{iPeriod}, [0 0], 'k', 'LineWidth',2)
    end
  end
else
  plot([commonPeriod(1) min([commonPeriod(2) frameTimes(end)])], [0 0], 'k', 'LineWidth',2)
end


%% Plot population firing rate
if min(size(spkDB)) == 1
  close(fH1)
  return
end

d = designFilterLP(1.5, 2, 0.5, 65, srData);
artifact = round((2000/srData)/dtInterp);

prPositiveFilt = filtfilt(d,full(prPositive));
prNegativeFilt = filtfilt(d,full(prNegative));
%prPositiveFilt(prPositiveFilt < 0) = 0;
%prNegativeFilt(prNegativeFilt < 0) = 0;
prFilt = sum(concatenateMat(prPositiveFilt, prNegativeFilt),1);
time = (1:numel(prFilt))./srData;
prPositiveFilt = interp1(time, prPositiveFilt, interpTimes);
prNegativeFilt = interp1(time, prNegativeFilt, interpTimes);
prFilt = interp1(time, prFilt, interpTimes);
traceGap = 1/100;

% Population firing rate
shiftFactor = max(pupilArea) + max(pupilArea)*traceGap;
scaleFactor = max(pupilArea)/max(prFilt(artifact+1:end-artifact));
data2plot = shiftFactor + scaleFactor.*(prFilt(artifact+1:end-artifact)) -...
  min(scaleFactor.*(prFilt(artifact+1:end-artifact)));
plot(interpTimes(artifact+1:end-artifact), data2plot, 'k', 'LineWidth',1);
text(rightLabelXPos, mean([min(data2plot) max(data2plot)]), 'PR', 'Color','k');

% Positive population firing rate
shiftFactor = max(data2plot) + max(data2plot)*traceGap;
scaleFactor = max(pupilArea)/max(prPositiveFilt(artifact+1:end-artifact));
data2plot = shiftFactor + scaleFactor.*(prPositiveFilt(artifact+1:end-artifact)) -...
  min(scaleFactor.*(prPositiveFilt(artifact+1:end-artifact)));
plot(interpTimes(artifact+1:end-artifact), data2plot, 'Color',matlabColours(9) , 'LineWidth',1);
text(rightLabelXPos, mean([min(data2plot) max(data2plot)]), 'PR+', 'Color',matlabColours(9));

% Negative population firing rate
shiftFactor = max(data2plot) + max(data2plot)*traceGap;
scaleFactor = max(pupilArea)/max(prNegativeFilt(artifact+1:end-artifact));
data2plot = shiftFactor + scaleFactor.*(prNegativeFilt(artifact+1:end-artifact)) -...
  min(scaleFactor.*(prNegativeFilt(artifact+1:end-artifact)));
plot(interpTimes(artifact+1:end-artifact), data2plot, 'Color',matlabColours(7), 'LineWidth',1);
text(rightLabelXPos, mean([min(data2plot) max(data2plot)]), 'PR-', 'Color',matlabColours(7));


%% Plot combined unit firing rate
urPositiveFilt = filtfilt(d,full(sum(urPositive,1)));
urNegativeFilt = filtfilt(d,full(sum(urNegative,1)));
urPositiveFilt = interp1(time, urPositiveFilt, interpTimes);
urNegativeFilt = interp1(time, urNegativeFilt, interpTimes);
%urPositiveFilt(urPositiveFilt < 0) = 0;
%urNegativeFilt(urNegativeFilt < 0) = 0;

% Positive combined unit firing rate
shiftFactor = max(data2plot) + max(data2plot)*traceGap;
scaleFactor = max(pupilArea)/max(urPositiveFilt(artifact+1:end-artifact));
data2plot = shiftFactor + scaleFactor.*(urPositiveFilt(artifact+1:end-artifact)) -...
  min(scaleFactor.*(urPositiveFilt(artifact+1:end-artifact)));
plot(interpTimes(artifact+1:end-artifact), data2plot, 'g', 'LineWidth',1);
text(rightLabelXPos, mean([min(data2plot) max(data2plot)]), 'UR+', 'Color','g');

% Negative combined unit firing rate
shiftFactor = max(data2plot) + max(data2plot)*traceGap;
scaleFactor = max(pupilArea)/max(urNegativeFilt(artifact+1:end-artifact));
data2plot = shiftFactor + scaleFactor.*(urNegativeFilt(artifact+1:end-artifact)) -...
  min(scaleFactor.*(urNegativeFilt(artifact+1:end-artifact)));
plot(interpTimes(artifact+1:end-artifact), data2plot, 'r', 'LineWidth',1);
text(rightLabelXPos, mean([min(data2plot) max(data2plot)]), 'UR-', 'Color','r');


%% Plot specific unit firing rate
if numel(series) > 1
  filenameInit = 'Pupil area and spiking in areas:';
  for dbCount = 1:numel(series)
    if strcmpi(animal, series{dbCount}(1:numel(animal)))
      seriesNameShort = seriesFromEntry(series{dbCount});
    else
      seriesNameShort = series{dbCount};
    end
    [~, ~, area] = determineAreaFromSeries(seriesNameShort);
    if dbCount == numel(series)
      filenameInit = [filenameInit ' ' area];
    else
      filenameInit = [filenameInit ' ' area ' +'];
    end
  end
  filenameInit = [filenameInit ' (' animal '_s' seriesNameShort(1:14) ')'];
else
  if strcmpi(animal, series{1}(1:numel(animal)))
    seriesNameShort = seriesFromEntry(series{1});
  else
    seriesNameShort = series{1};
  end
  [~, ~, area] = determineAreaFromSeries(seriesNameShort);
  filenameInit = ['Pupil area and spiking in ' area ' (' animal '_s' seriesNameShort ')'];
end
filenameInit = [filenameInit ' positivePercentMUAs ' num2str(percentageMUAs)];
filenameInit = [filenameInit ' positivePercentUnits ' num2str(percentageUnits)];
if ~exist(outputDir, 'file')
  mkdir(outputDir);
end

ur = [urPositive; urNegative];
urUnits = [unitsPos; unitsNeg];
urCorr = [rSpearmanUnitsPos; rSpearmanUnitsNeg];
[~, iSort3] = sort(urCorr, 'descend');
shiftFactor = max(data2plot) + max(data2plot)*traceGap;
unitCounter = 0;
for iUnit = torow(iSort3)
  unitCounter = unitCounter + 1;
  
  % unitActivity = full(ur(iUnit,:));
  % scaleFactor = max(pupilArea)/max(unitActivity);
  % data2plot = shiftFactor + scaleFactor.*unitActivity - min(scaleFactor.*unitActivity);
  % plot(interpTimes, data2plot, 'k', 'LineWidth',1);
  % text(rightLabelXPos, mean(unitActivity), 'unit x', 'Color','m');
  
  unitActivity = ur(iUnit,:);
  unitActivityFilt = filtfilt(d,full(unitActivity));
  %unitActivityFilt(unitActivityFilt < 0) = 0;
  unitActivityFilt = interp1(time, unitActivityFilt, interpTimes);
  scaleFactor = max(pupilArea)/max(unitActivityFilt(artifact+1:end-artifact));
  data2plot = shiftFactor + scaleFactor.*(unitActivityFilt(artifact+1:end-artifact)) -...
    min(scaleFactor.*(unitActivityFilt(artifact+1:end-artifact)));
  if urCorr(iUnit) >= 0
    p = plot(interpTimes(artifact+1:end-artifact), data2plot, 'Color',matlabColours(13), 'LineWidth',1);
    t = text(rightLabelXPos, mean([min(data2plot) max(data2plot)]), ['+u' num2str(urUnits(iUnit))], 'Color',matlabColours(13));
  elseif urCorr(iUnit) < 0
    p = plot(interpTimes(artifact+1:end-artifact), data2plot, 'm', 'LineWidth',1);
    t = text(rightLabelXPos, mean([min(data2plot) max(data2plot)]), ['-u' num2str(urUnits(iUnit))], 'Color','m');
  else
    error('Spearman correlation coeffiecient is NaN');
  end
  
  xlim(ops.period);
  ylim([0 max(data2plot) + max(data2plot)*traceGap]);
  
  
  %% Tidy and save the figure
  set(gca, 'YTick',[]);
  set(gca, 'YTicklabel',[]);
  ylabel('Normalised activity');
  title('Pupil area and spiking');
  if urCorr(iUnit) >= 0
    filename = [filenameInit '  unit rank ' num2str(unitCounter) 'positive  unit id ' num2str(urUnits(iUnit))];
  else
    filename = [filenameInit '  unit rank ' num2str(unitCounter) 'negative  unit id ' num2str(urUnits(iUnit))];
  end
  filename = [rankStr ' ' filename];
  set(fH1, 'Name',filename);
  filename = strrep(filename, ' ', '_');
  filename = strrep(filename, ':', '_');
  filename = strrep(filename, '.', 'p');
  filename = [outputDir filesep filename];
  %savefig(fH1, filename, 'compact');
  print(fH1, [filename '.png'],'-dpng','-r300');
  
  delete(p);
  delete(t);
end
hold off
grid minor
grid off


if plotRasters
  
  
  %% Reduce rasters
  [~, iSort3] = sort(rSpearman, 'descend'); %#ok<*UNRCH>
  iSort3 = iSort3';
  spkDBSorted = spkDB(iSort3,:);
  if fraction < 0.5
    %   nTopUnits = round(numel(rSpearman)*fraction/2);
    %   spkDBSorted = spkDBSorted([1:nTopUnits end-nTopUnits+1:end],:);
    nPositiveUnits = round(numel(spkDBUnitsPos)*fraction);
    nNegativeUnits = round(numel(spkDBUnitsNeg)*fraction);
    spkDBSorted = spkDBSorted([1:nPositiveUnits end-nNegativeUnits+1:end],:);
  end
  
  
  %% Plot unit+MUAs raster
  if numel(series) > 1
    opt.title = 'Population spiking in areas:';
    for dbCount = 1:numel(series)
      if strcmpi(animal, series{dbCount}(1:numel(animal)))
        seriesNameShort = seriesFromEntry(series{dbCount});
      else
        seriesNameShort = series{dbCount};
      end
      [~, ~, area] = determineAreaFromSeries(seriesNameShort);
      if dbCount == numel(series)
        opt.title = [opt.title ' ' area];
      else
        opt.title = [opt.title ' ' area ' +'];
      end
    end
    opt.title = [opt.title ' (' animal '_s' seriesNameShort(1:14) ')'];
  else
    if strcmpi(animal, series{1}(1:numel(animal)))
      seriesNameShort = seriesFromEntry(series{1});
    else
      seriesNameShort = series{1};
    end
    [~, ~, area] = determineAreaFromSeries(seriesNameShort);
    opt.title = ['Population spiking in ' area ' (' animal '_s' seriesNameShort ')'];
  end
  opt.title = [rankStr ' ' opt.title];
  if fraction < 0.5
    opt.dividingLine = nPositiveUnits;
  else
    opt.dividingLine = numel(spkDBUnitsPos);
  end
  opt.percentageStr = num2str(percentageMUAs);
  opt.xLim = xLim;
  opt.yLabel = 'Units + MUAs';
  if strcmp(type, 'regular')
    fH2 = rasterPlot(spkDBSorted, dsTime, opt);
  elseif strcmp(type, 'compressed')
    fH2 = rasterPlot(logical(spkDBSorted), dsTime, opt);
  end
  opt.title = [opt.title ' positivePercentage ' opt.percentageStr];
  title(opt.title, 'Interpreter','None');
  filename = strrep(opt.title, ' ', '_');
  filename = strrep(filename, ':', '_');
  filename = strrep(filename, '.', 'p');
  filename = [outputDir filesep filename];
  savefig(fH2, filename, 'compact');
  print(gcf, [filename '.png'],'-dpng','-r300');
  
  
  %% Plot unit raster
  if isempty(spk) || min(size(spk)) == 1
    fH3 = [];
  else
    [~, iSort3] = sort(rSpearmanUnits, 'descend');
    iSort3 = iSort3';
    spkSorted = spk(iSort3,:);
    if fraction < 0.5
      %     nTopUnits = round(numel(rSpearmanUnits)*fraction/2);
      %     spkSorted = spkSorted([1:nTopUnits end-nTopUnits+1:end],:);
      nPositiveUnits = round(numel(unitsPos)*fraction);
      nNegativeUnits = round(numel(unitsNeg)*fraction);
      spkSorted = spkSorted([1:nPositiveUnits end-nNegativeUnits+1:end],:);
    end
    
    opt.title = strrep(opt.title, 'Population', 'Unit');
    if fraction < 0.5
      opt.dividingLine = nPositiveUnits;
    else
      opt.dividingLine = numel(unitsPos);
    end
    opt.percentageStr = num2str(percentageUnits);
    opt.yLabel = 'Units';
    if strcmp(type, 'regular')
      fH3 = rasterPlot(spkSorted, dsTime, opt);
    elseif strcmp(type, 'compressed')
      fH3 = rasterPlot(logical(spkSorted), dsTime, opt);
    end
    opt.title = [opt.title ' positivePercentage ' opt.percentageStr];
    title(opt.title, 'Interpreter','None');
    filename = strrep(opt.title, ' ', '_');
    filename = strrep(filename, ':', '_');
    filename = strrep(filename, '.', 'p');
    filename = [outputDir filesep filename];
    savefig(fH3, filename, 'compact');
    print(gcf, [filename '.png'],'-dpng','-r300');
  end
  
  fH = {fH1, fH2, fH3};
else
  fH = fH1;
end