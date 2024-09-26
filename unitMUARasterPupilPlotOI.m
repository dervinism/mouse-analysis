function [fH, percentages] = unitMUARasterPupilPlotOI(dataStruct, series, ops)
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
  [~, ~, area] = determineAreaFromSeries(seriesNameShort);
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
    clear spk_positive spk_negative units_positiveArea units_negativeArea
    spk = concatenateMat(spk, spk_temp);
    clear spk_temp
    rSpearmanUnits = [rSpearmanUnits; rSpearmanUnits_positive'; rSpearmanUnits_negative']; %#ok<*AGROW>
    rSpearmanUnitsPos = [rSpearmanUnitsPos; rSpearmanUnits_positive'];
    rSpearmanUnitsNeg = [rSpearmanUnitsNeg; rSpearmanUnits_negative'];
    units = [units; units_positive; units_negative];
    unitsPos = [unitsPos; units_positive];
    [units_positiveArea{1:numel(units_positive)}] = deal(area);
    unitsPosArea = [unitsPosArea; units_positiveArea'];
    unitsNeg = [unitsNeg; units_negative];
    [units_negativeArea{1:numel(units_negative)}] = deal(area);
    unitsNegArea = [unitsNegArea; units_negativeArea'];
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
    [unitsPosArea{1:numel(unitsPos)}] = deal(area);
    unitsPosArea = unitsPosArea';
    unitsNeg = units_negative;
    [unitsNegArea{1:numel(unitsNeg)}] = deal(area);
    unitsNegArea = unitsNegArea';
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


%% Get population firing rate
if min(size(spkDB)) == 1
  close(fH1)
  return
end

d = designFilterLP(1.5, 2, 0.5, 65, srData);

prPositiveFilt = filtfilt(d,full(prPositive));
prNegativeFilt = filtfilt(d,full(prNegative));
%prPositiveFilt(prPositiveFilt < 0) = 0;
%prNegativeFilt(prNegativeFilt < 0) = 0;
prFilt = sum(concatenateMat(prPositiveFilt, prNegativeFilt),1);
time = (1:numel(prFilt))./srData;
%prPositiveFilt = interp1(time, prPositiveFilt, interpTimes);
%prNegativeFilt = interp1(time, prNegativeFilt, interpTimes);
%prFilt = interp1(time, prFilt, interpTimes);


%% Get unit firing rate
urPositiveFilt = filtfilt(d,full(sum(urPositive,1)));
urNegativeFilt = filtfilt(d,full(sum(urNegative,1)));
%urPositiveFilt = interp1(time, urPositiveFilt, interpTimes);
%urNegativeFilt = interp1(time, urNegativeFilt, interpTimes);
%urPositiveFilt(urPositiveFilt < 0) = 0;
%urNegativeFilt(urNegativeFilt < 0) = 0;


%% Get specific unit firing rates
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
urUnitsArea = [unitsPosArea; unitsNegArea];
urCorr = [rSpearmanUnitsPos; rSpearmanUnitsNeg];
[~, iSort3] = sort(urCorr, 'descend');
[~, iSort4] = sort(rSpearman, 'descend');

urTh = zeros(numel(ops.unitIDsTh), numel(prFilt));
for iUnit = 1:numel(ops.unitIDsTh)
  ind = find(urUnits == ops.unitIDsTh(iUnit));
  disp(['unit ' num2str(ops.unitIDsTh(iUnit)) ' in ' urUnitsArea{ind(1)}...
    '. All inds: ' num2str(ind')]);
  urTh(iUnit,:) = filtfilt(d,full(ur(ind(1),:)));
end

urCx = zeros(numel(ops.unitIDsCx), numel(prFilt));
for iUnit = 1:numel(ops.unitIDsCx)
  ind = find(urUnits == ops.unitIDsCx(iUnit));
  disp(['unit ' num2str(ops.unitIDsCx(iUnit)) ' in ' urUnitsArea{ind(1)}...
    '. All inds: ' num2str(ind')]);
  urCx(iUnit,:) = filtfilt(d,full(ur(ind(1),:)));
end

urHp = zeros(numel(ops.unitIDsHp), numel(prFilt));
for iUnit = 1:numel(ops.unitIDsHp)
  ind = find(urUnits == ops.unitIDsHp(iUnit));
  disp(['unit ' num2str(ops.unitIDsHp(iUnit)) ' in ' urUnitsArea{ind(1)}...
    '. All inds: ' num2str(ind')]);
  urHp(iUnit,:) = filtfilt(d,full(ur(ind(1),:)));
end

urDG = zeros(numel(ops.unitIDsDG), numel(prFilt));
for iUnit = 1:numel(ops.unitIDsDG)
  ind = find(urUnits == ops.unitIDsDG(iUnit));
  disp(['unit ' num2str(ops.unitIDsDG(iUnit)) ' in ' urUnitsArea{ind(1)}...
    '. All inds: ' num2str(ind')]);
  urDG(iUnit,:) = filtfilt(d,full(ur(ind(1),:)));
end


%% Get pupil area
commonPeriod = combinePeriods(period, eyePeriod, srData);
[pupilArea, frameTimes] = pupilFilt(dataStruct.eyeData.([animal '_s' seriesNameShort(1:min([numel(seriesNameShort) 14]))]),...
  srData, prFilt, 2*srData, commonPeriod, srData);
% interpTimes = frameTimes(1:3:end);
% dtInterp = mean(interpTimes(2:end) - interpTimes(1:end-1));


%% Exclude invalid periods
[inds1, spk] = determineInds(commonPeriod, srDataNew, spk);
dsTime = dsTime(inds1);
[inds2, prFilt] = determineInds(commonPeriod, srData, prFilt);
time = time(inds2);

data2plot.y{1} = pupilArea;
data2plot.x{1} = frameTimes;
data2plot.y{2} = spkDB(iSort4,inds1); %spk(iSort3,:);
data2plot.x{2} = dsTime;
data2plot.y{3} = prFilt;
data2plot.x{3} = time;
data2plot.y{4} = urPositiveFilt(inds2);
data2plot.x{4} = time;
data2plot.y{5} = urNegativeFilt(inds2);
data2plot.x{5} = time;
for iUnit = 1:numel(ops.unitIDsTh)
  data2plot.y{numel(data2plot.y)+1} = urTh(iUnit,inds2);
  data2plot.x{numel(data2plot.x)+1} = time;
end
for iUnit = 1:numel(ops.unitIDsCx)
  data2plot.y{numel(data2plot.y)+1} = urCx(iUnit,inds2);
  data2plot.x{numel(data2plot.x)+1} = time;
end
for iUnit = 1:numel(ops.unitIDsHp)
  data2plot.y{numel(data2plot.y)+1} = urHp(iUnit,inds2);
  data2plot.x{numel(data2plot.x)+1} = time;
end
for iUnit = 1:numel(ops.unitIDsDG)
  data2plot.y{numel(data2plot.y)+1} = urDG(iUnit,inds2);
  data2plot.x{numel(data2plot.x)+1} = time;
end

artifact = 2; %sec
artifact = [max([frameTimes(1) dsTime(1) time(1)]+artifact) min([frameTimes(end) dsTime(end) time(end)]-artifact)];
artifact = [max([artifact(1) ops.period(1)]) min([artifact(end) ops.period(end)])];
for iData = 1:numel(data2plot.y)
  data2plot.y{iData} = data2plot.y{iData}(:,data2plot.x{iData} >= artifact(1) & data2plot.x{iData} <= artifact(2));
  data2plot.x{iData} = data2plot.x{iData}(data2plot.x{iData} >= artifact(1) & data2plot.x{iData} <= artifact(2));
end


%% Plot the data
% data2plot.colours = {matlabColours(6),[0 0 0],[0 0 0],matlabColours(15),...
%   matlabColours(10),matlabColours(2),matlabColours(2).*(2/3),matlabColours(11),...
%   matlabColours(11).*(2/3),matlabColours(4),matlabColours(4).*(2/3)};
% data2plot.lineWidths = {2,1,1,1,1,1,1,1,1,1,1};
% data2plot.titles = {'Pupil area','Unit raster','PR','UR+','UR-',...
%   'LGN neuron 1','LGN neuron 2','V1 neuron 1','V2 neuron 2','CA neuron 1','CA neuron 2'};
% data2plot.scaleBars = {[],100,[],[],[],[],[],[],[],[],[]};
data2plot.colours = {matlabColours(6),[0 0 0],[0 0 0],matlabColours(15),...
  matlabColours(10),matlabColours(1),matlabColours(1).*(2/3),matlabColours(13),...
  matlabColours(13).*(2/3),matlabColours(3),matlabColours(3).*(2/3),...
  matlabColours(4),matlabColours(4).*(2/3),matlabColours(12),matlabColours(12).*(2/3)};
data2plot.lineWidths = {2,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
data2plot.titles = {'Pupil area','Unit raster','PR','UR+','UR-',...
  'VB + neuron','VB - neuron','S1 + neuron','S1 - neuron','RSC + neuron','RSC - neuron',...
  'CA + neuron','CA - neuron','DG + neuron','DG - neuron'};
data2plot.scaleBars = {[],100,[],[],[],[],[],[],[],[],[],[],[],[],[]};
fH = plotMultipleDataTypes(data2plot);

ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 25, 4/3, 2, [0.005 0], 'out',...
            'on', 'k', {'Time (s)'}, xlim, xticks,...
            'off', 'w', {}, ylim, []);
set(gca, 'XTickLabel',{'0','100','200','300','400','500','600'});
set(gcf,'color','w');
% figName = ['Multiplot: ' animal '_s' seriesNameShort(1:min([numel(seriesNameShort) 14]))];
figName = ['Multiplot_units+MUAs: ' animal '_s' seriesNameShort(1:min([numel(seriesNameShort) 14]))];
set(gcf,'Name',figName);

if ~exist(ops.outputDir, 'file')
  mkdir(ops.outputDir);
end
figFileName = strrep(figName, ':','_');
figFileName = strrep(figFileName, ' ','_');
figSize = 30;
label = [0.3 2.8];
margin = [0.5 0.1];
% width = ((182.221-31.832)/(279-152))*figSize-label(1)-margin(1);
width = ((182.221-31.832)/(279-133))*figSize-label(1)-margin(1);
height = figSize-label(2)-margin(2);
paperSize = resizeFig(gcf, ax1, width, height, label, margin, 0);
savefig(gcf, [ops.outputDir filesep figFileName '.fig'], 'compact');
exportFig(gcf, [ops.outputDir filesep figFileName '.png'],'-dpng','-r300', paperSize);
exportFig(gcf, [ops.outputDir filesep figFileName '.tif'],'-dtiff','-r1200', paperSize);
exportFig(gcf, [ops.outputDir filesep figFileName '.eps'],'-depsc','-r1200', paperSize);



%% Local function
function fH = plotMultipleDataTypes(data2plot)

traceGap = 1/100;
xMin = data2plot.x{1}(1);
xRange = 0;
for iData = numel(data2plot.y):-1:1
  xMin = min([xMin data2plot.x{iData}(1)]);
  xRange = max([xRange data2plot.x{iData}(end)-data2plot.x{iData}(1)]);
end
for iData = numel(data2plot.y):-1:1
  
  % Baseline
  if iData == numel(data2plot.y)
    baseline = 0;
    %scaleFactor = 1;
    fH = figure; hold on
  else
    %baseline = max(dataY) + max(data2plot.y{iData+1})*traceGap;
    baseline = max(dataY) + traceGap;
    %scaleFactor = max((data2plot.y{iData} - min(data2plot.y{iData})))/...
    %  max((data2plot.y{numel(data2plot.y)} - min(data2plot.y{numel(data2plot.y)})));
  end
  
  % Data
  if min(size(data2plot.y{iData})) > 1
    gridX = data2plot.x{iData};
    gridY = ylim;
    scaleFactor = ((gridY(2)-gridY(1))/size(data2plot.y{iData},1))*(1/2);
    gridY = baseline+(1:size(data2plot.y{iData},1))*scaleFactor;
    p = pcolor(gridX, gridY, flipud(double(logical(full(data2plot.y{iData})))));
    p.EdgeColor = 'none';
    colormap(flipud(gray));
    dataY = gridY(end);
    data2plot.y{iData} = gridY(end)-baseline;
  else
    scaleFactor = 1/...
      max((data2plot.y{iData} - min(data2plot.y{iData})));
    dataY = baseline + scaleFactor.*(data2plot.y{iData} - min(data2plot.y{iData}));
    plot(data2plot.x{iData},dataY, 'Color',data2plot.colours{iData},...
      'LineWidth',data2plot.lineWidths{iData});
  end
  
  % Scale bars
  yRange = max(dataY)-baseline;
  if ~isempty(data2plot.scaleBars{iData})
    plot([xMin+xRange+0.01*xRange xMin+xRange+0.01*xRange],...
      [baseline+0.01*yRange baseline+0.01*yRange+scaleFactor*data2plot.scaleBars{iData}],...
      'Color',data2plot.colours{iData}, 'LineWidth',4);
    xLim = xlim;
    xLim(2) = xMin+xRange+0.01*xRange+0.005*xRange;
  end
  
  % Titles
  if ~isempty(data2plot.titles{iData})
    midY = baseline+yRange/2;
    text(xMin-0.2*xRange,midY, data2plot.titles{iData}, 'Color',data2plot.colours{iData});
  end
end
yLim = ylim;
yLim(2) = max(dataY) + traceGap;
ylim(yLim);
xlim(xLim);
hold off