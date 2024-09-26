function fH = rasterPlotOI(dataStruct, series, ops)
% fH = rasterPlotOI(dataStruct, series, ops)
%
% Function produces two raster plots for a given data and series of
% interest. The first plot contains rasters of units and MUAs while the
% second one contains only unit rasters.
% Input: dataStruct
%        series - a cell with series names of interest. If the cell array
%                 has more than one element, rasters from multiple series
%                 are combined into a single raster plot.
%        ops - options structure variable with following fields:
%          normType - 'regular' for a regular raster and 'normalised' for a
%                      normalised one. Each neuron is divided by its mean.
%                      Default is 'regular'.
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
% Output: fH - a cell array of figure handles.


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
if nargin < 3 || ~isfield(ops, 'normType')
  opt.type = 'regular';
else
  opt.type = ops.normType;
end


fH = [];
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
  
  
  %% Get MUA data
  dbStruct = dataStruct.seriesData.(seriesName);
  srData = dbStruct.conf.samplingParams.srData;
  
  if isfield(dataStruct, 'seriesData_positive') && isfield(dataStruct.seriesData_positive, seriesName)
    dbStruct_positive = dataStruct.seriesData_positive.(seriesName);
    if ~isempty(dbStruct_positive.popData.spkDB)
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
dsTime = dsTime - dsTime(1)/2;


%% Plot MUA raster
if min(size(spkDB)) == 1
  return
end
[~, iSort3] = sort(rSpearman, 'descend');
iSort3 = iSort3';
spkDBSorted = spkDB(iSort3,:);
if fraction < 0.5
%   nTopUnits = round(numel(rSpearman)*fraction/2);
%   spkDBSorted = spkDBSorted([1:nTopUnits end-nTopUnits+1:end],:);
  nPositiveUnits = round(numel(spkDBUnitsPos)*fraction);
  nNegativeUnits = round(numel(spkDBUnitsNeg)*fraction);
  spkDBSorted = spkDBSorted([1:nPositiveUnits end-nNegativeUnits+1:end],:);
end

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
if fraction < 0.5
  opt.dividingLine = nPositiveUnits;
else
  opt.dividingLine = numel(spkDBUnitsPos);
end
opt.percentageStr = num2str(round(numel(spkDBUnitsPos)/(numel(spkDBUnitsPos) + numel(spkDBUnitsNeg))*100, 2));
opt.xLim = xLim;
opt.yLabel = 'Units + MUAs';
if strcmp(type, 'regular')
  fH1 = rasterPlot(spkDBSorted, dsTime, opt);
elseif strcmp(type, 'compressed')
  fH1 = rasterPlot(logical(spkDBSorted), dsTime, opt);
end
opt.title = [opt.title ' positivePercentage ' opt.percentageStr];
title(opt.title, 'Interpreter','None');
filename = strrep(opt.title, ' ', '_');
filename = strrep(filename, ':', '_');
filename = strrep(filename, '.', 'p');
if ~exist(outputDir, 'file')
  mkdir(outputDir);
end
filename = [outputDir filesep filename];
savefig(fH1, filename, 'compact');
print(gcf, [filename '.png'],'-dpng','-r300');


%% Plot unit raster
if isempty(spk) || min(size(spk)) == 1
  fH2 = [];
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
  opt.percentageStr = num2str(round(numel(unitsPos)/(numel(unitsPos) + numel(unitsNeg))*100, 2));
  opt.yLabel = 'Units';
  if strcmp(type, 'regular')
    fH2 = rasterPlot(spkSorted, dsTime, opt);
  elseif strcmp(type, 'compressed')
    fH2 = rasterPlot(logical(spkSorted), dsTime, opt);
  end
  opt.title = [opt.title ' positivePercentage ' opt.percentageStr];
  title(opt.title, 'Interpreter','None');
  filename = strrep(opt.title, ' ', '_');
  filename = strrep(filename, ':', '_');
  filename = strrep(filename, '.', 'p');
  filename = [outputDir filesep filename];
  savefig(fH2, filename, 'compact');
  print(gcf, [filename '.png'],'-dpng','-r300');
end

fH = {fH1, fH2};