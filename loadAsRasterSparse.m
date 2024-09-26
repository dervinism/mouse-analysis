% A helper function to AnPSD_load.
% Inputs: dirname
% `       basefilename
%         shanks - which shanks/tetrodes to load, e.g. 2:2:6
%         shCh
%         sr - sampling rate (in Hz), 30000 by default
%         binSize - binSize of the raster (in seconds, 1 by default)
%         opt -
%             .removeLastBin - true by default (because such bin typically spans beyond the end of recording
%             .selectedIntervals - (empty by default) a 2-column matrix with intervals (specified in seconds), only data in the intervals is returned
%             .selectedUnits - if provided only the specified units
%             (instead of all the units) will be included in the output, option can be provided only when loading from no more than one shank(!)
%             .templates - true or false. If true, will load templates
%             instead of clusters. False is default.
% Outputs: [raster, units] where
%         raster is units by time (2D) matrix
%         units is units-by-2 matrix specifying the shank & id of each unit
% (Note: Very similar to loadAsRaster but uses a sparse matrix for output)
function [raster, units, chanMap] = loadAsRasterSparse(dirname, baseFilename, probeFile, shanks, shCh, chOI, sr, binSize, zeroPad, opt)
if nargin < 7
  sr = 3e4;
end
if nargin < 8
  binSize = 1;
end
if nargin < 9
  zeroPad = 0;
end
if nargin < 10 || ~isfield(opt, 'removeLastBin')
  opt.removeLastBin = true;
end
if nargin < 10 || ~isfield(opt, 'selectedIntervals')
  opt.selectedIntervals = [];
end
if nargin < 10 || ~isfield(opt, 'selectedUnits')
  opt.selectedUnits = 'all';
end
if nargin < 10 || ~isfield(opt, 'templates')
  opt.templates = false;
end

if isempty(opt.selectedUnits)
  raster = [];
  units = [];
  chanMap = [];
  return
end

% ESTIMATE OUTPUT DATA DIMENSIONS
totalUnits = 0;
maxSample = 0;
clu = struct([]);
res = struct([]);
resT = struct([]);
tmpl = struct([]);
chanMap0 = struct([]);
for s = shanks
  if exist([baseFilename, '.clu.', num2str(s)], 'file') && exist([baseFilename, '.res.', num2str(s)], 'file')
    clu{s} = load([baseFilename, '.clu.', num2str(s)]);
    clu{s} = clu{s}(2:end);
    uClu = double(unique(clu{s}));
    resT{s} = load([baseFilename, '.res.', num2str(s)]);
    positions = unitPos(dirname, shCh, probeFile);
    unitsOI = positions(logical(sum(positions(:,1) == uClu',2)) & logical(sum(positions(:,2) == chOI,2)),1)';
    chanMap0{s} = positions(logical(sum(positions(:,1) == uClu',2)) & logical(sum(positions(:,2) == chOI,2)),:);
    res{s} = resT{s}(logical(sum(clu{s} == unitsOI,2)));
    clu{s} = clu{s}(logical(sum(clu{s} == unitsOI,2)));
  else
    [clu{s}, res{s}, ~, chanMap0{s}] = resCluFromKilosort(dirname, s, shCh, chOI, probeFile);
    clu{s} = clu{s}(2:end);
  end
  if opt.templates
    if exist([baseFilename, '.tmpl.', num2str(s)], 'file')
      tmpl{s} = load([baseFilename, '.tmpl.', num2str(s)]);
      uClu = double(unique(tmpl{s}));
      positions = unitPos(dirname, shCh, probeFile);
      unitsOI = positions(logical(sum(positions(:,1) == uClu',2)) & logical(sum(positions(:,2) == chOI,2)),1)';
      chanMap0{s} = positions(logical(sum(positions(:,1) == uClu',2)) & logical(sum(positions(:,2) == chOI,2)),:);
      resT{s} = resT{s}(logical(sum(tmpl{s} == unitsOI,2)));
      tmpl{s} = tmpl{s}(logical(sum(tmpl{s} == unitsOI,2)));
    else
      [~, resT{s}, tmpl{s}, chanMap0{s}] = resCluFromKilosort(dirname, s, shCh, chOI, probeFile);
    end
  else
    resT{s} = res{s};
    tmpl{s} = clu{s};
  end
  if ~strcmpi(opt.selectedUnits, 'all')
    res{s} = res{s}(clu{s}>1);
    clu{s} = clu{s}(clu{s}>1);
  end
  resT{s} = resT{s}(tmpl{s}>0);
  tmpl{s} = tmpl{s}(tmpl{s}>0);
  totalUnits = length(setdiff(unique(tmpl{s}), 0)) + totalUnits;
  assert(numel(res{s}) == numel(clu{s}))
  assert(numel(resT{s}) == numel(tmpl{s}))
  if ~isempty(resT{s})
    maxSample = max(maxSample, resT{s}(end));
  end
end
emptyTmplCount = 0;
for s = 1:numel(tmpl)
  if isempty(tmpl{s})
    emptyTmplCount = emptyTmplCount + 1;
  end
  if s == numel(tmpl) && emptyTmplCount == numel(tmpl)
    raster = [];
    units = [];
    chanMap = [];
    return
  end
end
binSize = round(sr*binSize); % translate to samples
assert(binSize >= 1, 'binSize less than sampling rate allows?!')
if numel(shanks) == 1 && isfield(opt, 'selectedUnits') && ~isempty(opt.selectedUnits) && ~strcmpi(opt.selectedUnits, 'all')
  totalUnits = length(setdiff(intersect(unique(tmpl{s}), opt.selectedUnits), 0));
elseif numel(shanks) > 1 && isfield(opt, 'selectedUnits') && ~isempty(opt.selectedUnits) && ~strcmpi(opt.selectedUnits, 'all')
  error('opt.selectedUnits cannot be used when multiple shanks specfied')
elseif ~strcmpi(opt.selectedUnits, 'all')
  opt.selectedUnits = [];
end

% INITIALISE OUTPUT VARIABLES
rasterUnit = zeros(1, ceil(maxSample/binSize)*binSize);
raster = sparse(totalUnits, max(ceil(maxSample/binSize),zeroPad+1));
units = zeros(totalUnits, 2);
for s = 1:numel(chanMap0)
  if ~isempty(chanMap0{s})
    chanMap = zeros(totalUnits, size(chanMap0{s},2));
    break
  end
end

% MARK SELECTED INTERVALS
if ~isempty(opt.selectedIntervals)
  I = false(1, ceil(maxSample/binSize)*binSize);
  for i = 1:size(opt.selectedIntervals, 1)
    if ceil(opt.selectedIntervals(i,1)*sr) > ceil(maxSample/binSize)*binSize
      continue
    end
    I(ceil(opt.selectedIntervals(i,1)*sr):min(ceil(maxSample/binSize)*binSize, ...
      floor(opt.selectedIntervals(i,2)*sr))) = true;
  end
  opt.removeLastBin = false; % we just did it anyway in the line above.
end

% POPULATE OUTPUT VARIABLES
unitCount = 0;
for s = shanks
  u_clu = unique(tmpl{s});
  
  for j = 1:length(u_clu)
    if ~isempty(opt.selectedUnits) && ~any(opt.selectedUnits == u_clu(j)) && ~strcmpi(opt.selectedUnits, 'all')
      continue % this unit is of no interest
    else
      unitCount = unitCount+1;
    end
    rasterUbin = rasterUnit;
    if ~isempty(opt.selectedIntervals)
      resU = resT{s}(tmpl{s} == u_clu(j) && I);
    else
      resU = resT{s}(tmpl{s} == u_clu(j));
    end
    n = numel(resU);
    rasterUbin(resU) = 1:numel(resU);
    resU = unique(resU);
    rasterUbin(resU) = diff([0 rasterUbin(resU)]);
    assert(n == sum(rasterUbin(resU)))
    units(unitCount, :) = [s u_clu(j)];
    if binSize > 1
      jRaster = squeeze(sum(reshape(rasterUbin, binSize, [])));
    else
      jRaster = rasterUbin;
    end
    if zeroPad >= length(jRaster)
      jRaster = [jRaster zeros(1,zeroPad-length(jRaster)+1)];
    end
    raster(unitCount,:) = sparse(jRaster); %#ok<*SPRIX>
    rowChanMap = logical(chanMap0{s}(:,1) == u_clu(j));
    chanMap(unitCount,:) = chanMap0{s}(rowChanMap,:); %#ok<*AGROW>
    assert(u_clu(j) == chanMap0{s}(rowChanMap,1));
  end
end
raster = raster(1:unitCount,:);
units = units(1:unitCount,:);
chanMap = chanMap(1:unitCount,:);

if opt.removeLastBin
  raster = raster(:, 1:end-1);
end
