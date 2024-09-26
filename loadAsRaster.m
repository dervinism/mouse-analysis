% A helper function to AnPSD_load and loadAsRaster.
% Inputs: basefilename
%         shanks - which shanks/tetrodes to load, e.g. 2:2:6
%         sr - sampling rate (in Hz), 30000 by default
%         binSize - binSize of the raster (in seconds, 1 by default)
%         opt -
%             .removeLastBin - true by default (because such bin typically spans beyond the end of recording
%             .selectedIntervals - (empty by default) a 2-column matrix with intervals (specified in seconds), only data in the intervals is returned
%             .selectedUnits - if provided only the specified units
%             (instead of all the units) will be included in the output, option can be provided only when loading from no more than one shank(!)
% Outputs: [raster, units] where
%         raster is units by time (2D) matrix
%         units is units-by-2 matrix specifying the shank & id of each unit
% (Note: Very similar to res2raster)
function [raster, units] = loadAsRaster(basefilename, shanks, sr, binSize, opt)
if nargin < 3
  sr = 3e4;
end
if nargin < 4
  binSize = 1;
end
if nargin < 5 || ~isfield(opt, 'removeLastBin')
  opt.removeLastBin = true;
end
if nargin < 5 || ~isfield(opt, 'selectedIntervals')
  opt.selectedIntervals = [];
end

totalUnits = 0;
maxSample = 0;
for s = shanks
  clu = load([basefilename, '.clu.', num2str(s)]);
  clu = clu(2:end);
  totalUnits = length(setdiff(unique(clu), 0)) + totalUnits;
  res = load([basefilename, '.res.', num2str(s)]);
  assert(numel(res) == numel(clu))
  maxSample = max(maxSample, res(end));
end
binSize = round(sr*binSize); % translate to samples
assert(binSize >= 1, 'binSize less than sampling rate allows?!')
if numel(shanks) == 1 && isfield(opt, 'selectedUnits') && ~isempty(opt.selectedUnits)
  totalUnits = length(setdiff(intersect(unique(clu), opt.selectedUnits), 0));
elseif numel(shanks) > 1 && isfield(opt, 'selectedUnits') && ~isempty(opt.selectedUnits)
  error('opt.selectedUnits cannot be used when multiple shanks specfied')
else
  opt.selectedUnits = [];
end

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

% raster is going to be a cell array of matrices
units = zeros(totalUnits, 2);

totalUnits = 0;
for s = shanks
  if ~numel(shanks) == 1
    clu = load([basefilename, '.clu.', num2str(s)]);
    clu = clu(2:end);
    res = load([basefilename, '.res.', num2str(s)]);
  end
  res = res(clu>0);
  clu = clu(clu>0);
  u_clu = unique(clu);
  
  for j = 1:length(u_clu)
    if ~isempty(opt.selectedUnits) && ~any(opt.selectedUnits == u_clu(j))
      continue % this unit is of no interest
    end
    totalUnits = totalUnits + 1;
    raster{totalUnits} = false(1, ceil(maxSample/binSize)*binSize);%logical(sparse(1, ceil(maxSample/binSize)*binSize));
    raster{totalUnits}(res(clu == u_clu(j))) = true;
    units(totalUnits, :) = [s u_clu(j)];
    
    if ~isempty(opt.selectedIntervals)
      raster{totalUnits} = raster{totalUnits}(I,:);
      raster{totalUnits} = raster{totalUnits}(1:floor(size(raster,1)/binSize)*binSize, :); % we have to make sure its an integer number of binSize windows
    end
    if binSize > 1
      raster{totalUnits} = squeeze(sum(reshape(raster{totalUnits}, binSize, [])));
    end
  end % loop on units  
end % loop on shanks
assert(totalUnits == numel(raster))
raster = cell2mat(raster');
if opt.removeLastBin
  raster = raster(:, 1:end-1);
end


