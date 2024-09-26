function [inds, truncatedData] = determineInds(period, sr, data)
% Determine valid indices for spiking and LFP data

maxInd = size(data,2);
if iscell(period)
  inds = [];
  for i = 1:numel(period)
    inds = [inds round(period{i}(1)*sr):round(period{i}(end)*sr)]; %#ok<*AGROW>
  end
else
  inds = round(period(1)*sr):round(period(end)*sr);
end
inds = inds(inds >= 1 & inds <= maxInd);
truncatedData = data(:,inds);