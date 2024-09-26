function [eMap, hDim, vDim] = electrodeMap(probe, sites, connectedSites)
% A helper function to AnPSD_segs.

siteIDs = sites;
siteIDs(~connectedSites) = NaN;
if strcmp(probe, 'A1x32_Edge_5mm_20_177')
  hDim = 1;
  vDim = numel(sites);
  eMap = siteIDs';
elseif strcmp(probe, 'imecP3opt3') || strcmp(probe, 'Neuropixels')
  hDim = 4;
  vDim = round(numel(sites)/2);
  basicMotive = [1 0 1 0 0 1 0 1];
  eMapExt = repmat(basicMotive,1,vDim/2);
  eMap1(logical(eMapExt)) = siteIDs';
  eMap1(~logical(eMapExt)) = NaN;
  eMap1 = reshape(eMap1,hDim,vDim);
  eMap1 = rot90(eMap1',2);
  eMap2 = eMap1(:,[2 1 4 3]);
%   eMap3 = fliplr(eMap1);
%   eMap4 = fliplr(eMap2);
  eMap = eMap2;
end