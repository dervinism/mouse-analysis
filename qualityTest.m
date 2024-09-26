function [qualityUnits, qualityUnitInd, qualityCluDist, qualityRefractCont] = qualityTest(qualitySeries, units, cluDist, refractCont, sortStr)
% [qualityUnits, qualityUnitInd, qualityCluDist, qualityRefractCont] = qualityTest(qualitySeries, units, cluDist, refractCont, sortStr)
%
% Function finds quality units for a given series.
% Input: qualitySeries - The file containing series' unit quality
%                        information ending in .qua.shankID.mat.
%        units - a vector with series units cluster IDs.
%        cluDist - cluster isolation distance (optional), default is 0.2.
%        refractCont - ACG violations as % from ACG baseline (optional),
%                      default is 20.
%        sortStr - if True, will sort units. Otherwise the original order
%                  will be preserved. Default is True.
% cluDist and refractCont are calculated by createQualityFileKilosort.
% Output: qualityUnits - units that pass the quality test.
%         qualityUnitInd - quality unit indices.
%         qualityCluDist - quality unit cluster isolation distances.
%         qualityRefractCont - quality unit ISI violations.

if nargin < 5
  sortStr = true;
end
if nargin < 4
  refractCont = 0.2;
end
if nargin < 3
  cluDist = 20;
end

if iscell(qualitySeries)
  unitQtemp = [];
  for iQS = 1:numel(qualitySeries)
    if ~isempty(qualitySeries{iQS})
      load(qualitySeries{iQS});
      unitQ(:,1) = unitQ(:,1)+1000*(iQS-1);
      unitQtemp = [unitQtemp; unitQ];
    end
  end
  unitQ = unitQtemp;
elseif ~isempty(qualitySeries)
  load(qualitySeries); %#ok<*LOAD>
else
  unitQ = [];
end

qualityUnits = [];
qualityUnitInd = [];
qualityCluDist = [];
qualityRefractCont = [];
if ~isempty(units)
  uCount = 0;
  for u = torow(units)
    uCount = uCount + 1;
    ind = find(u == unitQ(:,1)); %#ok<*IDISVAR,*NODEF>
    if ind > size(unitQ,1)
      ind = ind - size(unitQ,1);
    end
    if (isinf(cluDist) && isinf(refractCont)) || (unitQ(ind,2) >= cluDist && isinf(refractCont))...
        || (isinf(cluDist) && unitQ(ind,6) <= refractCont)...
        || (unitQ(ind,2) >= cluDist && unitQ(ind,6) <= refractCont) % quality check
      qualityUnits = [qualityUnits; u]; %#ok<*AGROW>
      qualityUnitInd = [qualityUnitInd; uCount];
      qualityCluDist = [qualityCluDist; unitQ(ind,2)];
      qualityRefractCont = [qualityRefractCont; unitQ(ind,6)];
    end
  end
end

if sortStr
  [qualityUnits, sortOrder] = sort(qualityUnits);
  qualityUnitInd = qualityUnitInd(sortOrder);
  qualityCluDist = qualityCluDist(sortOrder);
  qualityRefractCont = qualityRefractCont(sortOrder);
end