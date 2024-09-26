function [qualityUnits, qualityUnitInd, qualityCluDist, qualityRefractCont] = qualityTest2(unitMetadata, cluDist, refractCont, sortStr)
% [qualityUnits, qualityUnitInd, qualityCluDist, qualityRefractCont] = qualityTest2(unitMetadata, cluDist, refractCont, sortStr)
%
% Function finds quality units for a given series.
% Input: unitMetadata - unit metadata including unit quality info.
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

if nargin < 4
  sortStr = true;
end
if nargin < 3
  refractCont = 0.2;
end
if nargin < 2
  cluDist = 20;
end

units = unitMetadata(:,1);
unitQ = zeros(numel(units),6);
unitQ(:,1) = units;
unitQ(:,2) = unitMetadata(:,7);
unitQ(:,6) = unitMetadata(:,6);

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