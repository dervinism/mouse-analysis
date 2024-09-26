clear all
load('waveforms.mat')
qualityFile = dir('*.qua.1.mat');
qualityFile = qualityFile.name;
load(qualityFile)
uInds = logical(chanMap(:,3) > 1);
sum(uInds)
chanMap(uInds,1)
sum(unitQ(:,2) > 20 & unitQ(:,6) <= 0.2)
for u = 1:size(unitQ,1)
  if unitQ(u,2) > 20 && unitQ(u,6) <= 0.2
    unitQ(u,1)
  end
end


chanMap1 = chanMap;

load('waveforms.mat')
unitCounter = 0;
for i = 1:numel(chanMap1(:,3))
  if sum(chanMap1(i,3) == chanMap(:,3))
    unitCounter = unitCounter + 1;
  end
end
unitCounter