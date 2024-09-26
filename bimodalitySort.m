%load('R:\CSN\md406\runNSG\runNSG_M191128_A_MD\M191128_A_MD.mat')
%load('R:\CSN\md406\runNSG\runNSG_M191119_A_MD\M191119_A_MD.mat')
spikeMat = full(dataStruct.seriesData.M191128_A_MD_s201912122054361.shankData.shank1.spk);
%spikeMat = full(dataStruct.seriesData.M191119_A_MD_s2019120221020556.shankData.shank1.spk);
downsampledMat = downsampleRasterMatrix(spikeMat, 400, 0.2);
ops.iPC = 1:35;
ops.nC = 2;
[isort1, isort2, Sm] = mapTmap(downsampledMat,ops);

% Mean subtracted spiking
%figure; imagesc(zscore(downsampledMat-mean(mean(downsampledMat))));
figure;
h = pcolor(logical(flipud(downsampledMat)));
h.EdgeColor = 'none';
colormap(flipud(gray));
ax1 = axesProperties('Unit spiking in S1', 1, 'normal', 'off', 'w', 'Calibri', 20, 1, 2, [0.01 0.025], 'out', 'off', 'k',...
  'Time (s)', xlim, [1 400 800 1200], 'off', 'k', 'Unsorted units', ylim, fliplr([64 44 24 4]));
ax1.XTickLabel = {'0','400','800','1200'};
ax1.YTickLabel = {'60','40','20','0'};
set(gcf,'color','white');
filename = 'Unit_spiking_in_S1';
hgsave(gcf, filename);
print(gcf, [filename '.png'],'-dpng','-r300');

% Fully sorted
figure; imagesc(Sm);

% Time unsorted
newInd1(isort2) = 1:size(Sm,2);
Sm_timeUnsorted = Sm(:,newInd1);
figure; im = imagesc(Sm_timeUnsorted); hold on

% Find the dividing line
zscoreSum = sum(abs(Sm),2);
[~, dividingMember] = min(zscoreSum);
positiveSm = Sm(:,logical(Sm(dividingMember>0,:)));
if sum(sum(positiveSm(1:dividingMember-1,:)))/(dividingMember-1) >...
    sum(sum(positiveSm(dividingMember+1:end,:)))/(size(positiveSm,1)-dividingMember-1)
  firstGroup = 1:dividingMember;
  secondGroup = dividingMember+1:size(zscoreSum,1);
  dividingLine = dividingMember+0.5;
else
  firstGroup = 1:dividingMember-1;
  secondGroup = dividingMember:size(zscoreSum,1);
  dividingLine = dividingMember-0.5;
end
xLim = xlim;
%plot([xLim(1) xLim(2)], [dividingLine dividingLine], ':r'); hold off
%plot([xLim(1) xLim(2)], [dividingMember dividingMember], ':r'); hold off
%plot([xLim(1) xLim(2)], [46.5 46.5], ':r', 'Linewidth',2); hold off
ax1 = axesProperties('Manifold embedding algorithm', 1, 'normal', 'off', 'w', 'Calibri', 20, 1, 2, [0.01 0.025], 'out', 'off', 'k',...
  'Time (s)', xlim, [1 400 800 1200], 'off', 'k', 'Sorted units', ylim, [0.5 20 40 60]);
ax1.XTickLabel = {'0','400','800','1200'};
ax1.YTickLabel = {'0','20','40','60'};
set(gcf,'color','white');
colorbar;
filename = 'Manifold_embedding_of_S1';
hgsave(gcf, filename);
print(gcf, [filename '.png'],'-dpng','-r300');

figure; plot(zscoreSum);

% Fully unsorted
newInd2(isort1) = 1:size(Sm,1);
Sm_fullyUnsorted = Sm_timeUnsorted(newInd2,:);
figure; imagesc(Sm_fullyUnsorted);

clear newInd1 newInd2

% Sorted spiking
sortedRaster = downsampledMat(isort1,:);
imagesc(sortedRaster);

% z-scores sorted raster
zSortedRaster = sortedRaster-max(max(sortedRaster))/2;
imagesc(zSortedRaster);
colorbar;


[pcaCoef, explained, PCs, nPCs, prob] = pcaGeneric(downsampledMat);
PC1 = PCs(1,:);
[rPearson, pvalPearson] = corrMulti(PC1, downsampledMat, 'Pearson');
[rSpearman, pvalSpearman] = corrMulti(PC1, downsampledMat, 'Spearman');

[sortedR, iSort3] = sort(rSpearman, 'descend');
iSort3 = iSort3';

sortedRaster = downsampledMat(iSort3,:);

dividingMember = numel(sortedR) - find(sortedR < 0, 1, 'first') + 1;

% Sorted spiking
%figure; imagesc(zscore(downsampledMat-mean(mean(downsampledMat))));
figure;
h = pcolor(logical(flipud(downsampledMat(iSort3,:))));
h.EdgeColor = 'none';
colormap(flipud(gray));
ax1 = axesProperties('Unit spiking in S1', 1, 'normal', 'off', 'w', 'Calibri', 20, 1, 2, [0.01 0.025], 'out', 'off', 'k',...
  'Time (s)', xlim, [1 400 800 1200], 'off', 'k', 'Sorted units', ylim, fliplr([64 44 24 4]));
ax1.XTickLabel = {'0','400','800','1200'};
ax1.YTickLabel = {'60','40','20','0'};
set(gcf,'color','white');
filename = 'Unit_spiking_in_S1_after_sorting';
hgsave(gcf, filename);
print(gcf, [filename '.png'],'-dpng','-r300');

% Sorted spiking
%figure; imagesc(zscore(downsampledMat-mean(mean(downsampledMat))));
figure;
h = pcolor(logical(flipud(downsampledMat(iSort3,:))));
h.EdgeColor = 'none';
colormap(flipud(gray));
hold on; plot([xLim(1) xLim(2)], [dividingMember dividingMember], ':r', 'Linewidth',2); hold off
ax1 = axesProperties('Unit spiking in S1', 1, 'normal', 'off', 'w', 'Calibri', 20, 1, 2, [0.01 0.025], 'out', 'off', 'k',...
  'Time (s)', xlim, [1 400 800 1200], 'off', 'k', 'Sorted units', ylim, fliplr([64 44 24 4]));
ax1.XTickLabel = {'0','400','800','1200'};
ax1.YTickLabel = {'60','40','20','0'};
set(gcf,'color','white');
filename = 'Unit_spiking_in_S1_after_splitting';
hgsave(gcf, filename);
print(gcf, [filename '.png'],'-dpng','-r300');


%% PCA
% [pcaCoef, explained, PCs, nPCs, prob] = pcaGeneric(downsampledMat);
% PC1 = PCs(1,:);
% [rPearson, pvalPearson] = corrMulti(PC1, downsampledMat, 'Pearson');
% [rSpearman, pvalSpearman] = corrMulti(PC1, downsampledMat, 'Spearman');
% 
% [sortedR, iSort3] = sort(rSpearman, 'descend');
% iSort3 = iSort3';
% 
% sortedRaster = downsampledMat(iSort3,:);
% %figure; imagesc(zscore(sortedRaster - mean(mean(sortedRaster))));
% %figure; pcolor(zscore(sortedRaster - mean(mean(sortedRaster))));
% rescalledRaster = sortedRaster-1;
% rescalledRaster(rescalledRaster >= 0) = 0;
% rescalledRaster = rescalledRaster.*(-1);
% 
% h = pcolor(rescalledRaster);
% h.EdgeColor = 'none';
% colormap(gray)
% 
% h = pcolor(Sm_timeUnsorted);
% h.EdgeColor = 'none';
% colormap(gray)
% 
% % sortedSm = Sm_fullyUnsorted(iSort3,:);
% % figure; imagesc(zscore(sortedSm - mean(mean(sortedSm))));