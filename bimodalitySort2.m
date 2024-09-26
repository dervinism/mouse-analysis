load('R:\CSN\md406\runNSG\runNSG_M191119_A_MD\M191119_A_MD.mat')

fnsData = fieldnames(dataStruct.seriesData);
for dbCount = 1:numel(fnsData)
  spikeMat = full(dataStruct.seriesData.(fnsData{dbCount}).shankData.shank1.spk);
  downsampledMat = downsampleRasterMatrix(spikeMat, 400, 0.2);
  
  [pcaCoef, explained, PCs, nPCs, prob] = pcaGeneric(downsampledMat);
  PC1 = PCs(1,:);
  [rPearson, pvalPearson] = corrMulti(PC1, downsampledMat, 'Pearson');
  [rSpearman, pvalSpearman] = corrMulti(PC1, downsampledMat, 'Spearman');
  
  [sortedR, iSort3] = sort(rSpearman, 'descend');
  iSort3 = iSort3';
  
  sortedRaster = downsampledMat(iSort3,:);
  
  % Downsampled spiking
  figure;
  h = pcolor(logical(flipud(downsampledMat)));
  %h = pcolor(flipud(downsampledMat));
  h.EdgeColor = 'none';
  colormap(flipud(gray));
  ax1 = axesProperties('Unit spiking in S1', 1, 'normal', 'off', 'w', 'Calibri', 20, 1, 2, [0.01 0.025], 'out', 'off', 'k',...
    'Time (s)', xlim, [1 400 800 1200], 'off', 'k', 'Unsorted units', ylim, fliplr([64 44 24 4]));
  ax1.XTickLabel = {'0','400','800','1200'};
  ax1.YTickLabel = {'60','40','20','0'};
  set(gcf,'color','white');
  filename = 'Unit_spiking_in_S1';
%   hgsave(gcf, filename);
%   print(gcf, [filename '.png'],'-dpng','-r300');
  
  % Sorted downsampled spiking
  figure;
  h = pcolor(logical(flipud(downsampledMat(iSort3,:))));
  %h = pcolor(flipud(downsampledMat(iSort3,:)));
  h.EdgeColor = 'none';
  colormap(flipud(gray));
  ax1 = axesProperties('Unit spiking in S1', 1, 'normal', 'off', 'w', 'Calibri', 20, 1, 2, [0.01 0.025], 'out', 'off', 'k',...
    'Time (s)', xlim, [1 400 800 1200], 'off', 'k', 'Unsorted units', ylim, fliplr([64 44 24 4]));
  ax1.XTickLabel = {'0','400','800','1200'};
  ax1.YTickLabel = {'60','40','20','0'};
  set(gcf,'color','white');
  filename = 'Unit_spiking_in_S1';
%   hgsave(gcf, filename);
%   print(gcf, [filename '.png'],'-dpng','-r300');
  
  % Sorted spiking
  figure;
  h = pcolor(logical(flipud(spikeMat(iSort3,:))));
  %h = pcolor(flipud(downsampledMat(iSort3,:)));
  h.EdgeColor = 'none';
  colormap(flipud(gray));
  ax1 = axesProperties('Unit spiking in S1', 1, 'normal', 'off', 'w', 'Calibri', 20, 1, 2, [0.01 0.025], 'out', 'off', 'k',...
    'Time (s)', xlim, [1 400 800 1200], 'off', 'k', 'Unsorted units', ylim, fliplr([64 44 24 4]));
  ax1.XTickLabel = {'0','400','800','1200'};
  ax1.YTickLabel = {'60','40','20','0'};
  set(gcf,'color','white');
  filename = 'Unit_spiking_in_S1';
%   hgsave(gcf, filename);
%   print(gcf, [filename '.png'],'-dpng','-r300');
end