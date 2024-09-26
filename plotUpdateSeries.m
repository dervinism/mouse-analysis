function plotUpdateSeries(fH, fLegend, data, conditions, yLabel, titleStr, legendPos)
% Series figure updating. A helper function of globalPCA2.

for iCond = 1:numel(conditions)
  gLegend.lines = [];
  gLegend.text = {};
  figure(fH{iCond}); hold on
  if ~isempty(data{iCond})
    % Mean
    dataMean = sum(data{iCond},1)./size(data{iCond},1);
    p = plot(dataMean, 'r:o', 'LineWidth',2, 'MarkerSize',2, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    uistack(p,'bottom');
    fLegend{iCond}.lines = [fLegend{iCond}.lines p];
    fLegend{iCond}.text{numel(fLegend{iCond}.text)+1} = 'Mean';
    
    % 95% confidence interval
    sd = std(data{iCond},1,1, 'omitnan');
    SEM = sd ./ sqrt(size(data{iCond},1));
    CI95 = zeros(2,numel(dataMean));
    meanCI95 = zeros(2,numel(dataMean));
    for i = 1:numel(dataMean)
      CI95(:,i) = (tinv([0.025 0.975], size(data{iCond},1)-1))';
      meanCI95(:,i) = bsxfun(@times, SEM(i), CI95(:,i));
    end
    pC1 = plot(dataMean+meanCI95(1,:), 'r:o', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    uistack(pC1,'bottom');
    pC2 = plot(dataMean+meanCI95(2,:), 'r:o', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    uistack(pC2,'bottom');
    fLegend{iCond}.lines = [fLegend{iCond}.lines pC2];
    fLegend{iCond}.text{numel(fLegend{iCond}.text)+1} = '95% conf';
  end
  hold off
  
  % Tidy and save the figure
  figName = [conditions{iCond} ': ' titleStr];
  figFileName = [conditions{iCond} '_' strrep(titleStr, ' ', '_')];
  figFileName = strrep(figFileName, ':', '_');
  tidySaveFig(fH{iCond}, yLabel, fLegend{iCond}, figName, figFileName, legendPos)
end
end