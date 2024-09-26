function cohFreqProfilePlotUpdate(figsC, plotLegendsC, txtLegendsC, areas, conditions, FOI, areaCohFOI, areaCohFOIindividual, areaRecCount, areaRecSignificant, areaRecSignificantFOI, options)

if ~exist([options.mainFolder filesep options.coherenceFrequencyProfilesSubfolder], 'file')
  mkdir([options.mainFolder filesep options.coherenceFrequencyProfilesSubfolder]);
end

for iCond = 1:numel(conditions)
  for iArea = 1:numel(areas)
    figure(figsC{iCond}{iArea}); hold on
    if ~isempty(areaCohFOIindividual{iCond}{iArea})
      areaCohFOI{iCond}{iArea} = areaCohFOI{iCond}{iArea} ./ areaRecSignificantFOI{iCond}{iArea};
      p = semilogx(FOI,areaCohFOI{iCond}{iArea}, 'r:o', 'LineWidth',2, 'MarkerSize',5, 'MarkerEdgeColor','r', 'MarkerFaceColor','r'); % Mean
      uistack(p,'bottom');
      for iF = 1:numel(FOI)
        %text(FOI(iF),0.9, [num2str(areaRecSignificantFOI{iCond}{iArea}(iF)) '/' num2str(areaRecCount{iCond}{iArea})])
      end
      plotLegendsC{iCond}{iArea} = [plotLegendsC{iCond}{iArea} p];
      txtLegendsC{iCond}{iArea}{numel(txtLegendsC{iCond}{iArea})+1} = 'Mean FOI';
      
      cohStd = std(areaCohFOIindividual{iCond}{iArea}, 'omitnan');
      cohSEM = cohStd ./ sqrt(areaRecSignificantFOI{iCond}{iArea});
      CI95 = zeros(2,numel(areaRecSignificantFOI{iCond}{iArea}));
      cohCI95 = zeros(2,numel(areaRecSignificantFOI{iCond}{iArea}));
      for i = 1:numel(areaRecSignificantFOI{iCond}{iArea})
        CI95(:,i) = (tinv([0.025 0.975], areaRecSignificantFOI{iCond}{iArea}(i)-1))';
        cohCI95(:,i) = bsxfun(@times, cohSEM(i), CI95(:,i));
      end
      % Upper 95% confidence limit of the mean
      pC1 = semilogx(FOI,areaCohFOI{iCond}{iArea}+cohCI95(1,:), 'r:o', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
      uistack(pC1,'bottom');
      % Lower 95% confidence limit of the mean
      pC2 = semilogx(FOI,areaCohFOI{iCond}{iArea}+cohCI95(2,:), 'r:o', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
      uistack(pC2,'bottom');
      plotLegendsC{iCond}{iArea} = [plotLegendsC{iCond}{iArea} pC2];
      txtLegendsC{iCond}{iArea}{numel(txtLegendsC{iCond}{iArea})+1} = '95% conf';
    end
    hold off
    
    % Tidy up the figures
    xlabel('Frequency (Hz)');
    xlim([10e-3 - 0.002   10e2])
    ylabel('Coherence');
    legend(plotLegendsC{iCond}{iArea},txtLegendsC{iCond}{iArea}, 'Interpreter','None');
    title([areas{iArea} ': ' conditions{iCond}...
      '   Significant recordings: ' num2str(areaRecSignificant{iCond}{iArea}) '/' num2str(areaRecCount{iCond}{iArea})])
    
    ylimits = ylim;
    ylimits = [0 min([1 ylimits(2)])];
    xlimits = xlim;
    set(gcf,'color','w');
    legend boxoff
    ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 12, 4/3, 2, [0.005 0], 'out',...
      'on', 'k', {'Frequency (Hz)'}, xlimits, xticks,...
      'on', 'k', {'Coherence'}, ylimits, 0:0.1:1);
    
    set(gcf, 'Name',[areas{iArea} ': ' conditions{iCond}]);
    figFileName = [options.mainFolder filesep options.coherenceFrequencyProfilesSubfolder filesep...
      areas{iArea} options.compString '_' conditions{iCond} '_coherence'];
    hgsave(gcf, figFileName);
    
    label = [1.55 1.55];
    margin = [0.3 0.1];
    width = 1*options.figSize-label(1)-margin(1);
    height = 1*options.figSize-label(2)-margin(2);
    paperSize = resizeFig(gcf, ax1, width, height, label, margin, 0);
    exportFig(gcf, [figFileName '.png'],'-dpng','-r300', paperSize);
  end
end