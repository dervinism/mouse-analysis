function phaseFreqProfilePlotUpdate(figs, plotLegends, txtLegends, areas, conditions, FOI, areaPhaseFOIindividual, areaRecCount, areaRecSignificant, areaRecSignificantFOI, options)

% Initialise mean figures
for iArea = 1:numel(areas)
  figCondMean{iArea} = figure; %#ok<*AGROW>
  semilogx([10e-4 10e2],[0 0], 'k:');
  plotCondMean{iArea} = [];
  txtCondMean{iArea} = {};
end

for iCond = 1:numel(conditions) % Loop through conditions
  for iArea = 1:numel(areas) % Loop through areas
    if ~isempty(areaPhaseFOIindividual{iCond}{iArea})
      if options.fullRun
        figure(figs{iCond}{iArea}); hold on
      end
      
      % Get phase values
      phaseMean = NaN(1,numel(areaRecSignificantFOI{iCond}{iArea}));
      phaseStd = NaN(1,numel(areaRecSignificantFOI{iCond}{iArea}));
      phaseCI95 = NaN(1,numel(areaRecSignificantFOI{iCond}{iArea}));
      phaseCI95Fisher = NaN(1,numel(areaRecSignificantFOI{iCond}{iArea}));
      for i = 1:numel(areaRecSignificantFOI{iCond}{iArea})
        if sum(~isnan(areaPhaseFOIindividual{iCond}{iArea}(:,i)))
          iPhaseFOI = areaPhaseFOIindividual{iCond}{iArea}(:,i);
          phaseMean(i) = circmean(iPhaseFOI(~isnan(iPhaseFOI)));
          phaseStd(i) = circ_std(iPhaseFOI(~isnan(iPhaseFOI)), [], [], 1);
          phaseCI95(i) = circ_confmean(iPhaseFOI(~isnan(iPhaseFOI)), 0.05);
          phaseCI95Fisher(i) = circ_confmeanFisher(iPhaseFOI(~isnan(iPhaseFOI)), 0.05);
        end
      end
      
      % Update figures of individual phase frequency profiles with a mean phase frequency profile curve
      if options.fullRun
        p = semilogx(FOI(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)),...
          phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)),...
          'r:o', 'LineWidth',2, 'MarkerSize',5, 'MarkerEdgeColor','r', 'MarkerFaceColor','r'); % Mean
        uistack(p,'bottom');
        plotLegends{iCond}{iArea} = [plotLegends{iCond}{iArea} p];
        txtLegends{iCond}{iArea}{numel(txtLegends{iCond}{iArea})+1} = 'Mean FOI';
        
        pC1 = semilogx(FOI(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)),...
          phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)) + phaseCI95Fisher(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)),...
          'r:o', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','r', 'MarkerFaceColor','r'); % Upper 95% confidence limitof the mean
        uistack(pC1,'bottom');
        
        pC2 = semilogx(FOI(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)),...
          phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)) - phaseCI95Fisher(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)),...
          'r:o', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','r', 'MarkerFaceColor','r'); % Lower 95% confidence limitof the mean
        uistack(pC2,'bottom');
        plotLegends{iCond}{iArea} = [plotLegends{iCond}{iArea} pC2];
        txtLegends{iCond}{iArea}{numel(txtLegends{iCond}{iArea})+1} = '95% conf';
        hold off
      end
      
      % Display mean phase profiles on a separate figure as well
      figure(figCondMean{iArea}); hold on
      if iCond == 1
        p = semilogx(FOI(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)),...
          phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)),...
          'g:o', 'LineWidth',2, 'MarkerSize',5, 'MarkerEdgeColor','g', 'MarkerFaceColor','g'); % Mean during wakefulness
        txtCondMean{iArea}{numel(txtCondMean{iArea})+1} = 'Wakefulness';
      elseif iCond == 2
        p = semilogx(FOI(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)),...
          phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)),...
          'r:o', 'LineWidth',2, 'MarkerSize',5, 'MarkerEdgeColor','r', 'MarkerFaceColor','r'); % Mean during anaesthesia
        txtCondMean{iArea}{numel(txtCondMean{iArea})+1} = 'Anaesthesia';
        uistack(p,'bottom');
      end
      plotCondMean{iArea} = [plotCondMean{iArea} p];
      
      if iCond == 1
        pC = ciplot(phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)) - phaseCI95Fisher(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)),...
          phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)) + phaseCI95Fisher(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)),...
          FOI(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)), 'g', 0.1); % 95% confidence interval of the mean during wakefulness
      elseif iCond == 2
        pC = ciplot(phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)) - phaseCI95Fisher(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)),...
          phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)) + phaseCI95Fisher(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)),...
          FOI(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)), 'r', 0.05); % 95% confidence interval of the mean during anaesthesia
      end
      uistack(pC,'bottom');
      hold off
    end
    
    if ~exist([options.mainFolder filesep options.phaseFrequencyProfilesSubfolder], 'file')
      mkdir([options.mainFolder filesep options.phaseFrequencyProfilesSubfolder]);
    end
    
    % Tidy up the figures with individual profiles
    if options.fullRun
      figure(figs{iCond}{iArea});
      xlabel('Frequency (Hz)');
      xlim([10e-3 - 0.002   10e2])
      ylabel('Phase (rad)');
      %ylim([-(6/5)*pi (6/5)*pi])
      legend(plotLegends{iCond}{iArea},txtLegends{iCond}{iArea}, 'Interpreter','None');
      title([areas{iArea} ': ' conditions{iCond}...
        '   Significant recordings: ' num2str(areaRecSignificant{iCond}{iArea}) '/' num2str(areaRecCount{iCond}{iArea})])
      
      ylimits = ylim;
      xlimits = xlim;
      set(gcf,'color','w');
      legend boxoff
      ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 12, 4/3, 2, [0.005 0], 'out',...
        'on', 'k', {'Frequency (Hz)'}, xlimits, xticks,...
        'on', 'k', {'Phase (rad)'}, ylimits, [-pi 0 pi]);
      ax1.YTickLabel = {'-\pi','0','\pi'};
      
      set(gcf, 'Name',[areas{iArea} ': ' conditions{iCond}]);
      figFileName = [options.mainFolder filesep options.phaseFrequencyProfilesSubfolder filesep...
        areas{iArea} options.compString '_' conditions{iCond} '_phase'];
      hgsave(gcf, figFileName);
      
      label = [1.35 1.55];
      margin = [0.3 0.1];
      width = 1*options.figSize-label(1)-margin(1);
      height = 1*options.figSize-label(2)-margin(2);
      paperSize = resizeFig(gcf, ax1, width, height, label, margin, 0);
      exportFig(gcf, [figFileName '.png'],'-dpng','-r300', paperSize);
    end
    
    % Tidy up the figures with mean profiles
    if numel(conditions) == 1 || iCond == 2
      figure(figCondMean{iArea}); hold on
      xlabel('Frequency (Hz)');
      xlim([10e-3 - 0.002   10e2])
      ylabel('Phase (rad)');
      %ylim([-(6/5)*pi (6/5)*pi])
      legend(plotCondMean{iArea},txtCondMean{iArea}, 'Interpreter','None');
      title(['Means only: ' areas{iArea} options.compString])
      
      ylimits = [-pi pi]; %ylim;
      xlimits = xlim;
      set(gcf,'color','w');
      legend boxoff
      ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 12, 4/3, 2, [0.005 0], 'out',...
        'on', 'k', {'Frequency (Hz)'}, xlimits, xticks,...
        'on', 'k', {'Phase (rad)'}, ylimits, [-pi 0 pi]);
      ax1.YTickLabel = {'-\pi','0','\pi'};
      
      set(gcf, 'Name',['Means only: ' areas{iArea}]);
      figFileName = [options.mainFolder filesep options.phaseFrequencyProfilesSubfolder filesep...
        'Means_only__' areas{iArea} options.compString];
      hgsave(gcf, figFileName);
      
      label = [1.35 1.55];
      margin = [0.3 0.1];
      width = 1*options.figSize-label(1)-margin(1);
      height = 1*options.figSize-label(2)-margin(2);
      paperSize = resizeFig(gcf, ax1, width, height, label, margin, 0);
      exportFig(gcf, [figFileName '.png'],'-dpng','-r300', paperSize);
    end
  end
end