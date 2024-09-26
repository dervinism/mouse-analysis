function phaseFreqProfilePlotUpdate_units(areas, conditions, FOI, areaPhaseFOIindividual, areaPhaseFOIindividualPR, options)
  
  % Initialise figures
  for iArea = 1:numel(areas)
    figCondMean{iArea} = figure; %#ok<*AGROW>
    semilogx([10e-4 10e2],[0 0], 'k:');
    plotCondMean{iArea} = [];
    txtCondMean{iArea} = {};
  end
  
  for iCond = 1:numel(conditions) % Loop through conditions
    for iArea = 1:numel(areas) % Loop through areas
      if ~isempty(areaPhaseFOIindividual{iCond}{iArea})
        
        % Get phase values
        phaseMean = NaN(1,numel(FOI));
        phaseStd = NaN(1,numel(FOI));
        phaseCI95 = NaN(1,numel(FOI));
        phaseCI95Fisher = NaN(1,numel(FOI));
        phaseMean_PR = NaN(1,numel(FOI));
        phaseCI95Fisher_PR = NaN(1,numel(FOI));
        for i = 1:numel(FOI)
          if sum(~isnan(areaPhaseFOIindividual{iCond}{iArea}(:,i)))
            iPhaseFOI = areaPhaseFOIindividual{iCond}{iArea}(:,i);
            phaseMean(i) = circmean(iPhaseFOI(~isnan(iPhaseFOI)));
            phaseStd(i) = circ_std(iPhaseFOI(~isnan(iPhaseFOI)), [], [], 1);
            phaseCI95(i) = circ_confmean(iPhaseFOI(~isnan(iPhaseFOI)), 0.05);
            phaseCI95Fisher(i) = circ_confmeanFisher(iPhaseFOI(~isnan(iPhaseFOI)), 0.05);
            if ~isempty(areaPhaseFOIindividualPR) && ~isempty(areaPhaseFOIindividualPR{iCond}{iArea})
              iPhaseFOI_PR = areaPhaseFOIindividualPR{iCond}{iArea}(:,i);
              phaseMean_PR(i) = circmean(iPhaseFOI_PR(~isnan(iPhaseFOI_PR)));
              phaseCI95Fisher_PR(i) = circ_confmeanFisher(iPhaseFOI_PR(~isnan(iPhaseFOI_PR)), 0.05);
            end
          end
        end
        
        % Plot mean curves
        figure(figCondMean{iArea}); hold on
        if iCond == 1
          %           p1 = semilogx(FOI(~isnan(phaseMean) & ~isnan(phaseCI95)),...
          %             phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95)),...
          %             'g:o', 'LineWidth',2, 'MarkerSize',5,...
          %             'MarkerEdgeColor','g', 'MarkerFaceColor','g'); % Unit mean (parametric)
          p1 = semilogx(FOI(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)),...
            phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)),...
            'g:o', 'LineWidth',2, 'MarkerSize',5, 'MarkerEdgeColor','g', 'MarkerFaceColor','g'); % Unit mean
          txtCondMean{iArea}{numel(txtCondMean{iArea})+1} = 'Wakefulness';
          if ~isempty(areaPhaseFOIindividualPR)
            p2 = semilogx(FOI(~isnan(phaseMean_PR) & ~isnan(phaseCI95Fisher_PR)),...
              phaseMean_PR(~isnan(phaseMean_PR) & ~isnan(phaseCI95Fisher_PR)),...
              'g:', 'LineWidth',1); % Population mean
            txtCondMean{iArea}{numel(txtCondMean{iArea})+1} = 'Wakefulness PR';
          end
        elseif iCond == 2
          %           p1 = semilogx(FOI(~isnan(phaseMean) & ~isnan(phaseCI95)),...
          %             phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95)),...
          %             'r:o', 'LineWidth',2, 'MarkerSize',5,...
          %             'MarkerEdgeColor','r', 'MarkerFaceColor','r'); % Unit mean (parametric)
          p1 = semilogx(FOI(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)),...
            phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)),...
            'r:o', 'LineWidth',2, 'MarkerSize',5, 'MarkerEdgeColor','r', 'MarkerFaceColor','r'); % Unit mean
          txtCondMean{iArea}{numel(txtCondMean{iArea})+1} = 'Anaesthesia';
          uistack(p1,'bottom');
          if ~isempty(areaPhaseFOIindividualPR)
            p2 = semilogx(FOI(~isnan(phaseMean_PR) & ~isnan(phaseCI95Fisher_PR)),...
              phaseMean_PR(~isnan(phaseMean_PR) & ~isnan(phaseCI95Fisher_PR)),...
              'r:', 'LineWidth',1); % Population mean
            txtCondMean{iArea}{numel(txtCondMean{iArea})+1} = 'Anaesthesia PR';
          end
        end
        if ~isempty(areaPhaseFOIindividualPR)
          plotCondMean{iArea} = [plotCondMean{iArea} p1 p2];
          uistack(p2,'bottom');
        else
          plotCondMean{iArea} = [plotCondMean{iArea} p1];
        end
        
        % Plot 95% confidence intervals
        if iCond == 1
          %           pC = ciplot(phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95)) - phaseCI95(~isnan(phaseMean) & ~isnan(phaseCI95)),...
          %             phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95)) + phaseCI95(~isnan(phaseMean) & ~isnan(phaseCI95)),...
          %             FOI(~isnan(phaseMean) & ~isnan(phaseCI95)), 'g', 0.1); % Parametric
          pC = ciplot(phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)) - phaseCI95Fisher(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)),...
            phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)) + phaseCI95Fisher(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)),...
            FOI(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)), 'g', 0.1); % Non-parametric
        elseif iCond == 2
          %           pC = ciplot(phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95)) - phaseCI95(~isnan(phaseMean) & ~isnan(phaseCI95)),...
          %             phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95)) + phaseCI95(~isnan(phaseMean) & ~isnan(phaseCI95)),...
          %             FOI(~isnan(phaseMean) & ~isnan(phaseCI95)), 'r', 0.05); % Parametric
          pC = ciplot(phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)) - phaseCI95Fisher(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)),...
            phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)) + phaseCI95Fisher(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)),...
            FOI(~isnan(phaseMean) & ~isnan(phaseCI95Fisher)), 'r', 0.05); % Non-parametric
        end
        uistack(pC,'bottom');
        hold off
      end
      
      % Tidy and save the figure
      if numel(conditions) == 1 || iCond == 2
        figure(figCondMean{iArea}); hold on
        xlabel('Frequency (Hz)');
        xlim([10e-3 - 0.002   10e2])
        ylabel('Phase (rad)');
        %ylim([-(6/5)*pi (6/5)*pi])
        legend(plotCondMean{iArea},txtCondMean{iArea}, 'Interpreter','None');
        title(['Means only: ' areas{iArea}])
        
        ylimits = [-pi pi]; %ylim;
        xlimits = xlim;
        set(gcf,'color','w');
        legend boxoff
        ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 12, 4/3, 2, [0.005 0], 'out',...
          'on', 'k', {'Frequency (Hz)'}, xlimits, xticks,...
          'on', 'k', {'Phase (rad)'}, ylimits, [-pi 0 pi]);
        ax1.YTickLabel = {'-\pi','0','\pi'};
        
        set(gcf, 'Name',['Means only: ' areas{iArea}]);
        figFileName = [options.mainFolder filesep options.phaseFrequencyProfilesSubfolder filesep 'Means_only__'...
          areas{iArea} options.compString];
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
end