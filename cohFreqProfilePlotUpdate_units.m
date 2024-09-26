function cohFreqProfilePlotUpdate_units(areas, conditions, FOI, areaCohFOIindividual, areaCohFOIindividualPR, options)

% Initialise figures
  for iArea = 1:numel(areas)
    figCondMean{iArea} = figure; %#ok<*AGROW>
    semilogx([10e-4 10e2],[0 0], 'k:');
    plotCondMean{iArea} = [];
    txtCondMean{iArea} = {};
  end
  
  for iCond = 1:numel(conditions) % Loop through conditions
    for iArea = 1:numel(areas) % Loop through areas
      if ~isempty(areaCohFOIindividual{iCond}{iArea})
        
        % Get phase values
        cohMean = NaN(1,numel(FOI));
        cohStd = NaN(1,numel(FOI));
        cohSEM = NaN(2,numel(FOI));
        CI95 = NaN(2,numel(FOI));
        cohCI95 = NaN(2,numel(FOI));
        cohMean_PR = NaN(1,numel(FOI));
        cohStd_PR = NaN(1,numel(FOI));
        cohSEM_PR = NaN(2,numel(FOI));
        CI95_PR = NaN(2,numel(FOI));
        cohCI95_PR = NaN(2,numel(FOI));
        for i = 1:numel(FOI)
          if sum(~isnan(areaCohFOIindividual{iCond}{iArea}(:,i)))
            cohMean(i) = mean(areaCohFOIindividual{iCond}{iArea}(:,i), 'omitnan');
            cohStd(i) = std(areaCohFOIindividual{iCond}{iArea}(:,i), 'omitnan');
            cohSEM(:,i) = cohStd(i) ./ sqrt(sum(~isnan(areaCohFOIindividual{iCond}{iArea}(:,i))));
            CI95(:,i) = (tinv([0.025 0.975], sum(~isnan(areaCohFOIindividual{iCond}{iArea}(:,i)))-1))';
            cohCI95(:,i) = bsxfun(@times, cohSEM(i), CI95(:,i));
            if ~isempty(areaCohFOIindividualPR) && ~isempty(areaCohFOIindividualPR{iCond}{iArea})
              cohMean_PR(i) = mean(areaCohFOIindividualPR{iCond}{iArea}(:,i), 'omitnan');
              cohStd_PR(i) = std(areaCohFOIindividualPR{iCond}{iArea}(:,i), 'omitnan');
              cohSEM_PR(:,i) = cohStd_PR(i) ./ sqrt(sum(~isnan(areaCohFOIindividualPR{iCond}{iArea}(:,i))));
              CI95_PR(:,i) = (tinv([0.025 0.975], sum(~isnan(areaCohFOIindividualPR{iCond}{iArea}(:,i)))-1))';
              cohCI95_PR(:,i) = bsxfun(@times, cohSEM_PR(i), CI95_PR(:,i));
            end
          end
        end
        
        % Plot mean curves
        figure(figCondMean{iArea}); hold on
        if iCond == 1 && sum(~isnan(cohMean) & ~isnan(cohCI95(1,:)))
          %           p1 = semilogx(FOI(~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
          %             cohMean(~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
          %             'g:o', 'LineWidth',2, 'MarkerSize',5,...
          %             'MarkerEdgeColor','g', 'MarkerFaceColor','g'); % Unit mean (parametric)
          p1 = semilogx(FOI(~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
            cohMean(~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
            'g:o', 'LineWidth',2, 'MarkerSize',5, 'MarkerEdgeColor','g', 'MarkerFaceColor','g'); % Unit mean
          txtCondMean{iArea}{numel(txtCondMean{iArea})+1} = 'Wakefulness';
          if ~isempty(areaCohFOIindividualPR) && sum(~isnan(cohMean_PR) & ~isnan(cohCI95_PR(1,:)))
            p2 = semilogx(FOI(~isnan(cohMean_PR) & ~isnan(cohCI95_PR(1,:))),...
              cohMean_PR(~isnan(cohMean_PR) & ~isnan(cohCI95_PR(1,:))),...
              'g:', 'LineWidth',1); % Population mean
            txtCondMean{iArea}{numel(txtCondMean{iArea})+1} = 'Wakefulness PR';
          end
        elseif iCond == 2 && sum(~isnan(cohMean) & ~isnan(cohCI95(1,:)))
          %           p1 = semilogx(FOI(~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
          %             cohMean(~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
          %             'r:o', 'LineWidth',2, 'MarkerSize',5,...
          %             'MarkerEdgeColor','r', 'MarkerFaceColor','r'); % Unit mean (parametric)
          p1 = semilogx(FOI(~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
            cohMean(~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
            'r:o', 'LineWidth',2, 'MarkerSize',5, 'MarkerEdgeColor','r', 'MarkerFaceColor','r'); % Unit mean
          txtCondMean{iArea}{numel(txtCondMean{iArea})+1} = 'Anaesthesia';
          uistack(p1,'bottom');
          if ~isempty(areaCohFOIindividualPR) && sum(~isnan(cohMean_PR) & ~isnan(cohCI95_PR(1,:)))
            p2 = semilogx(FOI(~isnan(cohMean_PR) & ~isnan(cohCI95_PR(1,:))),...
              cohMean_PR(~isnan(cohMean_PR) & ~isnan(cohCI95_PR(1,:))),...
              'r:', 'LineWidth',1); % Population mean
            txtCondMean{iArea}{numel(txtCondMean{iArea})+1} = 'Anaesthesia PR';
          end
        end
        if ~isempty(areaCohFOIindividualPR) && sum(~isnan(cohMean) & ~isnan(cohCI95(1,:))) &&...
            sum(~isnan(cohMean_PR) & ~isnan(cohCI95_PR(1,:)))
          plotCondMean{iArea} = [plotCondMean{iArea} p1 p2];
          uistack(p2,'bottom');
        elseif sum(~isnan(cohMean) & ~isnan(cohCI95(1,:)))
          plotCondMean{iArea} = [plotCondMean{iArea} p1];
        end
        
        % Plot 95% confidence intervals
        if iCond == 1 && sum(~isnan(cohMean) & ~isnan(cohCI95(1,:)))
          %           pC = ciplot(cohMean(~isnan(cohMean) & ~isnan(cohCI95(1,:))) + cohCI95(1,~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
          %             cohMean(~isnan(cohMean) & ~isnan(cohCI95(2,:))) + cohCI95(2,~isnan(cohMean) & ~isnan(cohCI95(2,:))),...
          %             FOI(~isnan(cohMean) & ~isnan(cohCI95(1,:))), 'g', 0.1); % Parametric
          pC = ciplot(cohMean(~isnan(cohMean) & ~isnan(cohCI95(1,:))) + cohCI95(1,~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
            cohMean(~isnan(cohMean) & ~isnan(cohCI95(2,:))) + cohCI95(2,~isnan(cohMean) & ~isnan(cohCI95(2,:))),...
            FOI(~isnan(cohMean) & ~isnan(cohCI95(1,:))), 'g', 0.2); % Non-parametric
          uistack(pC,'bottom');
        elseif iCond == 2 && sum(~isnan(cohMean) & ~isnan(cohCI95(1,:)))
          %           pC = ciplot(cohMean(~isnan(cohMean) & ~isnan(cohCI95(1,:))) + cohCI95(1,~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
          %             cohMean(~isnan(cohMean) & ~isnan(cohCI95(2,:))) + cohCI95(2,~isnan(cohMean) & ~isnan(cohCI95(2,:))),...
          %             FOI(~isnan(cohMean) & ~isnan(cohCI95(1,:))), 'r', 0.05); % Parametric
          pC = ciplot(cohMean(~isnan(cohMean) & ~isnan(cohCI95(1,:))) + cohCI95(1,~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
            cohMean(~isnan(cohMean) & ~isnan(cohCI95(1,:))) + cohCI95(2,~isnan(cohMean) & ~isnan(cohCI95(2,:))),...
            FOI(~isnan(cohMean) & ~isnan(cohCI95(1,:))), 'r', 0.1); % Non-parametric
          uistack(pC,'bottom');
        end
        hold off
      end
      
      % Tidy and save the figure
      if (numel(conditions) == 1 && ~isempty(areaCohFOIindividual{1}{iArea})) ||...
          (iCond == 2 && (~isempty(areaCohFOIindividual{1}{iArea}) || ~isempty(areaCohFOIindividual{2}{iArea})))
        figure(figCondMean{iArea}); hold on
        xlabel('Frequency (Hz)');
        xlim([10e-3 - 0.002   10e2])
        ylabel('Coherence');
        if ~isempty(plotCondMean{iArea})
          assert(numel(plotCondMean{iArea}) == numel(txtCondMean{iArea}));
          legend(plotCondMean{iArea},txtCondMean{iArea}, 'Interpreter','None');
        else
          continue
        end
        title(['Means only: ' areas{iArea}])
        
        ylimits = ylim;
        ylim([0 min([1 ylimits(2)])]);
        xlimits = xlim;
        set(gcf,'color','w');
        legend boxoff
        ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 12, 4/3, 2, [0.005 0], 'out',...
          'on', 'k', {'Frequency (Hz)'}, xlimits, xticks,...
          'on', 'k', {'Coherence'}, ylim, 0:0.1:1);
        %ax1.YTickLabel = {'-\pi','0','\pi'};
        
        set(gcf, 'Name',['Means only: ' areas{iArea}]);
        figFileName = [options.mainFolder filesep options.coherenceFrequencyProfilesSubfolder filesep 'Means_only__'...
          areas{iArea} options.compString];
        
        label = [1.55 2.5];
        margin = [0.3 0.1];
        width = 1*options.figSize-label(1)-margin(1);
        height = 1*options.figSize-label(2)-margin(2);
        paperSize = resizeFig(gcf, ax1, width, height, label, margin, 0);
        hgsave(gcf, figFileName);
        exportFig(gcf, [figFileName '.png'],'-dpng','-r300', paperSize);
      end
    end
  end
end