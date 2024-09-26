function fH = meanPhaseProfilesSummary2(conditions, areas, areas_reverse, AOI, FOI, areaPhaseFOIindividual_reg, areaPhaseFOIindividual_reverse)
% fH = meanPhaseProfilesSummary2(conditions, areas, areas_reverse, AOI, FOI,
% areaPhaseFOIindividual_reg, areaPhaseFOIindividual_reverse, drawCI).
%
% Function produces figures with mean phase frequency profiles. The
% input variables are stored in globalUnits_ca.mat and
% globalUnits_ca_reverse.mat if you want to get figures for units. If you
% are interested in the population rate comparisons, load data from
% area2areaCohMats.mat and area2areaCohMats_reverse.mat. Function outputs
% figure handles.

% Get combination matrix
areaCode = zeros(1,numel(AOI));
for iArea = 1:numel(AOI)
  areaCode(iArea) = determineArea(AOI{iArea}); %#ok<*AGROW>
end
combos = perms(areaCode);


% Initialise figures
fH = {}; % figure handles
pH = {}; % plot handles
lH = {}; % legend text entries
for iCond = 1:numel(conditions)
  for iCombos = 1:size(combos,1)
    figCondMean{iCombos} = figure;
    semilogx([10e-4 10e2],[0 0], 'k:');
    plotCondMean{iCombos} = [];
    txtCondMean{iCombos} = {};
  end
  fH{iCond} = figCondMean;
  pH{iCond} = plotCondMean;
  lH{iCond} = txtCondMean;
end


for iCond = 1:numel(conditions)
  for iCombos = 1:size(combos,1)
    combo = combos(iCombos,:);
    phaseMean = NaN(numel(combo),numel(FOI));
    phaseCI95Fisher = NaN(numel(combo),numel(FOI));
    areaNameCombo = {};
    for iComp = 1:numel(combo)
      
      
      % Obtain mean phase and the 95% confidence interval for all comparisons
      if iComp == numel(combo)
        [~, ~, areaName1] = determineArea(combo(end));
        [~, ~, areaName2] = determineArea(combo(1));
      else
        [~, ~, areaName1] = determineArea(combo(iComp+1));
        [~, ~, areaName2] = determineArea(combo(iComp));
      end
      areaName = [areaName1 'Vs' areaName2];
      if strcmpi(areaName, 'RSCVsCA')
        areaName = 'RSCVsCA1';
      elseif strcmpi(areaName, 'CAVsRSC')
        areaName = 'CA1VsRSC';
      end
      areaInd = [];
      assert(numel(areas) == numel(areas_reverse));
      for i = 1:numel(areas)
        if strcmpi(areaName, areas{i})
          areaInd = i;
          areaPhaseFOIindividual = areaPhaseFOIindividual_reg;
          break
        elseif strcmpi(areaName, areas_reverse{i})
          areaInd = i;
          areaPhaseFOIindividual = areaPhaseFOIindividual_reverse;
          break
        elseif i == numel(areas)
          error('Area comparison does not exist.');
        end
      end
      if ~isempty(areaPhaseFOIindividual{iCond}{areaInd})
        for f = 1:numel(FOI)
          if sum(~isnan(areaPhaseFOIindividual{iCond}{areaInd}(:,f)))
            iPhaseFOI = areaPhaseFOIindividual{iCond}{areaInd}(:,f);
            phaseMean(iComp,f) = circmean(iPhaseFOI(~isnan(iPhaseFOI)));
            phaseCI95Fisher(iComp,f) = circ_confmeanFisher(iPhaseFOI(~isnan(iPhaseFOI)), 0.05);
          end
        end
      end
      areaNameCombo{iComp} = areaName;
    end
    
    
    % Plot the mean phase profile summary figure
    if iCond == 1
      colour = [0 255 0]./255; % green
    elseif iCond == 2
      colour = [255 0 0]./255; % red
    end
    lineType = {':o', '-.o', '--o', '-o'};
    figure(fH{iCond}{iCombos})
    hold on
    for iComp = 1:numel(combo)
      inds = ~isnan(phaseMean(iComp,:)) & ~isnan(phaseCI95Fisher(iComp,:));
      jitteredFOI = FOI(inds)+(iComp-2.5)*0.04.*FOI(inds);
      
      % Phase
      phaseMeanComp = phaseMean(iComp,inds);
      pPhase = semilogx(jitteredFOI, phaseMeanComp, lineType{iComp}, 'LineWidth',2, 'MarkerSize',5,...
        'Color',colour, 'MarkerEdgeColor',colour, 'MarkerFaceColor',colour);
      pH{iCond}{iCombos} = [pH{iCond}{iCombos} pPhase];
      lH{iCond}{iCombos}{numel(lH{iCond}{iCombos})+1} = ['\Phi: ' strrep(areaNameCombo{iComp}, 'Vs', ' v ')];
      
      % 95% confidence interval
      phaseCI95FisherComp = phaseCI95Fisher(iComp,inds);
      for iBarCI = 1:numel(jitteredFOI)
        plot([jitteredFOI(iBarCI) jitteredFOI(iBarCI)],...
          [phaseMeanComp(iBarCI) - phaseCI95FisherComp(iBarCI) phaseMeanComp(iBarCI) + phaseCI95FisherComp(iBarCI)],...
          '-', 'Color',colour, 'LineWidth',2);
      end
      colour = colour./1.5;
    end
    ylimits = ylim;
    if ylimits(1) > -pi
      ylimits(1) = -pi;
    end
    if ylimits(2) < pi
      ylimits(2) = pi;
    end
    ylim(ylimits);
    %plot([0.02 0.02], ylim, 'k:')
    plot([0.2 0.2], ylim, 'k:')
    plot([4 4], ylim, 'k:')
    hold off
    
    % Tidy the figure
    xlabel('Frequency (Hz)');
    xlim([10e-3 - 0.002   10e2])
    ylabel('Phase (rad)');
    legend(pH{iCond}{iCombos},lH{iCond}{iCombos}, 'Interpreter','tex');
    figTitle = ['Means only: ' conditions{iCond} ' '];
    for iComp = 1:numel(combo)
      figTitle = [figTitle areaNameCombo{iComp} '+'];
    end
    figTitle = figTitle(1:end-1);
    title(figTitle);
    
    ylimits = ylim;
    xlimits = xlim;
    set(gcf,'color','w');
    legend boxoff
    ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 12, 4/3, 2, [0.005 0], 'out',...
      'on', 'k', {'Frequency (Hz)'}, xlimits, xticks,...
      'on', 'k', {'Phase (rad)'}, ylimits, [-pi 0 pi]);
    ax1.YTickLabel = {'-\pi','0','\pi'};
    
    set(gcf, 'Name',figTitle);
    figFileName = strrep(figTitle, ' ', '_');
    figFileName = strrep(figFileName, ':', '_');
    figFileName = strrep(figFileName, '.', 'p');
    hgsave(gcf, figFileName);
    
    label = [1.35 1.55];
    margin = [0.3 0.1];
    width = 1*15-label(1)-margin(1);
    height = 1*15-label(2)-margin(2);
    paperSize = resizeFig(gcf, ax1, width, height, label, margin, 0);
    exportFig(gcf, [figFileName '.png'],'-dpng','-r300', paperSize);
    
    % Save portions of the figure separately
    xlim([0.01 0.2])
    ylim(ylimits);
    figFileNameIS = [figFileName '_infraslow.png'];
    hgsave(gcf, figFileNameIS);
    exportFig(gcf, figFileNameIS,'-dpng','-r300', paperSize);
    xlim([0.2 4])
    ylim(ylimits);
    figFileNameS = [figFileName '_slow.png'];
    hgsave(gcf, figFileNameS);
    exportFig(gcf, figFileNameS,'-dpng','-r300', paperSize);
    xlim([4 160])
    ylim(ylimits);
    figFileNameF = [figFileName '_fast.png'];
    hgsave(gcf, figFileNameF);
    exportFig(gcf, figFileNameF,'-dpng','-r300', paperSize);
  end
end