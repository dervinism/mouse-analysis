function fH = meanPhaseProfilesSummary3(conditions, areas, areas_reverse, AOI, freq, FOI, areaPhaseFOIindividual_reg, areaPhaseFOIindividual_reverse, outputFolder, opt)
% fH = meanPhaseProfilesSummary3(conditions, areas, areas_reverse, AOI, freq, FOI, areaPhaseFOIindividual_reg, areaPhaseFOIindividual_reverse, outputFolder, opt)
%
% Function produces figures with mean phase frequency profiles. The
% input variables are stored in globalUnits_ca.mat and
% globalUnits_ca_reverse.mat if you want to get figures for units. If you
% are interested in the population rate comparisons, load data from
% area2areaCohMats.mat and area2areaCohMats_reverse.mat. Function outputs
% figure handles.

if nargin < 10
  %opt.xlim = [0.04 30];
  opt.xlim = [0.01 30];
  opt.ylim = [-pi pi];
  opt.legendLoc = 'NorthEast';
  opt.subpop = 'all';
else
  if ~isfield(opt, 'xlim')
    %opt.xlim = [0.04 30];
    opt.xlim = [0.008 30];
  end
  if ~isfield(opt, 'ylim')
    opt.ylim = [-pi pi];
  end
  if ~isfield(opt, 'legendLoc')
    opt.legendLoc = 'NorthEast';
  end
  if ~isfield(opt, 'subpop')
    opt.subpop = 'all';
  end
end

% Get combination matrix
if iscell(AOI{1}) && numel(AOI{1}) > 1 && iscell(AOI{2}) && numel(AOI{2}) > 1
  nComp = numel(AOI{1});
  areaCode = zeros(nComp,2);
  for iComp = 1:nComp
    area1Code = determineArea(AOI{1}{iComp});
    area2Code = determineArea(AOI{2}{iComp});
    areaCode(iComp,:) = [area1Code(1) area2Code(1)];
  end
elseif iscell(AOI{1}) && numel(AOI{1}) > 1
  nComp = numel(AOI{1});
  areaCode = zeros(nComp,2);
  for iComp = 1:nComp
    area1Code = determineArea(AOI{1}{iComp});
    area2Code = determineArea(AOI{2});
    areaCode(iComp,:) = [area1Code(1) area2Code(1)];
  end
else
  nComp = numel(AOI{2});
  areaCode = zeros(nComp,2);
  for iComp = 1:nComp
    area1Code = determineArea(AOI{1});
    area2Code = determineArea(AOI{2}{iComp});
    areaCode(iComp,:) = [area1Code(1) area2Code(1)];
  end
end


% Initialise figures
fH = {}; % figure handles
pH = {}; % plot handles
lH = {}; % legend text entries
for iCond = 1:3
  fH{iCond} = figure; %#ok<*AGROW>
  semilogx([10e-4 10e2],[0 0], 'k:');
  pH{iCond} = [];
  lH{iCond} = {};
end


for iCond = 1:numel(conditions)
  if numel(conditions) == 1 && strcmp(conditions{iCond}, 'anaesthesia')
    iCond = 2; %#ok<*FXSET>
  end
  phaseMean = NaN(nComp,numel(freq));
  phaseCI95Fisher = NaN(nComp,numel(freq));
  areaNameComp = {};
  
  
  % Obtain mean phase and the 95% confidence interval for all comparisons
  for iComp = 1:nComp
    if iscell(AOI{1}) && numel(AOI{1}) > 1 && iscell(AOI{2}) && numel(AOI{2}) > 1
      areaName = [AOI{1}{iComp} 'Vs' AOI{2}{iComp}];
      areaName1 = [AOI{1}{iComp}(2:end) 'Vs' AOI{2}{iComp}];
      areaName2 = [AOI{1}{iComp} 'Vs' AOI{2}{iComp}(2:end)];
      areaName3 = [AOI{1}{iComp}(2:end) 'Vs' AOI{2}{iComp}(2:end)];
    elseif iscell(AOI{1}) && numel(AOI{1}) > 1
      areaName = [AOI{1}{iComp} 'Vs' AOI{2}];
      areaName1 = [AOI{1}{iComp}(2:end) 'Vs' AOI{2}];
      areaName2 = [AOI{1}{iComp} 'Vs' AOI{2}(2:end)];
      areaName3 = [AOI{1}{iComp}(2:end) 'Vs' AOI{2}(2:end)];
    else
      areaName = [AOI{1} 'Vs' AOI{2}{iComp}];
      areaName1 = [AOI{1}(2:end) 'Vs' AOI{2}{iComp}];
      areaName2 = [AOI{1} 'Vs' AOI{2}{iComp}(2:end)];
      areaName3 = [AOI{1}(2:end) 'Vs' AOI{2}{iComp}(2:end)];
    end
    areaInd = [];
    assert(numel(areas) == numel(areas_reverse));
    for i = 1:numel(areas)
      if strcmpi(areaName, areas{i}) || strcmpi(areaName1, areas{i}) || strcmpi(areaName2, areas{i}) || strcmpi(areaName3, areas{i})
        areaInd = i;
        areaPhaseFOIindividual = areaPhaseFOIindividual_reg;
        break
      elseif strcmpi(areaName, areas_reverse{i}) || strcmpi(areaName1, areas_reverse{i}) || strcmpi(areaName2, areas_reverse{i}) || strcmpi(areaName3, areas_reverse{i})
        areaInd = i;
        areaPhaseFOIindividual = areaPhaseFOIindividual_reverse;
        break
      elseif i == numel(areas)
        error('Area comparison does not exist.');
      end
    end
    if ~isempty(areaPhaseFOIindividual{iCond}{areaInd})
      for f = 1:numel(freq)
        if sum(~isnan(areaPhaseFOIindividual{iCond}{areaInd}(:,f)))
          iPhaseFOI = areaPhaseFOIindividual{iCond}{areaInd}(:,f);
          phaseMean(iComp,f) = circmean(iPhaseFOI(~isnan(iPhaseFOI)));
          phaseCI95Fisher(iComp,f) = circ_confmeanFisher(iPhaseFOI(~isnan(iPhaseFOI)), 0.05);
        end
      end
    end
    areaNameComp{iComp} = areaName;
  end
  
  
  % Plot the mean phase profile summary figure
  if iCond == 1
    colour = [0 255 0]./255; % green
  elseif iCond == 2
    colour = [255 0 0]./255; % red
  end
  %lineType = {':o', '--o', '-o'};
  lineType = {':', '-.', '--', '-', '-'};
  figure(fH{iCond})
  hold on
  
  for iComp = 1:nComp
    inds = ~isnan(phaseMean(iComp,:)) & ~isnan(phaseCI95Fisher(iComp,:));
    jitteredFreq = freq(inds)+(iComp-(nComp/2))*0.04.*freq(inds);
    jitteredFOI = FOI+(iComp-(nComp/2))*0.04.*FOI;
    %jitteredFreq = freq(inds);
    %jitteredFOI = FOI;
    
    indsF = ismember(jitteredFreq,jitteredFOI);
    
    % Phase
    phaseMeanComp = bestUnwrap(recentrePhase(bestUnwrap(phaseMean(iComp,inds)), 0));
    %phaseMeanComp = recentrePhase(bestUnwrap(phaseMean(iComp,inds)), pi-2.548);
    %phaseMeanComp = recentrePhase(bestUnwrap(phaseMean(iComp,inds)), pi-2.545);
    %cutoff = 1.5;
    %phaseMeanComp = [recentrePhase(bestUnwrap(phaseMean(iComp,inds(1:round(numel(inds)/cutoff)))),0) recentrePhase(bestUnwrap(phaseMean(iComp,[false(1,round(numel(inds)/cutoff)) inds(round(numel(inds)/cutoff)+1:end)])),-pi/2.7)];
    %cutoff = 1.8;
    %phaseMeanComp = [recentrePhase(bestUnwrap(phaseMean(iComp,inds(1:round(numel(inds)/cutoff)))),pi/4) recentrePhase(bestUnwrap(phaseMean(iComp,[false(1,round(numel(inds)/cutoff)) inds(round(numel(inds)/cutoff)+1:end)])),-pi/2.7)];
    %cutoff = 2;
    %phaseMeanComp = [recentrePhase(bestUnwrap(phaseMean(iComp,inds(1:round(numel(inds)/cutoff)))),0) recentrePhase(bestUnwrap(phaseMean(iComp,[false(1,round(numel(inds)/cutoff)) inds(round(numel(inds)/cutoff)+1:end)])),pi/2.7)];
    %phaseMeanComp = [recentrePhase(bestUnwrap(phaseMean(iComp,inds(1:round(numel(inds)/cutoff)))),0) recentrePhase(bestUnwrap(phaseMean(iComp,[false(1,round(numel(inds)/cutoff)) inds(round(numel(inds)/cutoff)+1:end)])),-pi/2.25)];
    %pPhase = semilogx(jitteredFreq, phaseMeanComp, lineType{iComp}, 'LineWidth',2, 'MarkerSize',5,...
    %  'Color',colour, 'MarkerEdgeColor',colour, 'MarkerFaceColor',colour);
    pPhase = semilogx(jitteredFreq, phaseMeanComp, lineType{iComp}, 'LineWidth',2, 'Color',colour);
    pH{iCond} = [pH{iCond} pPhase];
    areaNameComp{iComp} = strrep(areaNameComp{iComp}, 'l', '');
    areaNameComp{iComp} = strrep(areaNameComp{iComp}, 'r', '');
    if strcmpi(opt.subpop, 'all')
      lH{iCond}{numel(lH{iCond})+1} = ['\Phi: unit ' strrep(areaNameComp{iComp}, 'Vs', ' wrt ') ' PR'];
    elseif strcmpi(opt.subpop, 'positive')
      lH{iCond}{numel(lH{iCond})+1} = ['\Phi: unit ' strrep(areaNameComp{iComp}, 'Vs', '+ wrt ') ' PR+'];
    elseif strcmpi(opt.subpop, 'negative')
      lH{iCond}{numel(lH{iCond})+1} = ['\Phi: unit ' strrep(areaNameComp{iComp}, 'Vs', '- wrt ') ' PR-'];
    end
    
    % 95% confidence interval
    phaseCI95FisherComp = phaseCI95Fisher(iComp,inds);
    for iBarCI = 1:numel(jitteredFreq)
      if indsF(iBarCI)
        plot([jitteredFreq(iBarCI) jitteredFreq(iBarCI)],...
          [phaseMeanComp(iBarCI) - phaseCI95FisherComp(iBarCI) phaseMeanComp(iBarCI) + phaseCI95FisherComp(iBarCI)],...
          '-', 'Color',colour, 'LineWidth',1);
      end
    end
    colour = colour./1.5;
    ylimits = opt.ylim;
    if ylimits(1) > -pi
      ylimits(1) = -pi;
    end
    if ylimits(2) < pi
      ylimits(2) = pi;
    end
    ylim(ylimits);
    %plot([0.02 0.02], opt.ylim, 'k:')
    plot([0.2 0.2], opt.ylim, 'k:')
    plot([16 16], opt.ylim, 'k:')
  end
  hold off
  
  % Tidy the figure
  xlabel('Frequency (Hz)');
  xlim([10e-3 - 0.002   10e2])
  ylabel('Phase (rad)');
  legend(pH{iCond},lH{iCond}, 'Interpreter','tex', 'Location',opt.legendLoc);
  if numel(conditions) == 1 && strcmp(conditions{1}, 'anaesthesia')
    figTitle = ['Means only: ' conditions{1} ' '];
  else
    figTitle = ['Means only: ' conditions{iCond} ' '];
  end
  for iComp = 1:nComp
    figTitle = [figTitle areaNameComp{iComp} '+'];
  end
  figTitle = figTitle(1:end-1);
  title(figTitle);
  
  ylimits = opt.ylim;
  xlimits = opt.xlim; %xlim;
  set(gcf,'color','w');
  legend boxoff
  if abs(ylimits(2)-ylimits(1)) < pi/2
    ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 20, 4/3, 2, [0.005 0], 'out',...
      'on', 'k', {'Frequency (Hz)'}, xlimits, xticks,...
      'on', 'k', {'Phase (rad)'}, ylimits, [-pi/8 -pi/16 0 pi/16 pi/8]);
    ax1.YTickLabel = {'-\pi/8','-\pi/16','0','\pi/16','\pi/8'};
  elseif abs(ylimits(2)-ylimits(1)) < pi
    ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 20, 4/3, 2, [0.005 0], 'out',...
      'on', 'k', {'Frequency (Hz)'}, xlimits, xticks,...
      'on', 'k', {'Phase (rad)'}, ylimits, [-pi/4 -pi/8 0 pi/8 pi/4]);
    ax1.YTickLabel = {'-\pi/4','-\pi/8','0','\pi/8','\pi/4'};
  elseif abs(ylimits(2)-ylimits(1)) < 2*pi
    ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 20, 4/3, 2, [0.005 0], 'out',...
      'on', 'k', {'Frequency (Hz)'}, xlimits, xticks,...
      'on', 'k', {'Phase (rad)'}, ylimits, [-pi -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi]);
    ax1.YTickLabel = {'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'};
  else
    ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 20, 4/3, 2, [0.005 0], 'out',...
      'on', 'k', {'Frequency (Hz)'}, xlimits, xticks,...
      'on', 'k', {'Phase (rad)'}, ylimits, [-pi -pi/2 0 pi/2 pi]);
    ax1.YTickLabel = {'-\pi','-\pi/2','0','\pi/2','\pi'};
  end
  %ax1.XAxis.Exponent = 0;
  %xtickformat('%.2f')
  xticklabels({'0.01','0.1','1','10','100',})
  
  set(gcf, 'Name',figTitle);
  figFileName = strrep(figTitle, ' ', '_');
  figFileName = strrep(figFileName, ':', '_');
  figFileName = strrep(figFileName, '.', 'p');
%   hgsave(gcf, ['R:\CSN\Shared\Dynamics\Data\mouse_analysis_data\area2area_comparisons\units\meanPhaseSums'...
%     filesep figFileName]);
  if ~exist(outputFolder, 'file')
    mkdir(outputFolder);
  end
  hgsave(gcf, [outputFolder filesep figFileName]);
  
  if opt.ylim(1) == -pi/2 && opt.ylim(2) == pi/2
    label = [2.75 2.2];
    margin = [0.3 0.35];
  else
    label = [2.2 2.2];
    margin = [0.3 0.1];
  end
  width = 1*15-label(1)-margin(1);
  height = 1*15-label(2)-margin(2);
  paperSize = resizeFig(gcf, ax1, width, height, label, margin, 0);
  exportFig(gcf, [outputFolder filesep figFileName '.png'],'-dpng','-r300', paperSize);
  
  % Save portions of the figure separately
%   xlim([0.01 0.2])
%   ylim(ylimits);
%   figFileNameIS = [figFileName '_infraslow.png'];
%   hgsave(gcf, figFileNameIS);
%   exportFig(gcf, figFileNameIS,'-dpng','-r300', paperSize);
%   xlim([0.2 4])
%   ylim(ylimits);
%   figFileNameS = [figFileName '_slow.png'];
%   hgsave(gcf, figFileNameS);
%   exportFig(gcf, figFileNameS,'-dpng','-r300', paperSize);
%   xlim([4 160])
%   ylim(ylimits);
%   figFileNameF = [figFileName '_fast.png'];
%   hgsave(gcf, figFileNameF);
%   exportFig(gcf, figFileNameF,'-dpng','-r300', paperSize);
end