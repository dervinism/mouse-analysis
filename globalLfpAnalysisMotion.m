% Run this script to perform analyses comparing LFP ripple rate and total movement data.

cleanUp



%% INITIALISE PARAMETERS
params
lists



%% DISPLAY RIPPLE RATE FREQUENCY PHASE PROFILES FOR ALL ANIMALS AND ALL RECORDINGS
% LOOP THROUGH DB ENTRIES
recCount = 0;
recSignificant = 0;
for dbCount = 1:numel(fullSeries)
  strSep = strfind(fullSeries{dbCount},'s');
  series = fullSeries{dbCount}(strSep+1:end);
  animal = fullSeries{dbCount}(1:strSep-2);
  animalColour = animalColours(animal);
  if dbCount == 1 || ~strcmpi(animal,prevAnimal)
    try
      load([dataDir filesep animal filesep 'CAR' filesep animal '.mat'])
    catch
      load([dataDir filesep animal filesep animal '.mat'])
    end
    updateLegend = true;
  end
  prevAnimal = animal;
  
% LOAD THE CONTENTS OF THE DB STRUCTURE VARIABLE
  dbStruct = dataStruct.seriesData.(fullSeries{dbCount});
  if dbCount == 1
    FOI = dbStruct.FOI;
    FOI(FOI > 2) = [];
    totPhaseFOI = [];
  end
  if ~isfield(dbStruct,'lfpphaseCohDataMotion')
    disp(['No phase data for ' fullSeries{dbCount} '. Skipping to the next db entry...']);
    continue
  end
  recCount = recCount + 1;
  
% PLOT THE DATA
  phase = bestUnwrap(dbStruct.lfpphaseCohDataMotion.rippleRate{ch{dbCount}}.phase);
  phase(isnan(dbStruct.lfpphaseCohDataMotion.rippleRate{ch{dbCount}}.phase_confU)) = NaN;
  phase(isnan(dbStruct.lfpphaseCohDataMotion.rippleRate{ch{dbCount}}.phase_confL)) = NaN;
  phase(isnan(dbStruct.lfpphaseCohDataMotion.rippleRate{ch{dbCount}}.coh_conf)) = NaN;
  freq = dbStruct.lfpphaseCohDataMotion.rippleRate{ch{dbCount}}.freq;
  [~, endFreq] = min(abs(freq - 2));
  phase = phase(1:endFreq);
  freq = freq(1:endFreq);
%   if sum(isnan(phase)) > 0.6*numel(phase)
%     continue
%   end
  recSignificant = recSignificant + 1;
  if mean(phase, 'omitnan') > 2*pi
    phase = phase - 2*pi;
  end
  if mean(phase, 'omitnan') < -(0.15)*pi
    phase = phase + 2*pi;
  end
  if ~exist('figH','var')
    figH = figure('units', 'normalized', 'position', [0.002, .04, 1, .88], 'Visible', 'on');
    semilogx([10e-4 10e2],[0 0], 'k:'); hold on
    p = semilogx(freq,phase, 'Color',animalColour);
    pLegend = p;
    txtLegend = {animal};
    updateLegend = false;
  else
    if updateLegend
      p = semilogx(freq,phase, 'Color',animalColour);
      pLegend = [pLegend p];
      txtLegend{numel(txtLegend)+1} = animal;
      updateLegend = false;
    else
      semilogx(freq,phase, 'Color',animalColour);
    end
  end
  
  if isempty(totPhaseFOI)
  	totPhaseFOI = phaseCohFOI(FOI, freq, phase, ones(1,numel(phase)), ones(1,numel(phase)), ones(1,numel(phase)));
  else
    totPhaseFOI = [totPhaseFOI; phaseCohFOI(FOI, freq, phase, ones(1,numel(phase)), ones(1,numel(phase)), ones(1,numel(phase)))]; %#ok<*AGROW>
  end
end

% Mean and 95% confidence intervals
phaseMean = NaN(1,numel(FOI));
phaseStd = NaN(1,numel(FOI));
phaseCI95 = NaN(1,numel(FOI));
for i = 1:numel(FOI)
  iPhaseFOI = totPhaseFOI(:,i);
  if sum(~isnan(iPhaseFOI))
    phaseMean(i) = circmean(iPhaseFOI(~isnan(iPhaseFOI)));
    if phaseMean(i) < -0.5
      phaseMean(i) = phaseMean(i) + 2*pi;
    end
    phaseStd(i) = circ_std(iPhaseFOI(~isnan(iPhaseFOI)), [], [], 1);
    phaseCI95(i) = circ_confmean(iPhaseFOI(~isnan(iPhaseFOI)), 0.05);
  end
end

p = semilogx(FOI(~isnan(phaseMean)),phaseMean(~isnan(phaseMean)), 'r:o', 'LineWidth',2, 'MarkerSize',5, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
uistack(p,'bottom');
pLegend = [pLegend p];
txtLegend{numel(txtLegend)+1} = 'Mean FOI';

p = semilogx(FOI(~isnan(phaseMean) & ~isnan(phaseCI95)),...
  phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95)) + phaseCI95(~isnan(phaseMean) & ~isnan(phaseCI95)),...
  'r:o', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
uistack(p,'bottom');

p = semilogx(FOI(~isnan(phaseMean) & ~isnan(phaseCI95)),...
  phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95)) - phaseCI95(~isnan(phaseMean) & ~isnan(phaseCI95)),...
  'r:o', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
uistack(p,'bottom');
pLegend = [pLegend p];
txtLegend{numel(txtLegend)+1} = '95% Conf';

xlabel('Frequency (Hz)');
xlim([10e-3 - 0.002   10e-1 + 2])
ylabel('Phase (rad)');
legend(pLegend,txtLegend, 'Interpreter','None');
title(['Ripple rate vs total movement.' ' Significant recordings: ' num2str(recSignificant) '/' num2str(recCount)]);
set(figH, 'Name','Ripple rate vs total movement');

ylimits = ylim;
xlimits = xlim;
set(gcf,'color','w');
legend boxoff
ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 12, 4/3, 2, [0.005 0], 'out',...
    'on', 'k', {'Frequency (Hz)'}, [xlimits(1) 10e-1 + 1], xticks,...
    'on', 'k', {'Phase (rad)'}, [-0.25 ylimits(2)], [0 pi]);
ax1.XTickLabel = {'10^-^2','10^-^1','10^0'};
ax1.YTickLabel = {'0','\pi'};
hold off

figFileName = [dataDir filesep lfp2motionDir filesep 'rippleRate2motion'];
hgsave(figH, [figFileName '.fig']);

label = [1.35 1.55];
margin = [0.3 0.1];
width = 1*15-label(1)-margin(1);
height = 1*15-label(2)-margin(2);
paperSize = resizeFig(gcf, ax1, width, height, label, margin, 0);
exportFig(gcf, [figFileName '.png'],'-dpng','-r300', paperSize);

close(figH);