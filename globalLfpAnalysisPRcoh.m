% Run this script to perform analyses comparing LFP and pupil area data.

fclose all;
close all
clear
clc



%% INITIALISE PARAMETERS
dataDir = 'R:\CSN\Shared\Dynamics\Data\mouse_analysis_data';
lists

animals = {
  'M180313';
  'M180711_MD';
  'M180712_MD';
  'M180713_MD';
  'M180717_B_MD';
  'M180719_MD';
  'M180924_A_MD';
  'M180924_B_MD';
  'M181112_B_MD';
  'M190114_A_MD';
  'M190128_A_MD';
  'M190128_B_MD';
  'M190128_C_MD';
  'M190322_B_MD';
  'M190322_C_MD';
  'M190503_A_MD';
  'M190503_B_MD';
  'M190523_A_MD';
  'M190523_B_MD'};



%% DISPLAY RIPPLE RATE FREQUENCY COHERENCE PROFILES FOR ALL ANIMALS AND ALL RECORDINGS
% LOOP THROUGH DB ENTRIES
recCount = 0;
recSignificant = 0;
for dbCount = 1:numel(fullSeries)
  strSep = strfind(fullSeries{dbCount},'s');
  series = fullSeries{dbCount}(strSep+1:end);
  animal = fullSeries{dbCount}(1:strSep-2);
  %animalColour = animalColours(animal);
  animalColour = [0.85 0.85 0.85];
  if dbCount == 1 || ~strcmpi(animal,prevAnimal)
    load([dataDir filesep animal filesep 'CAR' filesep animal '.mat'])
    updateLegend = true;
  end
  prevAnimal = animal;
  
% LOAD THE CONTENTS OF THE DB STRUCTURE VARIABLE
  dbStruct = dataStruct.seriesData.(fullSeries{dbCount});
  if dbCount == 1
    FOI = dbStruct.FOI;
    FOI(FOI > 2) = [];
    recSignificantFOI = zeros(1,numel(FOI));
    totCohFOI = [];
  end
  if ~isfield(dbStruct,'lfpphaseCohDataPR')
    disp(['No coherence data for ' fullSeries{dbCount} '. Skipping to the next db entry...']);
    continue
  end
  recCount = recCount + 1;
  
% PLOT THE DATA
  coh = dbStruct.lfpphaseCohDataPR.rippleRate{ch{dbCount}}.coh;
  coh(isnan(dbStruct.lfpphaseCohDataPR.rippleRate{ch{dbCount}}.coh_conf)) = NaN;
  rateadjust_kappa = dbStruct.lfpphaseCohDataPR.rippleRate{ch{dbCount}}.rateadjust_kappa;
  coh = coh .* rateadjust_kappa;
  coh_confL = (coh - dbStruct.lfpphaseCohDataPR.rippleRate{ch{dbCount}}.coh_conf) .* rateadjust_kappa;
  coh_confL(coh_confL <= 0) = NaN;
  coh(isnan(coh_confL)) = NaN;
  if sum(isnan(coh)) == numel(coh)
    continue
  end
  freq = dbStruct.lfpphaseCohDataPR.rippleRate{ch{dbCount}}.freq;
  [~, endFreq] = min(abs(freq - 2));
%   coh1 = coh(1:endFreq);
%   if sum(isnan(coh1)) > 0.6*numel(coh1)
%     continue
%   end
  recSignificant = recSignificant + 1;
  if ~exist('figH','var')
    figH = figure('units', 'normalized', 'position', [0.002, .04, 1, .88], 'Visible', 'on');
    semilogx([10e-4 10e2],[0 0], 'k:'); hold on
    p = semilogx(freq,coh, 'Color',animalColour);
    pLegend = p;
    txtLegend = {animal};
    updateLegend = false;
  else
    if updateLegend
      p = semilogx(freq,coh, 'Color',animalColour);
      pLegend = [pLegend p];
      txtLegend{numel(txtLegend)+1} = animal;
      updateLegend = false;
    else
      semilogx(freq,coh, 'Color',animalColour);
    end
  end
  
  [~, cohFOI] = phaseCohFOI(FOI, freq, ones(1,numel(coh)), coh, coh_confL, ones(1,numel(coh)));
  recSignificantFOI = recSignificantFOI + ~isnan(cohFOI);
  if isempty(totCohFOI)
    totCohFOI = cohFOI;
  else
    totCohFOI = [totCohFOI; cohFOI]; %#ok<*AGROW>
  end
end


cohMean = sum(totCohFOI, 'omitnan') ./ recSignificantFOI;
p = semilogx(FOI,cohMean, 'r:o', 'LineWidth',2, 'MarkerSize',5, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
uistack(p,'bottom');
% for iF = 1:numel(FOI)
%   text(FOI(iF),0.9, [num2str(recSignificantFOI(iF)) '/' num2str(recCount)])
% end
pLegend = [pLegend p];
txtLegend{numel(txtLegend)+1} = 'Mean FOI';

cohStd = std(totCohFOI, 'omitnan');
cohSEM = cohStd ./ sqrt(recSignificantFOI);
CI95 = zeros(2,numel(recSignificantFOI));
cohCI95 = zeros(2,numel(recSignificantFOI));
for i = 1:numel(recSignificantFOI)
  CI95(:,i) = (tinv([0.025 0.975], recSignificantFOI(i)-1))';
  cohCI95(:,i) = bsxfun(@times, cohSEM(i), CI95(:,i));
end
pC1 = semilogx(FOI,cohMean+cohCI95(1,:), 'r:o', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
uistack(pC1,'bottom');
pC2 = semilogx(FOI,cohMean+cohCI95(2,:), 'r:o', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
uistack(pC1,'bottom');
pLegend = [pLegend pC2];
txtLegend{numel(txtLegend)+1} = '95% conf';
hold off
xlabel('Frequency (Hz)');
xlim([10e-3 - 0.002   10e-1 + 2])
ylabel('Coherence');
ylim([0 1]);
legend(pLegend,txtLegend, 'Interpreter','None');
title(['Ripple rate coherence with population firing rate.' ' Significant recordings: ' num2str(recSignificant) '/' num2str(recCount)]);

ylimits = ylim;
xlimits = xlim;
set(gcf,'color','w');
legend boxoff
ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 12, 4/3, 2, [0.005 0], 'out',...
  'on', 'k', {'Frequency (Hz)'}, xlimits, xticks,...
  'on', 'k', {'Coherence'}, ylimits, [0, 0.2, 0.4, 0.6, 0.8, 1]);
ax1.YTickLabel = {'0','0.2','0.4','0.6','0.8','1'};

set(gcf, 'Name','Ripple rate coherence with population firing rate');
figFileName = 'rippleRate2muaCoh';
hgsave(figH, figFileName);

label = [1.35 1.55];
margin = [0.3 0.1];
width = 1*15-label(1)-margin(1);
height = 1*15-label(2)-margin(2);
paperSize = resizeFig(gcf, ax1, width, height, label, margin, 0);
exportFig(gcf, [figFileName '.png'],'-dpng','-r300', paperSize);

close(figH);