% Run this script to perform analyses comparing LFP and pupil area data.

fclose all;
close all
clear
clc



%% INITIALISE PARAMETERS
dataDir = 'R:\CSN\Shared\Dynamics\Data\mouse_analysis_data';
lists

animals = {
  'M190114_A_MD';
  'M190128_A_MD';
  'M190128_B_MD';
  'M190128_C_MD';
  'M190322_C_MD';
  'M190503_B_MD';
  'M190523_A_MD';
  'M190523_B_MD'};

conditions = {
  'awake'};

areas = {
  'ThVsCA1';
  'S1VsCA1';
  'RSCVsCA1';
  'VBVsCA1'};

awake = {
  '20190122191103';
  '20190205085751';
  '20190208100404';
  '20190212095950';
  '20190215112507';
  '20190208152417';
  '20190218132050';
  '20190311150750';
  '20190214143752';
  '20190312133144';
  '20190404111519';
  '20190405102743';
  '20190522131946';
  '20190601133104';
  '20190602202843';
  '20190604174736';
  '20190606142912';
  '20190617184125';
  '20190621114337';
  '20190602180606';
  '20190604103515';
  '20190606105826';
  '20190618152633';
  '20190621174558'};
anaesthesia = {};
except = {};
    


%% DISPLAY POPULATION RATE FREQUENCY PHASE PROFILES RELATIVE TO PUPIL AREA FOR ALL ANIMALS AND ALL RECORDINGS
% LOOP THROUGH ANIMALS
for animal = 1:numel(animals)
  load([dataDir filesep animals{animal} filesep 'CAR' filesep animals{animal} '.mat'])
  animalColour = animalColours(animals(animal));
  fnsData = fieldnames(dataStruct.seriesData_ca);
  
  % INITIALISE FIGURES
  fnsData2 = fieldnames(dataStruct.seriesData);
  FOI = dataStruct.seriesData.(fnsData2{1}).FOI;
  if animal == 1
    [figs, plotLegends, txtLegends, figEmptyAnimal] = initFigs(conditions, areas);
    [figsC, plotLegendsC, txtLegendsC, figEmptyAnimalC] = initFigs(conditions, areas);
    areaRecCount = {};
    areaRecSignificant = {};
    areaRecSignificantFOI = {};
    areaCohFOI = {};
    areaCohFOIindividual = {};
    areaPhaseFOIindividual = {};
    for iCond = 1:numel(conditions)
      areaRecCountCond = {};
      areaRecSignificantCond = {};
      areaRecSignificantFOIcond = {};
      areaCohFOIcond = {};
      areaCohFOIindividualCond = {};
      areaPhaseFOIindividualCond = {};
      for iArea = 1:numel(areas)
        areaRecCountCond{iArea} = 0;
        areaRecSignificantCond{iArea} = 0;
        areaRecSignificantFOIcond{iArea} = zeros(1,numel(FOI));
        areaCohFOIcond{iArea} = zeros(1,numel(FOI));
        areaCohFOIindividualCond{iArea} = {};
        areaPhaseFOIindividualCond{iArea} = {};
      end
      areaRecCount{iCond} = areaRecCountCond;
      areaRecSignificant{iCond} = areaRecSignificantCond;
      areaRecSignificantFOI{iCond} = areaRecSignificantFOIcond;
      areaCohFOI{iCond} = areaCohFOIcond;
      areaCohFOIindividual{iCond} = areaCohFOIindividualCond;
      areaPhaseFOIindividual{iCond} = areaPhaseFOIindividualCond;
    end
  else
    figEmptyAnimal = {};
    figEmptyAnimalC = {};
    for iCond = 1:numel(conditions)
      figEmptyCond = {};
      for iArea = 1:numel(areas)
        figEmptyCond{iArea} = true; %#ok<*SAGROW>
      end
      figEmptyAnimal{iCond} = figEmptyCond;
      figEmptyAnimalC{iCond} = figEmptyCond;
    end
  end
  
  % INITIALISE COHERENCE MATRICES
  if animal == 1
    recMatInit = NaN(7,7,numel(FOI));
    cohMatInit = NaN(7,7,numel(FOI));
    cohConfMatInit = NaN(7,7,numel(FOI));
    mfrPR1MatInit = NaN(7,7,numel(FOI));
    mfrPR2MatInit = NaN(7,7,numel(FOI));
    for iCond = 1:numel(conditions)
      recMat{iCond} = {};
      cohMat{iCond} = {};
      cohConfMat{iCond} = {};
      mfrPR1Mat{iCond} = {};
      mfrPR2Mat{iCond} = {};
      seriesMat{iCond} = {};
      series1Mat{iCond} = {};
      series2Mat{iCond} = {};
      animalMat{iCond} = {};
    end
  end
  
  % LOOP THROUGH DB ENTRIES
  for dbCount = 1:numel(fnsData)
    dbStruct = dataStruct.seriesData_ca.(fnsData{dbCount});
      
    % DETERMINE IF SERIES PHASE AND COHERENCE DATA EXIST
    if exist('seriesName1', 'var')
      prevRec = seriesName1(1:14);
    else
      prevRec = '';
    end
    [seriesName1, seriesName2] = seriesNames(fnsData{dbCount});
    rec = seriesName1(1:14);
    if ~strcmpi(prevRec, rec)
      expandMatrices = true;
    end
    if ~isfield(dbStruct, 'lfpphaseCohDataPR') || isempty(dbStruct.lfpphaseCohDataPR)
      continue
    end
    
    % TEST FOR EXCEPTIONS
    for iExcept = 1:numel(except)
      if strcmpi(seriesName1, except{iExcept})
        breakClause = true;
        break
      elseif strcmpi(seriesName2, except{iExcept})
        breakClause = true;
        break
      else
        breakClause = false;
      end
    end
    if exist('breakClause', 'var') && breakClause
      continue
    end
    
    % DETERMINE WHICH FIGURES OWN THE SERIES
    areaID1 = seriesName1(15:end);
    areaID2 = seriesName2(15:end);
    if strcmpi(areaID1, '24') && strcmpi(areaID2, '6')
      area = 1;
    elseif strcmpi(areaID1, '1') && strcmpi(areaID2, '6')
      area = 2;
    elseif strcmpi(areaID1, '7') && strcmpi(areaID2, '6')
      area = 3;
    elseif strcmpi(areaID1, '2') && strcmpi(areaID2, '6')
      area = 4;
    else
      continue
    end
    for entry = 1:numel(awake)
      if strcmpi(awake{entry}, seriesName1(1:14)) && strcmpi(awake{entry}, seriesName2(1:14))
        condition = 1;
        iCh = ch{entry};
      end
    end
    for entry = 1:numel(anaesthesia)
      if strcmpi(anaesthesia{entry}, seriesName1(1:14)) && strcmpi(anaesthesia{entry}, seriesName2(1:14))
        condition = 2;
      end
    end
    
    % ASSIGN MATRICES
    areaRecCount{condition}{area} = areaRecCount{condition}{area} + 1;
    if expandMatrices
      recMat{condition}{numel(recMat{condition})+1} = recMatInit;
      cohMat{condition}{numel(cohMat{condition})+1} = cohMatInit;
      cohConfMat{condition}{numel(cohConfMat{condition})+1} = cohConfMatInit;
      seriesMat{condition}{numel(seriesMat{condition})+1} = fnsData{dbCount};
      series1Mat{condition}{numel(series1Mat{condition})+1} = seriesName1;
      series2Mat{condition}{numel(series2Mat{condition})+1} = seriesName2;
      animalMat{condition}{numel(animalMat{condition})+1} = animals{animal};
      expandMatrices = false;
    end
    
    % LOAD THE CONTENTS OF THE DB STRUCTURE VARIABLE
    % Phase
    phase = bestUnwrap(dbStruct.lfpphaseCohDataPR.rippleRate{iCh}.phase);
    phase(isnan(dbStruct.lfpphaseCohDataPR.rippleRate{iCh}.phase_confU)) = NaN;
    phase(isnan(dbStruct.lfpphaseCohDataPR.rippleRate{iCh}.phase_confL)) = NaN;
    phase(isnan(dbStruct.lfpphaseCohDataPR.rippleRate{iCh}.coh_conf)) = NaN;
    phase = -phase;
    freq = dbStruct.lfpphaseCohDataPR.rippleRate{iCh}.freq;
    [~, endFreq] = min(abs(freq - 2));
    phase1 = phase(1:endFreq);
    phase2 = phase(endFreq:end);
    if sum(isnan(phase1)) > 0.6*numel(phase1)
      continue
    end
    if mean(phase1, 'omitnan') > 1.5*pi || mean(phase2, 'omitnan') > 1.5*pi
      phase = phase - 2*pi;
    elseif mean(phase1, 'omitnan') < -1.5*pi || mean(phase2, 'omitnan') < -1.5*pi
      phase = phase + 2*pi;
    end
    
    % Coherence
    coh = dbStruct.lfpphaseCohDataPR.rippleRate{iCh}.coh;
    coh(isnan(dbStruct.lfpphaseCohDataPR.rippleRate{iCh}.coh_conf)) = NaN;
    rateadjust_kappa = dbStruct.lfpphaseCohDataPR.rippleRate{iCh}.rateadjust_kappa;
    coh_confL = (coh - dbStruct.lfpphaseCohDataPR.rippleRate{iCh}.coh_conf) .* rateadjust_kappa;
    coh = coh .* rateadjust_kappa;
    coh_confL(coh_confL <= 0) = NaN;
    coh(isnan(coh_confL)) = NaN;
    coh1 = coh(1:endFreq);
    if sum(isnan(coh1)) > 0.6*numel(coh1)
      continue
    end
    
    % COUNT RECORDINGS WITH SIGNIFICANT VALUES FOR A PARTICULAR CONDITION AND AREA
    areaRecSignificant{condition}{area} = areaRecSignificant{condition}{area} + 1;
    [phaseFOIinit, cohFOIinit, confFOI] = phaseCohFOI(FOI, freq, phase, coh, coh_confL, ones(1,numel(coh)));
    areaRecSignificantFOI{condition}{area} = areaRecSignificantFOI{condition}{area} + ~isnan(cohFOIinit);
    cohFOI = cohFOIinit;
    cohFOI(isnan(cohFOI)) = 0;
    areaCohFOI{condition}{area} = areaCohFOI{condition}{area} + cohFOI;
    if isempty(areaPhaseFOIindividual{condition}{area})
      areaPhaseFOIindividual{condition}{area} = phaseFOIinit;
    else
      areaPhaseFOIindividual{condition}{area} = [areaPhaseFOIindividual{condition}{area}; phaseFOIinit];
    end
    if isempty(areaCohFOIindividual{condition}{area})
      areaCohFOIindividual{condition}{area} = cohFOIinit;
    else
      areaCohFOIindividual{condition}{area} = [areaCohFOIindividual{condition}{area}; cohFOIinit];
    end
    
    % DRAW FIGURES FOR A PARTICULAR CONDITION AND AREA
    [figEmptyAnimal{condition}{area}, plotLegends{condition}{area}, txtLegends{condition}{area}] = plotFig(figs{condition}{area},...
      freq, phase, animalColour, figEmptyAnimal{condition}{area}, plotLegends{condition}{area}, txtLegends{condition}{area}, animals{animal});
    [figEmptyAnimalC{condition}{area}, plotLegendsC{condition}{area}, txtLegendsC{condition}{area}] = plotFig(figsC{condition}{area},...
      freq, coh, animalColour, figEmptyAnimalC{condition}{area}, plotLegendsC{condition}{area}, txtLegendsC{condition}{area}, animals{animal});
    
    % PRODUCE COHERENCE MATRICES FOR SPECIFIC CONDITIONS
    iRow = areaEntry(areaID1);
    iCol = areaEntry(areaID2);
    for iF = 1:numel(FOI)
      recMat{condition}{numel(recMat{condition})}(iCol,iRow,iF) = 1;
      recMat{condition}{numel(recMat{condition})}(iRow,iCol,iF) = 1;
      cohMat{condition}{numel(cohMat{condition})}(iCol,iRow,iF) = cohFOIinit(iF);
      cohMat{condition}{numel(cohMat{condition})}(iRow,iCol,iF) = cohFOIinit(iF);
      cohConfMat{condition}{numel(cohConfMat{condition})}(iCol,iRow,iF) = confFOI(iF);
      cohConfMat{condition}{numel(cohConfMat{condition})}(iRow,iCol,iF) = confFOI(iF);
      for iDiagonal = 1:size(cohMat{condition}{numel(cohMat{condition})},1)
        cohMat{condition}{numel(cohMat{condition})}(iDiagonal,iDiagonal,iF) = 1;
        cohConfMat{condition}{numel(cohConfMat{condition})}(iDiagonal,iDiagonal,iF) = 1;
      end
    end
  end
end

% CALCULATE MEAN OF MATRICES
for iCond = 1:numel(conditions)
  cohMatMean{iCond} = zeros(7,7,numel(FOI));
  tallyMat{iCond} = zeros(7,7,numel(FOI));
  cohConfMatMean{iCond} = zeros(7,7,numel(FOI));
  for iRec = 1:numel(cohMat{iCond})
    cohMatRec = cohMat{iCond}{iRec};
    cohMatRec(isnan(cohMatRec)) = 0;
    cohMatMean{iCond} = cohMatMean{iCond} + cohMatRec;
    tallyMat{iCond} = tallyMat{iCond} + ~isnan(cohMat{iCond}{iRec});
    cohConfMatRec = cohConfMat{iCond}{iRec};
    cohConfMatRec(isnan(cohConfMatRec)) = 0;
    cohConfMatMean{iCond} = cohConfMatMean{iCond} + cohConfMatRec;
  end
  cohMatMean{iCond} = cohMatMean{iCond} ./ tallyMat{iCond};
  cohConfMatMean{iCond} = cohConfMatMean{iCond} ./ tallyMat{iCond};
end

% SAVE MATRICES
comparisonAreas = {'S1', 'VB', 'Th', 'DG', 'CA1', 'Hp', 'RSC'};
areaRecCountMeaning = {'Th_vs_CA1', 'S1_vs_CA1', 'RSC_vs_CA1', 'VB_vs_CA1'};
save('rippleRate2areaCohMats', 'cohMatMean','cohConfMatMean','tallyMat','recMat','cohMat','cohConfMat','mfrPR1Mat','mfrPR2Mat',...
  'seriesMat','series1Mat','series2Mat','animalMat','conditions','areaRecCount','areaRecCountMeaning','FOI','comparisonAreas', '-v7.3');

% UPDATE AND SAVE THE FIGURES
% Phase
for iCond = 1:numel(conditions)
  for iArea = 1:numel(areas)
    figure(figs{iCond}{iArea}); hold on
    if ~isempty(areaPhaseFOIindividual{iCond}{iArea})
      phaseMean = NaN(1,numel(areaRecSignificantFOI{iCond}{iArea}));
      phaseStd = NaN(1,numel(areaRecSignificantFOI{iCond}{iArea}));
      phaseCI95 = NaN(1,numel(areaRecSignificantFOI{iCond}{iArea}));
      for i = 1:numel(areaRecSignificantFOI{iCond}{iArea})
        if sum(~isnan(areaPhaseFOIindividual{iCond}{iArea}(:,i)))
          iPhaseFOI = areaPhaseFOIindividual{iCond}{iArea}(:,i);
          phaseMean(i) = circmean(iPhaseFOI(~isnan(iPhaseFOI)));
          phaseStd(i) = circ_std(iPhaseFOI(~isnan(iPhaseFOI)), [], [], 1);
          phaseCI95(i) = circ_confmean(iPhaseFOI(~isnan(iPhaseFOI)), 0.05);
        end
      end
      phaseCI95(phaseCI95 == 0) = NaN;
      
      p = semilogx(FOI(~isnan(phaseMean)),phaseMean(~isnan(phaseMean)), 'r:o', 'LineWidth',2, 'MarkerSize',5, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
      uistack(p,'bottom');
      plotLegends{iCond}{iArea} = [plotLegends{iCond}{iArea} p];
      txtLegends{iCond}{iArea}{numel(txtLegends{iCond}{iArea})+1} = 'Mean FOI';
      
      pC1 = semilogx(FOI(~isnan(phaseMean)),...
        phaseMean(~isnan(phaseMean)) + phaseCI95(~isnan(phaseMean)),...
        'r:o', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
      uistack(pC1,'bottom');
      
      pC2 = semilogx(FOI(~isnan(phaseMean)),...
        phaseMean(~isnan(phaseMean)) - phaseCI95(~isnan(phaseMean)),...
        'r:o', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
      uistack(pC2,'bottom');
      plotLegends{iCond}{iArea} = [plotLegends{iCond}{iArea} pC2];
      txtLegends{iCond}{iArea}{numel(txtLegends{iCond}{iArea})+1} = '95% conf';
    end
    hold off
    xlabel('Frequency (Hz)');
    xlim([10e-3 - 0.002   10e2])
    ylabel('Phase (rad)');
    %ylim([-(6/5)*pi (6/5)*pi])
    legend(plotLegends{iCond}{iArea},txtLegends{iCond}{iArea}, 'Interpreter','None');
    strSep = strfind(areas{iArea}, 'Vs');
    areaReordered = [areas{iArea}(strSep+2:end) 'Vs' areas{iArea}(1:strSep-1)];
    title([areaReordered ': ' conditions{iCond}...
      '   Significant recordings: ' num2str(areaRecSignificant{iCond}{iArea}) '/' num2str(areaRecCount{iCond}{iArea})])
    
    ylimits = ylim;
    xlimits = xlim;
    set(gcf,'color','w');
    legend boxoff
    ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 12, 4/3, 2, [0.005 0], 'out',...
      'on', 'k', {'Frequency (Hz)'}, xlimits, xticks,...
      'on', 'k', {'Phase (rad)'}, ylimits, [- pi 0 pi]);
    ax1.YTickLabel = {'-\pi','0','\pi'};
    
    set(gcf, 'Name',[areaReordered ': ' conditions{iCond}]);
    figFileName = [areaReordered '_' conditions{iCond} '_phase'];
    hgsave(gcf, figFileName);
    
    label = [1.35 1.55];
    margin = [0.3 0.1];
    width = 1*15-label(1)-margin(1);
    height = 1*15-label(2)-margin(2);
    paperSize = resizeFig(gcf, ax1, width, height, label, margin, 0);
    exportFig(gcf, [figFileName '.png'],'-dpng','-r300', paperSize);
  end
end

% Coherence
for iCond = 1:numel(conditions)
  for iArea = 1:numel(areas)
    figure(figsC{iCond}{iArea}); hold on
    if ~isempty(areaCohFOIindividual{iCond}{iArea})
      areaCohFOI{iCond}{iArea} = areaCohFOI{iCond}{iArea} ./ areaRecSignificantFOI{iCond}{iArea};
      p = semilogx(FOI,areaCohFOI{iCond}{iArea}, 'r:o', 'LineWidth',2, 'MarkerSize',5, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
      uistack(p,'bottom');
      for iF = 1:numel(FOI)
        text(FOI(iF),0.9, [num2str(areaRecSignificantFOI{iCond}{iArea}(iF)) '/' num2str(areaRecCount{iCond}{iArea})])
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
      pC1 = semilogx(FOI,areaCohFOI{iCond}{iArea}+cohCI95(1,:), 'r:o', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
      uistack(pC1,'bottom');
      pC2 = semilogx(FOI,areaCohFOI{iCond}{iArea}+cohCI95(2,:), 'r:o', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
      uistack(pC2,'bottom');
      plotLegendsC{iCond}{iArea} = [plotLegendsC{iCond}{iArea} pC2];
      txtLegendsC{iCond}{iArea}{numel(txtLegendsC{iCond}{iArea})+1} = '95% conf';
    end
    hold off
    xlabel('Frequency (Hz)');
    xlim([10e-3 - 0.002   10e1 + 200])
    ylabel('Coherence');
    ylim([0 1])
    legend(plotLegendsC{iCond}{iArea},txtLegendsC{iCond}{iArea}, 'Interpreter','None');
    areaReordered = [areas{iArea}(strSep+2:end) 'Vs' areas{iArea}(1:strSep-1)];
    title([areaReordered ': ' conditions{iCond}...
      '   Significant recordings: ' num2str(areaRecSignificant{iCond}{iArea}) '/' num2str(areaRecCount{iCond}{iArea})])
    set(gcf, 'Name',[areaReordered ': ' conditions{iCond}]);
    hgsave(gcf, [areaReordered '_' conditions{iCond} '_coherence']);
  end
end
close all
    

function [figs, plotLegends, txtLegends, figEmptyAnimal] = initFigs(conditions, areas)

figs = {};
plotLegends = {};
txtLegends = {};
figEmptyAnimal = {};
for iCond = 1:numel(conditions)
  figCond = {};
  plotCond = {};
  txtCond = {};
  figEmptyCond = {};
  for iArea = 1:numel(areas)
    figCond{iArea} = figure; %#ok<*AGROW>
    semilogx([10e-4 10e2],[0 0], 'k:');
    plotCond{iArea} = [];
    txtCond{iArea} = {};
    figEmptyCond{iArea} = true;
  end
  figs{iCond} = figCond;
  plotLegends{iCond} = plotCond;
  txtLegends{iCond} = txtCond;
  figEmptyAnimal{iCond} = figEmptyCond;
end
end

function [figDrawingIndicator, plotLegend, txtLegend] = plotFig(figH,...
  freq, values, animalColour, figDrawingIndicator, plotLegend, txtLegend, animal)

figure(figH); hold on
p = semilogx(freq,values, 'Color',animalColour);
if figDrawingIndicator
  if isempty(plotLegend)
    plotLegend = p;
  else
    plotLegend = [plotLegend p];
  end
  txtLegend{numel(txtLegend)+1} = animal;
  figDrawingIndicator = false;
end
hold off
end

function entry = areaEntry(areaStr)

if strcmpi(areaStr, '1')
  entry = 1;
elseif strcmpi(areaStr, '2')
  entry = 2;
elseif strcmpi(areaStr, '24')
  entry = 3;
elseif strcmpi(areaStr, '5')
  entry = 4;
elseif strcmpi(areaStr, '6')
  entry = 5;
elseif strcmpi(areaStr, '56')
  entry = 6;
elseif strcmpi(areaStr, '7')
  entry = 7;
end
end