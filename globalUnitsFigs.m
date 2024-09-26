% A script for producing unit summary figures for all animals


clearvars -except allData

params
lists

repository = 'uol';
if strcmp(repository,'all')
  rootFolder = [dataDir filesep laDir];
  animals = animalsOI;
  xLim = [FOI(end) FOI(1)];
elseif strcmp(repository,'uol')
  rootFolder = [dataDir filesep laDir_uol];
  animals = animalsUOLOI;
  xLim = [FOI(end) FOI(1)];
elseif strcmp(repository,'allensdk')
  rootFolder = [dataDir filesep laDir_allensdk];
  animals = animalsAllensdk;
  xLim = [0.045 FOI(1)];
end
% if ~exist('allData', 'var')
%   if strcmp(repository,'all')
%     load([dataDir filesep 'allData.mat']);
%   elseif strcmp(repository,'uol')
%     load([dataDir filesep 'allData_uol.mat']);
%   elseif strcmp(repository,'allensdk')
%     load([dataDir filesep 'allData_allensdk.mat']);
%   end
% end


% User input
opt.spikingPlot = true;
opt.pupilPlot = true;
opt.mixedPlot = true;
opt.visibility = 'on';
opt.scaleColour = true;
opt.qualityCheck = false;
opt.refractCont = refractCont;
opt.cluDist = cluDist;
opt.uniformityTest = [true true true true false];
opt.modalityTest = false;
phaseCentre = pi/2;
edges = edges + phaseCentre;
if phaseCentre ~= 0
  opt.adjustAxes = true;
else
  opt.adjustAxes = false;
end

areaOI = 'lVB';
condOI = 'anaesthesia';

mainFolder = [rootFolder filesep areaOI filesep condOI];

%%
if ~exist(mainFolder, 'file')
  mkdir(mainFolder);
end
if exist([mainFolder filesep 'globalUnits.mat'], 'file')
  load([mainFolder filesep 'globalUnits.mat'])
else
  series = {};
  for animal = 1:numel(animals)
    disp(['processing animal ' animals{animal}]);
    if exist('allData', 'var')
      dataStruct = allData.(animals{animal}).dataStruct;
    else
      load([dataDir filesep animals{animal} filesep animals{animal} '.mat'])
    end
    fnsData = fieldnames(dataStruct.seriesData);
    for dbCount = 1:numel(fnsData)
      dbStruct = dataStruct.seriesData.(fnsData{dbCount});
      
      % Determine if series phase and coherence data exist
      seriesName = seriesFromEntry(dbStruct.db(dbCount).entryName);
      if ~isfield(dbStruct.popData, 'phaseCoh') || isempty(dbStruct.popData.phaseCoh)
        continue
      end
      
      % Test for exceptions
      if exceptionTest(except, seriesName)
        continue
      end
      
      % Determine if population rate > 0
      if firingRateTest(sum(dbStruct.popData.MUAsAll,1), dbStruct.conf.samplingParams.srData)
        continue
      end
      
      % Determine recording area
      [~, ~, area] = determineArea(seriesName);
      if strcmp(repository, 'allensdk') && ~strcmp(area, areaOI)
        continue
      elseif strcmp(repository, 'uol') && ~startsWith(area, areaOI)
        continue
      end
      
      % Determine recording condition (i.e., awake or anaesthesia)
      [breakClause, iCond] = series2condition(awake, anaesthesia, seriesName);
      if breakClause
        continue
      end
      if strcmp(condOI, 'awake') && iCond ~= 1
        continue
      elseif strcmp(condOI, 'anaesthesia') && iCond ~= 2
        continue
      end
      
      series{numel(series)+1} = fnsData{dbCount}; %#ok<*SAGROW>
    end
  end
  
  % Get concatenated data
  concatenatedData = concatenateSeries(dataDir, series, opt);
  save([mainFolder filesep 'globalUnits.mat'], 'concatenatedData', '-v7.3');
end


units = concatenatedData.units;
FOI = concatenatedData.FOI;
phaseFOI = zeros(numel(units),numel(FOI));
cohFOI = zeros(numel(units),numel(FOI));
coh_confFOI = zeros(numel(units),numel(FOI));
phaseFOI_pupil = zeros(numel(units),numel(FOI));
cohFOI_pupil = zeros(numel(units),numel(FOI));
coh_confFOI_pupil = zeros(numel(units),numel(FOI));
phaseFOI_half1 = zeros(numel(units),numel(FOI));
cohFOI_half1 = zeros(numel(units),numel(FOI));
coh_confFOI_half1 = zeros(numel(units),numel(FOI));
phaseFOI_half1_pupil = zeros(numel(units),numel(FOI));
cohFOI_half1_pupil = zeros(numel(units),numel(FOI));
coh_confFOI_half1_pupil = zeros(numel(units),numel(FOI));
phaseFOI_half2 = zeros(numel(units),numel(FOI));
cohFOI_half2 = zeros(numel(units),numel(FOI));
coh_confFOI_half2 = zeros(numel(units),numel(FOI));
phaseFOI_half2_pupil = zeros(numel(units),numel(FOI));
cohFOI_half2_pupil = zeros(numel(units),numel(FOI));
coh_confFOI_half2_pupil = zeros(numel(units),numel(FOI));
for u = 1:numel(units)
  if isempty(concatenatedData.pr.phaseCohFOI{u}) || isempty(concatenatedData.pr.phaseCohFOI{u}.phaseFOI)
    phaseFOI(u,:) = nan(size(phaseFOI(u,:)));
    cohFOI(u,:) = nan(size(cohFOI(u,:)));
    coh_confFOI(u,:) = nan(size(coh_confFOI(u,:)));
    phaseFOI_half1(u,:) = nan(size(phaseFOI_half1(u,:)));
    phaseFOI_half2(u,:) = nan(size(phaseFOI_half2(u,:)));
    cohFOI_half1(u,:) = nan(size(cohFOI_half1(u,:)));
    cohFOI_half2(u,:) = nan(size(cohFOI_half2(u,:)));
    coh_confFOI_half1(u,:) = nan(size(coh_confFOI_half1(u,:)));
    coh_confFOI_half2(u,:) = nan(size(coh_confFOI_half2(u,:)));
  else
    phaseFOI(u,:) = concatenatedData.pr.phaseCohFOI{u}.phaseFOI;
    cohFOI(u,:) = concatenatedData.pr.phaseCohFOI{u}.cohFOI;
    coh_confFOI(u,:) = concatenatedData.pr.phaseCohFOI{u}.coh_confFOI;
    cohFOI(u, cohFOI(u,:) - coh_confFOI(u,:) <= 0) = NaN;
    phaseFOI_half1(u,:) = concatenatedData.pr.phaseCohHalvesFOI{u}.half1.phase;
    phaseFOI_half2(u,:) = concatenatedData.pr.phaseCohHalvesFOI{u}.half2.phase;
    cohFOI_half1(u,:) = concatenatedData.pr.phaseCohHalvesFOI{u}.half1.coh;
    cohFOI_half2(u,:) = concatenatedData.pr.phaseCohHalvesFOI{u}.half2.coh;
    coh_confFOI_half1(u,:) = concatenatedData.pr.phaseCohHalvesFOI{u}.half1.coh_conf;
    coh_confFOI_half2(u,:) = concatenatedData.pr.phaseCohHalvesFOI{u}.half2.coh_conf;
    cohFOI_half1(u, cohFOI_half1(u,:) - coh_confFOI_half1(u,:) <= 0) = NaN;
    cohFOI_half2(u, cohFOI_half2(u,:) - coh_confFOI_half2(u,:) <= 0) = NaN;
  end
  if opt.pupilPlot || opt.mixedPlot
    if isempty(concatenatedData.pupil.unitData.phaseCohFOI{u}.phaseFOI)
      phaseFOI_pupil(u,:) = nan(size(phaseFOI_pupil(u,:)));
      cohFOI_pupil(u,:) = nan(size(cohFOI_pupil(u,:)));
      coh_confFOI_pupil(u,:) = nan(size(coh_confFOI_pupil(u,:)));
      cohFOI_pupil(u,:) = nan(size(cohFOI_pupil(u,:)));
      phaseFOI_half1_pupil(u,:) = nan(size(phaseFOI_half1_pupil(u,:)));
      phaseFOI_half2_pupil(u,:) = nan(size(phaseFOI_half2_pupil(u,:)));
      cohFOI_half1_pupil(u,:) = nan(size(cohFOI_half1_pupil(u,:)));
      cohFOI_half2_pupil(u,:) = nan(size(cohFOI_half2_pupil(u,:)));
      coh_confFOI_half1_pupil(u,:) = nan(size(coh_confFOI_half1_pupil(u,:)));
      coh_confFOI_half2_pupil(u,:) = nan(size(coh_confFOI_half2_pupil(u,:)));
    else
      phaseFOI_pupil(u,:) = concatenatedData.pupil.unitData.phaseCohFOI{u}.phaseFOI;
      cohFOI_pupil(u,:) = concatenatedData.pupil.unitData.phaseCohFOI{u}.cohFOI;
      coh_confFOI_pupil(u,:) = concatenatedData.pupil.unitData.phaseCohFOI{u}.coh_confFOI;
      cohFOI_pupil(u, cohFOI_pupil(u,:) - coh_confFOI_pupil(u,:) <= 0) = NaN;
      phaseFOI_half1_pupil(u,:) = concatenatedData.pupil.unitData.phaseCohHalvesFOI{u}.half1.phase;
      phaseFOI_half2_pupil(u,:) = concatenatedData.pupil.unitData.phaseCohHalvesFOI{u}.half2.phase;
      cohFOI_half1_pupil(u,:) = concatenatedData.pupil.unitData.phaseCohHalvesFOI{u}.half1.coh;
      cohFOI_half2_pupil(u,:) = concatenatedData.pupil.unitData.phaseCohHalvesFOI{u}.half2.coh;
      coh_confFOI_half1_pupil(u,:) = concatenatedData.pupil.unitData.phaseCohHalvesFOI{u}.half1.coh_conf;
      coh_confFOI_half2_pupil(u,:) = concatenatedData.pupil.unitData.phaseCohHalvesFOI{u}.half2.coh_conf;
      cohFOI_half1_pupil(u, cohFOI_half1_pupil(u,:) - coh_confFOI_half1_pupil(u,:) <= 0) = NaN;
      cohFOI_half2_pupil(u, cohFOI_half2_pupil(u,:) - coh_confFOI_half2_pupil(u,:) <= 0) = NaN;
    end
  end
end


%% SPIKING
if opt.spikingPlot
  % Phase frequency histograms and maps
  phase = recentrePhase(phaseFOI, phaseCentre);
  if ~isempty(phase)
    histPhaseOverF = zeros(numel(edges),numel(FOI));
    for f = 1:numel(FOI)
      histPhaseOverF(:,f) = [sum(isnan(phase(:,f))) histcounts(phase(:,f), edges)];
    end
    
    plotCount = 0;
    sbPhase = figure('units', 'normalized', 'position', [0.002, .04, 1, .88], 'Visible', 'on');
    for f = 1:numel(FOI)
      disp(['processing spiking frequency ' num2str(FOI(f)) ' Hz']);
      
      % Individual histos
      phaseHistPlotUnadjusted(edges, histPhaseOverF(:,f), 'Unit count', [mainFolder filesep 'SPIKING_'...
        areaOI '_' condOI '_unit_phase_at_' num2str(FOI(f)) '_Hz'], 'on', opt.adjustAxes); % Figure
      [distributionStats_spiking.pRayleigh{f}, distributionStats_spiking.zRayleigh{f}, distributionStats_spiking.pOmnibus{f},...
        distributionStats_spiking.mOmnibus{f}, distributionStats_spiking.pRao{f}, distributionStats_spiking.U_Rao{f},...
        distributionStats_spiking.pV{f}, distributionStats_spiking.vV{f}, distributionStats_spiking.pHR{f},...
        distributionStats_spiking.T_HR{f}, distributionStats_spiking.nModes{f}, distributionStats_spiking.excessMass{f},...
        distributionStats_spiking.U2_KDE{f}, distributionStats_spiking.pKDE{f}, distributionStats_spiking.modes{f},...
        distributionStats_spiking.dipHDT{f}, distributionStats_spiking.pHDT{f}] = phaseDistributionTests(...
        phase(:,f), edges, opt.modalityTest, opt.uniformityTest); % Stats
      
      % Summary histos
      if FOI(f) == 0.01 || FOI(f) == 0.03 || FOI(f) == 0.05 || FOI(f) == 0.1 || FOI(f) == 0.3 ||...
          FOI(f) == 0.5 || FOI(f) == 1 || FOI(f) == 4 || FOI(f) == 10
        figure(sbPhase);
        plotCount = plotCount + 1;
        subplot(3, 3, 10 - plotCount);
        if plotCount <= 2
          phaseHistSubplotUnadjusted(edges, histPhaseOverF(:,f), '', 'Phase (rad)', 'k', opt.adjustAxes);
          dispFOI(FOI(f), sum(histPhaseOverF(:,f)));
        elseif plotCount == 3
          phaseHistSubplotUnadjusted(edges, histPhaseOverF(:,f), 'Unit count', 'Phase (rad)', 'k', opt.adjustAxes);
          dispFOI(FOI(f), sum(histPhaseOverF(:,f)));
        elseif plotCount <= 5
          phaseHistSubplotUnadjusted(edges, histPhaseOverF(:,f), '', '', 'k', opt.adjustAxes);
          dispFOI(FOI(f), sum(histPhaseOverF(:,f)));
        elseif plotCount == 6
          phaseHistSubplotUnadjusted(edges, histPhaseOverF(:,f), 'Unit count', '', 'k', opt.adjustAxes);
          dispFOI(FOI(f), sum(histPhaseOverF(:,f)));
        elseif plotCount <= 8
          phaseHistSubplotUnadjusted(edges, histPhaseOverF(:,f), '', '', 'k', opt.adjustAxes);
          dispFOI(FOI(f), sum(histPhaseOverF(:,f)));
        elseif plotCount == 9
          phaseHistSubplotUnadjusted(edges, histPhaseOverF(:,f), 'Unit count', '', 'k', opt.adjustAxes);
          dispFOI(FOI(f), sum(histPhaseOverF(:,f)));
          set(gcf, 'Color','w');
          print(sbPhase, [mainFolder filesep 'SPIKING_' areaOI '_' condOI...
            '_unit_phase_at_p01_p03_p05_p1_p3_p5_1_4_10_Hz.png'], '-dpng', '-r300');
          hgsave([mainFolder filesep 'SPIKING_' areaOI '_' condOI...
            '_unit_phase_at_p01_p03_p05_p1_p3_p5_1_4_10_Hz.fig'])
          close(sbPhase);
        end
      end
    end
    
    % Phase frequency map
    cm = 'cold';
    h = phaseGraph(['SPIKING Phase map ' areaOI ' ' condOI], fliplr(FOI), xLim,...
      [0.01 0.1 1 10 100], edges, [edges(1) edges(end)], [-pi -pi/2 0 pi/2 pi]+phaseCentre,...
      fliplr(histPhaseOverF(2:end,:)), cm, 'on', opt.adjustAxes);
    set(gca, 'XTickLabel',{'0.01','0.1','1','10','100'});
    set(h,'color','white')
    if opt.scaleColour
      %rescaleColour
    end
    hgsave(h, [mainFolder filesep 'SPIKING_Phase_map_' areaOI '_' condOI]);
    print(h, [mainFolder filesep 'SPIKING_Phase_map_' areaOI '_' condOI '.png'],'-dpng','-r300');
    close(h);
  end
  
  % Phase correlations
  phase1 = recentrePhase(phaseFOI_half1, 0);
  phase2 = recentrePhase(phaseFOI_half2, 0);
  if ~isempty(phase1) && ~isempty(phase2)
    
    % Correlations and individual graphs
    [figPhase, rFOI_spiking, pvalFOI_spiking, nFOI_spiking] = halfCorrPlot_ca(phase1, phase2, phase, FOI, condOI, areaOI, mainFolder,...
      'Half-rec unit phase corr for ', true, 'SPIKING', '', opt);
    close all
  end
  
  % Summary subplots
  sbPhase = halfCorrSummary(phase1, phase2, FOI, rFOI_spiking, pvalFOI_spiking, nFOI_spiking, condOI, areaOI, mainFolder,...
    '_unit_half_phase_correlations_at_p01_p03_p05_p1_p3_p5_1_4_10_Hz', 'SPIKING', '', opt.adjustAxes);
  close(sbPhase)
  
  % Coherence histograms
  edges1 = 0:0.05:1;
  edges2 = 0:0.025:0.5;
  separator = 12;
  coh = cohFOI;
  histCohOverF1 = zeros(numel(edges1),numel(FOI(separator+1:end)));
  histCohOverF2 = zeros(numel(edges2),numel(FOI(1:separator)));
  for f = 1:size(histCohOverF1,2)
    histCohOverF1(:,f) = [sum(isnan(coh(:,separator+f))) histcounts(coh(:,separator+f), edges1)];
  end
  for f = 1:size(histCohOverF1,2)
    cohHistPlot(edges1, histCohOverF1(:,f), {'non-signif.','0','0.25','0.5','0.75','1'}, 'Unit count',...
      [mainFolder filesep 'SPIKING_coherence_at_' num2str(FOI(separator+f)) '_Hz'], opt.visibility);
  end
  for f = 1:size(histCohOverF2,2)
    histCohOverF2(:,f) = [sum(isnan(coh(:,f))) histcounts(coh(:,f), edges2)];
  end
  for f = 1:size(histCohOverF2,2)
    cohHistPlot(edges2, histCohOverF2(:,f), {'non-signif.','0','0.125','0.25','0.375','0.5'}, 'Unit count',...
      [mainFolder filesep 'SPIKING_coherence_at_' num2str(FOI(f)) '_Hz'], opt.visibility);
  end
end


%% PUPIL
if opt.pupilPlot
  % Phase frequency histograms and maps
  FOIind = 1:numel(FOI);
  FOIind(FOI > 2) = [];
  FOI_pupil = FOI(FOIind);
  phase = recentrePhase(phaseFOI_pupil(:,FOIind), phaseCentre);
  if ~isempty(phase)
    histPhaseOverF = zeros(numel(edges),numel(FOI_pupil));
    for f = 1:numel(FOI_pupil)
      histPhaseOverF(:,f) = [sum(isnan(phase(:,f))) histcounts(phase(:,f), edges)];
    end
    
    plotCount = 0;
    sbPhase = figure('units', 'normalized', 'position', [0.002, .04, 1, .88], 'Visible', 'on');
    for f = 1:numel(FOI_pupil)
      disp(['processing pupil frequency ' num2str(FOI_pupil(f)) ' Hz']);
      
      % Individual histos
      phaseHistPlotUnadjusted(edges, histPhaseOverF(:,f), 'Unit count', [mainFolder filesep 'PUPIL_'...
        areaOI '_' condOI '_unit_phase_at_' num2str(FOI_pupil(f)) '_Hz'], 'on', opt.adjustAxes); % Figure
      [distributionStats_pupil.pRayleigh{f}, distributionStats_pupil.zRayleigh{f}, distributionStats_pupil.pOmnibus{f},...
        distributionStats_pupil.mOmnibus{f}, distributionStats_pupil.pRao{f}, distributionStats_pupil.U_Rao{f},...
        distributionStats_pupil.pV{f}, distributionStats_pupil.vV{f}, distributionStats_pupil.pHR{f},...
        distributionStats_pupil.T_HR{f}, distributionStats_pupil.nModes{f}, distributionStats_pupil.excessMass{f},...
        distributionStats_pupil.U2_KDE{f}, distributionStats_pupil.pKDE{f}, distributionStats_pupil.modes{f},...
        distributionStats_pupil.dipHDT{f}, distributionStats_pupil.pHDT{f}] = phaseDistributionTests(...
        phase(:,f), edges, opt.modalityTest, opt.uniformityTest); % Stats
      
      % Summary histos
      if FOI_pupil(f) == 0.01 || FOI_pupil(f) == 0.03 || FOI_pupil(f) == 0.05 || FOI_pupil(f) == 0.1 || FOI_pupil(f) == 0.3 ||...
          FOI_pupil(f) == 0.5 || FOI_pupil(f) == 1 || FOI_pupil(f) == 4 || FOI_pupil(f) == 10
        figure(sbPhase);
        plotCount = plotCount + 1;
        subplot(3, 3, 10 - plotCount);
        if plotCount <= 2
          phaseHistSubplotUnadjusted(edges, histPhaseOverF(:,f), '', 'Phase (rad)', 'k', opt.adjustAxes);
          dispFOI(FOI_pupil(f), sum(histPhaseOverF(:,f)));
        elseif plotCount == 3
          phaseHistSubplotUnadjusted(edges, histPhaseOverF(:,f), 'Unit count', 'Phase (rad)', 'k', opt.adjustAxes);
          dispFOI(FOI_pupil(f), sum(histPhaseOverF(:,f)));
        elseif plotCount <= 5
          phaseHistSubplotUnadjusted(edges, histPhaseOverF(:,f), '', '', 'k', opt.adjustAxes);
          dispFOI(FOI_pupil(f), sum(histPhaseOverF(:,f)));
        elseif plotCount == 6
          phaseHistSubplotUnadjusted(edges, histPhaseOverF(:,f), 'Unit count', '', 'k', opt.adjustAxes);
          dispFOI(FOI_pupil(f), sum(histPhaseOverF(:,f)));
        elseif plotCount <= 8
          phaseHistSubplotUnadjusted(edges, histPhaseOverF(:,f), '', '', 'k', opt.adjustAxes);
          dispFOI(FOI_pupil(f), sum(histPhaseOverF(:,f)));
        elseif plotCount == 9
          phaseHistSubplotUnadjusted(edges, histPhaseOverF(:,f), 'Unit count', '', 'k', opt.adjustAxes);
          dispFOI(FOI_pupil(f), sum(histPhaseOverF(:,f)));
          set(gcf, 'Color','w');
          print(sbPhase, [mainFolder filesep 'PUPIL_' areaOI 'VsPupil_' condOI...
            '_unit_phase_at_p01_p03_p05_p1_p3_p5_1_4_10_Hz.png'], '-dpng', '-r300');
          hgsave([mainFolder filesep 'PUPIL_' areaOI 'VsPupil_' condOI...
            '_unit_phase_at_p01_p03_p05_p1_p3_p5_1_4_10_Hz.fig'])
          close(sbPhase);
        end
      end
    end
    
    % Phase frequency map
    cm = 'cold';
    h = phaseGraph(['PUPIL Phase map ' areaOI 'VsPupil ' condOI], fliplr(FOI_pupil), [xLim(1) FOI_pupil(1)],...
      [0.01 0.1 1 10 100], edges, [edges(1) edges(end)], [-pi -pi/2 0 pi/2 pi]+phaseCentre,...
      fliplr(histPhaseOverF(2:end,:)), cm, 'on');
    set(gca, 'YTickLabel',{'-\pi/2','0','\pi/2','\pi','3\pi/2'});
    set(gca, 'XTickLabel',{'0.01','0.1','1','10','100'});
    set(h,'color','white')
    hgsave(h, [mainFolder filesep 'PUPIL_Phase_map_' areaOI 'VsPupil_' condOI]);
    print(h, [mainFolder filesep 'PUPIL_Phase_map_' areaOI 'VsPupil_' condOI '.png'],'-dpng','-r300');
    close(h);
  end
  
  % Phase and coherence correlations
  phase1 = recentrePhase(phaseFOI_half1_pupil(:,FOIind), 0);
  phase2 = recentrePhase(phaseFOI_half2_pupil(:,FOIind), 0);
  if ~isempty(phase1) && ~isempty(phase2)
    
    % Correlations and individual graphs
    [figPhase, rFOI_pupil, pvalFOI_pupil, nFOI_pupil] = halfCorrPlot_ca(phase1, phase2, phase, FOI_pupil, condOI, areaOI, mainFolder,...
      'Half-rec unit phase corr for ', true, 'PUPIL', 'VsPupil', opt);
    close all
  end
  
  % Summary subplots
  sbPhase = halfCorrSummary(phase1, phase2, FOI_pupil, rFOI_pupil, pvalFOI_pupil, nFOI_pupil, condOI, areaOI, mainFolder,...
    '_unit_half_phase_correlations_at_p01_p03_p05_p1_p3_p5_1_4_10_Hz', 'PUPIL', 'VsPupil', opt.adjustAxes);
  close(sbPhase)
  
  % Coherence histograms
  edges1 = 0:0.05:1;
  edges2 = 0:0.025:0.5;
  separator = 12;
  coh = cohFOI_pupil(:,FOIind);
  histCohOverF1 = zeros(numel(edges1),numel(FOI_pupil(separator+1:end)));
  histCohOverF2 = zeros(numel(edges2),numel(FOI_pupil(1:separator)));
  for f = 1:size(histCohOverF1,2)
    histCohOverF1(:,f) = [sum(isnan(coh(:,separator+f))) histcounts(coh(:,separator+f), edges1)];
  end
  for f = 1:size(histCohOverF1,2)
    cohHistPlot(edges1, histCohOverF1(:,f), {'non-signif.','0','0.25','0.5','0.75','1'}, 'Unit count',...
      [mainFolder filesep 'PUPIL_coherence_at_' num2str(FOI_pupil(separator+f)) '_Hz'], opt.visibility);
  end
  for f = 1:size(histCohOverF2,2)
    histCohOverF2(:,f) = [sum(isnan(coh(:,f))) histcounts(coh(:,f), edges2)];
  end
  for f = 1:size(histCohOverF2,2)
    cohHistPlot(edges2, histCohOverF2(:,f), {'non-signif.','0','0.125','0.25','0.375','0.5'}, 'Unit count',...
      [mainFolder filesep 'PUPIL_coherence_at_' num2str(FOI_pupil(f)) '_Hz'], opt.visibility);
  end
end


%% MIXED
if opt.mixedPlot
  % Phase and coherence correlations
  opt.xLabel = 'Spiking';
  opt.yLabel = 'Pupil';
  FOIind = 1:numel(FOI);
  FOIind(FOI > 2) = [];
  FOI_pupil = FOI(FOIind);
  phase1 = recentrePhase(phaseFOI(:,FOIind), 0);
  phase2 = recentrePhase(phaseFOI_pupil(:,FOIind), 0);
  if ~isempty(phase1) && ~isempty(phase2)
    
    % Correlations and individual graphs
    [figPhase, rFOI_mixed, pvalFOI_mixed, nFOI_mixed] = halfCorrPlot_ca(phase1, phase2, phase1, FOI_pupil, condOI, areaOI, mainFolder,...
      'Spiking-pupil unit phase corr for ', true, 'MIXED', 'VsMixed', opt);
    close all
  end
  
  % Summary subplots
  sbPhase = halfCorrSummary(phase1, phase2, FOI_pupil, rFOI_mixed, pvalFOI_mixed, nFOI_mixed, condOI, areaOI, mainFolder,...
    '_unit_spiking_pupil_correlations_at_p01_p03_p05_p1_p3_p5_1_4_10_Hz', 'MIXED', 'VsMixed', opt.adjustAxes);
  close(sbPhase)
  
  % Coherence of population rate and pupil area
  initFigName = [mainFolder filesep 'MIXED_correlation_'];
  fH = figure('units', 'normalized', 'position', [0.002, .04, 1/2, .88/1.04], 'Visible', opt.visibility);
  for s =1:numel(concatenatedData.pupil.popData)
    scaleType = get(gca,'XScale');
    if strcmpi(scaleType,'linear')
      semilogx(concatenatedData.pupil.popData{s}.freq,...
        concatenatedData.pupil.popData{s}.coh .* concatenatedData.pupil.popData{s}.rateadjust_kappa)
    else
      hold on, semilogx(concatenatedData.pupil.popData{s}.freq,...
        concatenatedData.pupil.popData{s}.coh .* concatenatedData.pupil.popData{s}.rateadjust_kappa)
    end
  end
  title('Population_rate_coherence_with_pupil_area_for_all_series', 'Interpreter', 'none')
  xlabel('Frequency (Hz)');
  ylabel('Coherence');
  hgsave(fH, [mainFolder filesep 'MIXED_population_rate_coherence_with_pupil_area']);
  close(fH);
  
  % beta and pupil coherence correlations
  for i = 1:numel(FOI)
    figBeta(i) = figure('units', 'normalized', 'position', [0.002, .04, 1/2, .88/1.04], 'Visible', opt.visibility);
    plot(cohFOI_pupil(:,i)', concatenatedData.pr.beta', '.', 'MarkerSize',20)
    title(['Beta: ' num2str(FOI(i)) ' Hz']);
    xlabel('Coherence')
    ylabel('beta')
    [mixed_rMFR_beta, mixed_rhoMFR_beta] = corrStatesMFR(figBeta(i), [cohFOI_pupil(:,i)'; concatenatedData.pr.beta']);
  end
  
  for i = 1:numel(figBeta)
    h = figure(figBeta(i)); %#ok<*NASGU>
    figFileName = [initFigName '_beta_' num2str(FOI(i)) 'Hz'];
    figFileName = strrep(figFileName,'.','p');
    set(gcf, 'Name',figFileName);
    hgsave(gcf, figFileName);
    
    %   ax1 = findobj(h, 'type','axe');
    %   xLabel = get(get(ax1,'xlabel'),'string');
    %   xTick = get(ax1,'xtick');
    %   yLabel = get(get(ax1,'ylabel'),'string');
    %   ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 25, 4/3, 2, [0.005 0], 'out',...
    %   	'on', 'k', {xLabel}, [], xTick,...
    %     'on', 'k', {yLabel}, [], xTick);
    %   label = [4 2.6];
    %   margin = [1 1];
    %   width = 1.3*15-label(1)-margin(1);
    %   height = (1.3*15)-label(2)-margin(2);
    %   paperSize = resizeFig(h, ax1, width, height, label, margin, 0);
    %   exportFig(h, [figFileName '.png'],'-dpng','-r300', paperSize);
  end
  close(figBeta);
end


%% Save data
if opt.spikingPlot && opt.pupilPlot && opt.mixedPlot
  save([mainFolder filesep 'tests.mat'], 'distributionStats_spiking','rFOI_spiking','pvalFOI_spiking','nFOI_spiking',...
    'distributionStats_pupil','rFOI_pupil','pvalFOI_pupil','nFOI_pupil',...
    'rFOI_mixed','pvalFOI_mixed','nFOI_mixed','mixed_rMFR_beta','mixed_rhoMFR_beta', '-v7.3');
else
  disp('Distribution and correlation test statistics were not saved');
end