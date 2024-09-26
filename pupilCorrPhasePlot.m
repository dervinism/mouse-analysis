% clearvars -except repository subpop reverse qualityCheck allData
clearvars -except allData


%% INITIALISE PARAMETERS
params
lists

if ~exist('repository', 'var')
  repository = 'uol';
end
if ~exist('subpop', 'var')
  subpop = 'negative';
end
if strcmp(repository,'all')
  if strcmp(subpop, 'positive')
    rootFolder = [dataDir filesep area2pupilDir_positive];
  elseif strcmp(subpop, 'negative')
    rootFolder = [dataDir filesep area2pupilDir_negative];
  end
  animals = animalsOI;
  xLim = [FOI(end) FOI(1)];
elseif strcmp(repository,'uol')
  if strcmp(subpop, 'positive')
    rootFolder = [dataDir filesep area2pupilDir_uol_positive];
  elseif strcmp(subpop, 'negative')
    rootFolder = [dataDir filesep area2pupilDir_uol_negative];
  end
  animals = animalsUOLOI;
  xLim = [FOI(end) FOI(1)];
elseif strcmp(repository,'allensdk')
  if strcmp(subpop, 'positive')
    rootFolder = [dataDir filesep area2pupilDir_allensdk_positive];
  elseif strcmp(subpop, 'negative')
    rootFolder = [dataDir filesep area2pupilDir_allensdk_negative];
  end
  animals = animalsAllensdk;
  xLim = [0.045 FOI(1)];
  conditions = {'awake'};
end

fullRun = true;
if ~exist('qualityCheck', 'var')
  qualityCheck = true;
end
drawPhaseProfiles = true;
drawPhaseHistos = [false true true];

phaseCentre = pi/2;
edges = edges + phaseCentre;
if phaseCentre ~= 0
  opt.adjustAxes = true;
else
  opt.adjustAxes = false;
end

if qualityCheck
  mainFolder = [rootFolder filesep qualityUnitsFolder];
else
  mainFolder = [rootFolder filesep unitsFolder];
end


%% COMPUTE VARIABLES AND INITIALISE FIGURES NEEDED FOR DISPLAYING UNIT PHASE AND COHERENCE FREQUENCY PROFILES
if fullRun
  for animal = 1:numel(animals) % Loop through animals
    disp(['processing animal ' animals{animal}]);
    if exist('allData', 'var')
      dataStruct = allData.(animals{animal}).dataStruct;
    else
      load([dataDir filesep animals{animal} filesep animals{animal} '.mat'])
    end
    fnsData = fieldnames(dataStruct.seriesData);
    
    % Initialise variables
    FOI = dataStruct.seriesData.(fnsData{1}).conf.FOI;
    FOIind = 1:numel(FOI);
    FOIind(FOI > 2) = [];
    FOI(FOI > 2) = [];
    xLim = [xLim(1) FOI(1)];
    if animal == 1
      areaCohFOIindividual = {};
      areaCohConfFOIindividual = {};
      areaPhaseFOIindividual = {};
      for iCond = 1:numel(conditions)
        areaCohFOIindividualCond = {};
        areaCohConfFOIindividualCond = {};
        areaPhaseFOIindividualCond = {};
        for iArea = 1:numel(areas)
          areaCohFOIindividualCond{iArea} = {}; %#ok<*SAGROW>
          areaCohConfFOIindividualCond{iArea} = {};
          areaPhaseFOIindividualCond{iArea} = {};
        end
        areaCohFOIindividual{iCond} = areaCohFOIindividualCond;
        areaCohConfFOIindividual{iCond} = areaCohConfFOIindividualCond;
        areaPhaseFOIindividual{iCond} = areaPhaseFOIindividualCond;
      end
    end
    
    for dbCount = 1:numel(fnsData) % Loop through db entries
      dbStruct = dataStruct.seriesData.(fnsData{dbCount});
      seriesName = seriesFromEntry(fnsData{dbCount});
      
      % Determine if series pupil data exist
      if isempty(dbStruct.popData)
        continue
      end
      if ~isfield(dbStruct.popData, 'pupil') || (isfield(dbStruct.popData, 'pupil') && isempty(dbStruct.popData.pupil))
        continue
      end
      if isempty(dbStruct.popData.pupil.popData) || isempty(dbStruct.popData.pupil.popData.phase)
        continue
      end
      if strcmp(subpop, 'positive')
        dbStructSubpop = dataStruct.seriesData_positive.(fnsData{dbCount});
      elseif strcmp(subpop, 'negative')
        dbStructSubpop = dataStruct.seriesData_negative.(fnsData{dbCount});
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
      if strcmp(repository,'all')
        error('Only allensdk and uol repositories are supported currently.');
      elseif strcmp(repository,'uol')
        [~, ~, ~, ~, ~, area] = determineArea(seriesName);
      elseif strcmp(repository,'allensdk')
        area = determineArea(seriesName);
      end
      
      % Determine recording condition (i.e., awake or anaesthesia)
      [breakClause, iCond] = series2condition(awake, anaesthesia, seriesName);
      if breakClause
        continue
      end
      
      % Disqualify low quality units if needed
      units = [];
      for sh = 1:numel(fieldnames(dbStruct.shankData))
        units = [units; dbStruct.shankData.(['shank' num2str(sh)]).units]; %#ok<*AGROW>
      end
      if qualityCheck
        unitMetadata = [];
        for sh = 1:numel(fieldnames(dbStruct.shankData))
          unitMetadata = [unitMetadata; dbStruct.shankData.(['shank' num2str(sh)]).unitMetadata];
        end
        [~, qualityUnitInd] = qualityTest2(unitMetadata, cluDist, refractCont, false);
      else
        qualityUnitInd = 1:numel(units);
      end
      if isempty(qualityUnitInd)
        continue
      end
      
      for iAreaPlusAll = area % Loop through the main and pooled areas
        if numel(conditions) > 1
          condLoop = [iCond numel(conditions)];
        else
          condLoop = iCond;
        end
        for iCondPlusAll = condLoop % Loop through the main and pooled conditions
          close all
          % Load and store phase of units
          units = dbStruct.shankData.shank1.units;
          unitsSubpop = dbStructSubpop.shankData.shank1.units;
          inds1 = ismember(units, unitsSubpop)';
          inds2 = false(1,numel(units));
          inds2(qualityUnitInd) = true;
          phase = dbStruct.popData.pupil.unitData.phaseFOI(inds1 & inds2, FOIind);
          if isempty(areaPhaseFOIindividual{iCondPlusAll}{iAreaPlusAll})
            areaPhaseFOIindividual{iCondPlusAll}{iAreaPlusAll} = phase;
          else
            areaPhaseFOIindividual{iCondPlusAll}{iAreaPlusAll} = [areaPhaseFOIindividual{iCondPlusAll}{iAreaPlusAll}; phase];
          end
          
          % Load and store coherence
          coh = dbStruct.popData.pupil.unitData.cohFOI(inds1 & inds2, FOIind);
          cohConf = dbStruct.popData.pupil.unitData.coh_confFOI(inds1 & inds2, FOIind);
          if isempty(areaCohFOIindividual{iCondPlusAll}{iAreaPlusAll})
            areaCohFOIindividual{iCondPlusAll}{iAreaPlusAll} = coh;
            areaCohConfFOIindividual{iCondPlusAll}{iAreaPlusAll} = cohConf;
          else
            areaCohFOIindividual{iCondPlusAll}{iAreaPlusAll} = [areaCohFOIindividual{iCondPlusAll}{iAreaPlusAll}; coh];
            areaCohConfFOIindividual{iCondPlusAll}{iAreaPlusAll} = [areaCohConfFOIindividual{iCondPlusAll}{iAreaPlusAll}; cohConf];
          end
        end
      end
    end
  end
end

% Determine the file name and either save or load the data
if qualityCheck
  filename = [mainFolder filesep 'globalUnits_area2pupil_quality.mat'];
else
  filename = [mainFolder filesep 'globalUnits_area2pupil.mat']; %#ok<*UNRCH>
end
if fullRun
  if ~exist(mainFolder, 'file')
    mkdir(mainFolder);
  end
  save(filename, 'conditions','areas','FOI','areaPhaseFOIindividual','areaCohFOIindividual','areaCohConfFOIindividual');
else
  load(filename);
end


%% PLOT UNIT MEAN PHASE FREQUENCY PROFILES
if drawPhaseProfiles
  if ~exist([mainFolder filesep phaseFrequencyProfilesSubfolder], 'file')
    mkdir([mainFolder filesep phaseFrequencyProfilesSubfolder]);
  end
  
  % Load the population rate data
%   areaPhaseFOIindividualPR = load([rootFolder filesep MUAsFolder filesep 'area2pupilCohMats.mat']);
%   areaPhaseFOIindividualPR = areaPhaseFOIindividualPR.areaPhaseFOIindividual;
  
  options.figSize = figSize;
  options.mainFolder = mainFolder;
  options.phaseFrequencyProfilesSubfolder = phaseFrequencyProfilesSubfolder;
  options.compString = 'VsPupil';
  phaseFreqProfilePlotUpdate_units(areas, conditions, FOI, areaPhaseFOIindividual, [], options)
end
close all


%% PHASE HISTOGRAMS
if (drawPhaseHistos(1) || drawPhaseHistos(2)) && ~exist([mainFolder filesep histosSubfolder], 'file')
  mkdir([mainFolder filesep histosSubfolder]);
end
distributionStats = {};
phaseHistos = {};
for iCond = 1:numel(conditions)
  for iArea = 1:numel(areas)
    disp(['Processing unit data for ' conditions{iCond} ' ' areas{iArea}...
      ' (comparison # ' num2str((iCond-1)*numel(areas) + iArea) '/' num2str(numel(conditions)*numel(areas)) ')']);
    
    phase = recentrePhase(areaPhaseFOIindividual{iCond}{iArea},phaseCentre);
    if ~isempty(phase)
      histPhaseOverF = zeros(numel(edges),numel(FOI));
      for f = 1:numel(FOI)
        histPhaseOverF(:,f) = [sum(sum(isnan(phase(:,f))))...
          histcounts(phase(:,f), edges)];
      end
      if drawPhaseHistos(1) || drawPhaseHistos(2) || drawPhaseHistos(3)
        if drawPhaseHistos(2)
          plotCount = 0;
          sbPhase = figure('units', 'normalized', 'position', [0.002, .04, 1, .88], 'Visible', 'on');
        end
        for f = 1:numel(FOI)
          
          % Individual histos
          if drawPhaseHistos(1)
            phaseHistPlotUnadjusted(edges, histPhaseOverF(:,f), '# units', [mainFolder filesep histosSubfolder filesep...
              'PUPIL_' areas{iArea} 'VsPupil_' conditions{iCond} '_unit_phase_at_' num2str(FOI(f)) '_Hz'], 'on', opt.adjustAxes); % Figure
            [distributionStats.pRayleigh{iCond}{iArea}{f}, distributionStats.zRayleigh{iCond}{iArea}{f}, distributionStats.pOmnibus{iCond}{iArea}{f},...
              distributionStats.mOmnibus{iCond}{iArea}{f}, distributionStats.pRao{iCond}{iArea}{f}, distributionStats.U_Rao{iCond}{iArea}{f},...
              distributionStats.pV{iCond}{iArea}{f}, distributionStats.vV{iCond}{iArea}{f}, distributionStats.pHR{iCond}{iArea}{f},...
              distributionStats.T_HR{iCond}{iArea}{f}, distributionStats.nModes{iCond}{iArea}{f}, distributionStats.excessMass{iCond}{iArea}{f},...
              distributionStats.U2_KDE{iCond}{iArea}{f}, distributionStats.pKDE{iCond}{iArea}{f}, distributionStats.modes{iCond}{iArea}{f},...
              distributionStats.dipHDT{iCond}{iArea}{f}, distributionStats.pHDT{iCond}{iArea}{f}] = phaseDistributionTests(...
              phase(:,f), edges, true); % Stats
          end
          
          % Summary histos
          if drawPhaseHistos(2)
            if FOI(f) == 0.01 || FOI(f) == 0.02 || FOI(f) == 0.03 || FOI(f) == 0.1 || FOI(f) == 0.2 ||...
                FOI(f) == 1 || FOI(f) == 2 || FOI(f) == 4 || FOI(f) == 10
              figure(sbPhase);
              plotCount = plotCount + 1;
              subplot(3, 3, 10 - plotCount);
              if plotCount <= 2
                phaseHistSubplotUnadjusted(edges, histPhaseOverF(:,f), '', 'Phase (rad)', 'k', opt.adjustAxes);
                dispFOI(FOI(f), sum(histPhaseOverF(:,f)));
              elseif plotCount == 3
                phaseHistSubplotUnadjusted(edges, histPhaseOverF(:,f), '# units', 'Phase (rad)', 'k', opt.adjustAxes);
                dispFOI(FOI(f), sum(histPhaseOverF(:,f)));
              elseif plotCount <= 5
                phaseHistSubplotUnadjusted(edges, histPhaseOverF(:,f), '', '', 'k', opt.adjustAxes);
                dispFOI(FOI(f), sum(histPhaseOverF(:,f)));
              elseif plotCount == 6
                phaseHistSubplotUnadjusted(edges, histPhaseOverF(:,f), '# units', '', 'k', opt.adjustAxes);
                dispFOI(FOI(f), sum(histPhaseOverF(:,f)));
              elseif plotCount == 7
                phaseHistSubplotUnadjusted(edges, histPhaseOverF(:,f), '# units', '', 'k', opt.adjustAxes);
                dispFOI(FOI(f), sum(histPhaseOverF(:,f)));
                set(gcf, 'Color','w');
                print(sbPhase, [mainFolder filesep histosSubfolder filesep 'PUPIL_' areas{iArea} 'VsPupil_' conditions{iCond}...
                  '_unit_phase_at_p01_p02_p03_p1_p3_1_2_4_10_Hz.png'], '-dpng', '-r300');
                hgsave([mainFolder filesep histosSubfolder filesep 'PUPIL_' areas{iArea} 'VsPupil_' conditions{iCond}...
                  '_unit_phase_at_p01_p02_p03_p1_p3_1_2_4_10_Hz.fig'])
                close(sbPhase);
              end
            end
          end
        end
        
        % Phase frequency map
        if drawPhaseHistos(3)
          cm = 'cold';
          h = phaseGraph(['PUPIL Phase map ' areas{iArea} 'VsPupil ' conditions{iCond}], fliplr(FOI), xLim,...
            [0.01 0.1 1 10 100], edges, [edges(1) edges(end)], [-pi -pi/2 0 pi/2 pi]+phaseCentre,...
            fliplr(histPhaseOverF(2:end,:)), cm, 'on', opt.adjustAxes);
          set(gca, 'XTickLabel',{'0.01','0.1','1','10','100'});
          set(h,'color','white')
          if ~exist([mainFolder filesep mapsSubfolder], 'file')
            mkdir([mainFolder filesep mapsSubfolder]);
          end
          hgsave(h, [mainFolder filesep mapsSubfolder filesep 'PUPIL_Phase_map_' areas{iArea} 'VsPupil_' conditions{iCond}]);
          print(h, [mainFolder filesep mapsSubfolder filesep 'PUPIL_Phase_map_' areas{iArea} 'VsPupil_' conditions{iCond} '.png'],'-dpng','-r300');
          close(h);
        end
        
        phaseHistos{iCond}{iArea} = histPhaseOverF;
      end
    end
  end
end


%% STATS ON MEAN PHASE FREQUENCY PROFILES
% Mardia-Watson-Wheeler Uniform-Scores Test
fPEst = {};
fWTest = {};
strPMethod = {};
if numel(areaPhaseFOIindividual) > 1
  for iArea = 1:numel(areas)
    for f = 1:numel(FOI)
      if ~isempty(areaPhaseFOIindividual{1}{iArea}) && ~isempty(areaPhaseFOIindividual{2}{iArea})
        cvfX{1} = areaPhaseFOIindividual{1}{iArea}(:,f);
        cvfX{2} = areaPhaseFOIindividual{2}{iArea}(:,f);
        [~, fPEst{iArea}(f), fWTest{iArea}(f), strPMethod{iArea}{f}] = mardiatestn_circ_equal(cvfX);
      end
    end
  end
end

% Watson U2 test
pEst = {};
U2 = {};
pObs = {};
U2Obs = {};
if numel(areaPhaseFOIindividual) > 1
  for iArea = 1:numel(areas)
    for f = 1:numel(FOI)
      if ~isempty(areaPhaseFOIindividual{1}{iArea}) && ~isempty(areaPhaseFOIindividual{2}{iArea})
        a1 = areaPhaseFOIindividual{1}{iArea}(:,f);
        a2 = areaPhaseFOIindividual{2}{iArea}(:,f);
        [pEst{iArea}(f), U2{iArea}(f)] = watsons_U2_approx_p(a1, a2);
        [pObs{iArea}(f), U2Obs{iArea}(f)] = watsons_U2_perm_test(a1, a2, 1000);
      end
    end
  end
end


%% SAVE THE HISTOGRAM DATA
if ~exist([mainFolder filesep histosSubfolder], 'file')
  mkdir([mainFolder filesep histosSubfolder]);
end
if qualityCheck
  if drawPhaseHistos(1)
    save([mainFolder filesep histosSubfolder filesep 'phaseHistos_units_quality'], 'areaPhaseFOIindividual', 'phaseHistos', 'distributionStats',...
      'fPEst', 'fWTest', 'strPMethod', 'pEst', 'U2', 'pObs', 'U2Obs', 'areas', 'FOI', '-v7.3');
  else
    save([mainFolder filesep histosSubfolder filesep 'phaseHistos_units_quality'], 'areaPhaseFOIindividual', 'phaseHistos',...
      'fPEst', 'fWTest', 'strPMethod', 'pEst', 'U2', 'pObs', 'U2Obs', 'areas', 'FOI', '-v7.3');
  end
else
  if drawPhaseHistos(1)
    save([mainFolder filesep histosSubfolder filesep 'phaseHistos_units'], 'areaPhaseFOIindividual', 'phaseHistos', 'distributionStats',...
      'fPEst', 'fWTest', 'strPMethod', 'pEst', 'U2', 'pObs', 'U2Obs', 'areas', 'FOI', '-v7.3');
  else
    save([mainFolder filesep histosSubfolder filesep 'phaseHistos_units'], 'areaPhaseFOIindividual', 'phaseHistos',...
      'fPEst', 'fWTest', 'strPMethod', 'pEst', 'U2', 'pObs', 'U2Obs', 'areas', 'FOI', '-v7.3');
  end
end