% Run this script to produce unit phase location profiles and histograms
% for cross-area comparisons.
%
% The following data files are produced:
% dataDir\caDir\unitsFolder\phaseLocationProfilesSubfolder\globalUnitsLocs_ca.mat or
%   dataDir\caDir\unitsFolder\phaseLocationProfilesSubfolder\globalUnitsLocs_ca_reverse.mat or
%   dataDir\caDir\qualityUnitsFolder\phaseLocationProfilesSubfolder\globalUnitsLocs_ca_quality.mat or
%   dataDir\caDir\qualityUnitsFolder\phaseLocationProfilesSubfolder\globalUnitsLocs_ca_reverse_quality.mat.
%   These files contain phase and coherence frequency profiles, unit
%   channel, probe location, and area location information. Correlation
%   analysis and statistical test information regarding phase location
%   distributions is also saved in these files.
%
% Phase location profile figures and summary subplots are saved either in
%   dataDir\caDir\unitsFolder\phaseLocationProfilesSubfolder or
%   dataDir\caDir\qualityUnitsFolder\phaseLocationProfilesSubfolder.

cleanUp


%% INITIALISE PARAMETERS
params
lists

repository = 'uol';
if strcmp(repository,'all')
  rootFolder = [dataDir filesep caDir];
  animals = animalsOI;
elseif strcmp(repository,'uol')
  rootFolder = [dataDir filesep caDir_uol];
  animals = animalsUOLOI;
elseif strcmp(repository,'allensdk')
  rootFolder = [dataDir filesep caDir_allensdk];
  animals = animalsAllensdk;
end
areas = areas2compare;

fullRun = true;
reverse = false;
qualityCheck = false;
individualGraphs = true;
summaryGraphs = true;

if qualityCheck
  mainFolder = [rootFolder filesep qualityUnitsFolder filesep phaseLocationProfilesSubfolder];
else
  mainFolder = [rootFolder filesep unitsFolder filesep phaseLocationProfilesSubfolder];
end


%% COMPUTE VARIABLES NEEDED FOR DISPLAYING UNIT PHASE LOCATION PROFILES
if fullRun
  for animal = 1:numel(animals) % Loop through animals
    disp(['Loading data for animal ' animals{animal}]);
    load([dataDir filesep animals{animal} filesep animals{animal} '.mat'])
    fnsData_ca = fieldnames(dataStruct.seriesData_ca);
    
    % Initialise storage variables
    fnsData = fieldnames(dataStruct.seriesData);
    FOI = dataStruct.seriesData.(fnsData{1}).conf.FOI;
    if animal == 1
      areaCohFOIindividual = {};
      areaPhaseFOIindividual = {};
      areaProbeChFOIindividual = {};
      areaChFOIindividual = {};
      areaProbeLocFOIindividual = {};
      areaLocFOIindividual = {};
      probes = {};
      for iCond = 1:numel(conditions)
        areaCohFOIindividualCond = {};
        areaPhaseFOIindividualCond = {};
        areaProbeChFOIindividualCond = {};
        areaChFOIindividualCond = {};
        areaProbeLocFOIindividualCond = {};
        areaLocFOIindividualCond = {};
        probesCond = {};
        for iArea = 1:numel(areas)
          areaCohFOIindividualCond{iArea} = {}; %#ok<*SAGROW>
          areaPhaseFOIindividualCond{iArea} = {};
          areaProbeChFOIindividualCond{iArea} = {};
          areaChFOIindividualCond{iArea} = {};
          areaProbeLocFOIindividualCond{iArea} = {};
          areaLocFOIindividualCond{iArea} = {};
          probesCond{iArea} = {};
        end
        areaCohFOIindividual{iCond} = areaCohFOIindividualCond;
        areaPhaseFOIindividual{iCond} = areaPhaseFOIindividualCond;
        areaProbeChFOIindividual{iCond} = areaProbeChFOIindividualCond;
        areaChFOIindividual{iCond} = areaChFOIindividualCond;
        areaProbeLocFOIindividual{iCond} = areaProbeLocFOIindividualCond;
        areaLocFOIindividual{iCond} = areaLocFOIindividualCond;
        probes{iCond} = probesCond;
      end
    end
    
    for dbCount = 1:numel(fnsData_ca) % Loop through database entries
      dbStruct_ca = dataStruct.seriesData_ca.(fnsData_ca{dbCount});
      
      % Determine if series phase and coherence data exist
      if exist('seriesName1', 'var')
        prevRec = seriesName1(1:14);
      else
        prevRec = '';
      end
      [seriesName1, seriesName2] = seriesNames(fnsData_ca{dbCount});
      recording = seriesName1(1:14);
      if ~isfield(dbStruct_ca, 'popData') || ~isfield(dbStruct_ca.popData, 'phaseCoh') || isempty(dbStruct_ca.popData.phaseCoh)
        continue
      end
      
      % Test for exceptions
      if exceptionTest(except, seriesName1, seriesName2)
        continue
      end
      
      % Determine if population rate > 0
      [breakClause, spkDBmfr] = firingRateTest(sum(dataStruct.seriesData.([animals{animal} '_s' seriesName1]).popData.MUAsAll,1),...
        dataStruct.seriesData.([animals{animal} '_s' seriesName1]).conf.params.srData);
      if breakClause
        continue
      end
      [breakClause, PRmfr] = firingRateTest(sum(dataStruct.seriesData.([animals{animal} '_s' seriesName2]).popData.MUAsAll,1),...
        dataStruct.seriesData.([animals{animal} '_s' seriesName2]).conf.params.srData);
      if breakClause
        continue
      end
      
      % Determine area comparison and any grouped area comparisons
      [breakClause, comp, compNames, areasReverse] = series2Comparison(areas, seriesName1, seriesName2, reverse);
      if breakClause
        continue
      end
      
      % Determine recording condition (i.e., awake or anaesthesia)
      [breakClause, iCond] = series2condition(awake, anaesthesia, seriesName1, seriesName2);
      if breakClause
        continue
      end
      
      % Report an error for multishank recordings and remove S1 units obtained with Neuronexus probes
      if strcmp(dbStruct.conf.repository, 'allensdk')
        % do nothing
      else
        probeConfFile = findProbeConfFile(dbStruct.io.dataDir);
        load(probeConfFile);
      end
      probe = dbStruct.conf.probe;
      [~, ~, areaName] = determineArea(seriesName);
      if numel(fieldnames(dbStruct.shankData)) > 1
        error('The analysis script is currently not set up to work with multi-shank recordings.')
      elseif ~strcmpi(probe, 'Neuropixels') && strcmpi(areaName, 'lS1')
        continue
      end
      
      % Disqualify low quality units if needed
      units = [];
      for sh = 1:numel(dbStruct_ca.shankData)
        units = [units; dbStruct_ca.shankData.(['shank' num2str(sh)]).units]; %#ok<*AGROW>
      end
      if qualityCheck
        fullSeriesName1 = [animals{animal} '_s' seriesName1];
        qualitySeries = [dataStruct.seriesData.(fullSeriesName1).io.baseFilename '.qua.1.mat'];
        [~, qualityUnitInd] = qualityTest(dbStruct_ca, qualitySeries, units, cluDist, refractCont);
      else
        qualityUnitInd = 1:numel(units);
      end
      if isempty(qualityUnitInd)
        continue
      end
      units = units(qualityUnitInd);
      
      for iComp = 1:numel(comp) % Loop through non-grouped and grouped area comparisons
        area = comp(iComp);
        
        if numel(conditions) > 1
          condLoop = [iCond numel(conditions)];
        else
          condLoop = iCond;
        end
        for iCondPlusAll = condLoop % Loop through the main and pooled conditions
          
          % Load and store phase of units
          phase = dbStruct_ca.popData.phaseCoh.phaseFOI(qualityUnitInd,:);
          if isempty(areaPhaseFOIindividual{iCondPlusAll}{area})
            areaPhaseFOIindividual{iCondPlusAll}{area} = phase;
          else
            areaPhaseFOIindividual{iCondPlusAll}{area} = [areaPhaseFOIindividual{iCondPlusAll}{area}; phase];
          end
          
          % Load and store coherence of units
          coh = dbStruct_ca.popData.phaseCoh.cohFOI(qualityUnitInd,:);
          if isempty(areaCohFOIindividual{iCondPlusAll}{area})
            areaCohFOIindividual{iCondPlusAll}{area} = coh;
          else
            areaCohFOIindividual{iCondPlusAll}{area} = [areaCohFOIindividual{iCondPlusAll}{area}; coh];
          end
          
          % Compute and store unit locations
          if strcmp(dbStruct.conf.repository, 'allensdk')
            if isempty(areaProbeChFOIindividual{iCondPlusAll}{area})
              areaProbeChFOIindividual{iCondPlusAll}{area} = dbStruct.shankData.(['shank' num2str(sh)]).unitMetadata(qualityUnitInd,3); % Probe channel
              areaChFOIindividual{iCondPlusAll}{area} = areaProbeChFOIindividual{iCondPlusAll}{area} - min(dbStruct.dbSeries.chOI) + 1; % Relative probe channel
              areaProbeLocFOIindividual{iCondPlusAll}{area} = dbStruct.shankData.(['shank' num2str(sh)]).unitMetadata(qualityUnitInd,5); % Y-coordinate
              areaLocFOIindividual{iCondPlusAll}{area} = ...
                areaProbeLocFOIindividual{iCondPlusAll}{area} - ycoords(min(dbStruct.dbSeries.chOI))' + ycoords(1)/2; % Relative Y-coordinate (there is a problem here)
            else
              % do nothing
            end
          else
            for sh = 1:numel(fieldnames(dbStruct.shankData))
              if isempty(areaProbeChFOIindividual{iCondPlusAll}{area})
                areaProbeChFOIindividual{iCondPlusAll}{area} = dbStruct.shankData.(['shank' num2str(sh)]).unitMetadata(qualityUnitInd,3);
                areaChFOIindividual{iCondPlusAll}{area} = areaProbeChFOIindividual{iCondPlusAll}{area} - min(dbStruct.dbSeries.chOI) + 1;
                areaProbeLocFOIindividual{iCondPlusAll}{area} = ycoords(areaProbeChFOIindividual{iCondPlusAll}{area})';
                areaLocFOIindividual{iCondPlusAll}{area} = ...
                  areaProbeLocFOIindividual{iCondPlusAll}{area} - ycoords(min(dbStruct.dbSeries.chOI))' + ycoords(1)/2;
              else
                uCentres = dbStruct.shankData.(['shank' num2str(sh)]).unitMetadata(qualityUnitInd,3);
                areaProbeChFOIindividual{iCondPlusAll}{area} = [areaProbeChFOIindividual{iCondPlusAll}{area}; uCentres]; %#ok<*AGROW>
                areaChFOIindividual{iCondPlusAll}{area} = [areaChFOIindividual{iCondPlusAll}{area}; uCentres - min(dbStruct.dbSeries.chOI) + 1];
                areaProbeLocFOIindividual{iCondPlusAll}{area} = [areaProbeLocFOIindividual{iCondPlusAll}{area};...
                  dbStruct.shankData.(['shank' num2str(sh)]).unitMetadata(qualityUnitInd,5)];
                areaLocFOIindividual{iCondPlusAll}{area} = [areaLocFOIindividual{iCondPlusAll}{area};...
                  dbStruct.shankData.(['shank' num2str(sh)]).unitMetadata(qualityUnitInd,5) - ycoords(min(dbStruct.dbSeries.chOI))' + ycoords(1)/2];
              end
              for u = 1:numel(dbStruct.shankData.(['shank' num2str(sh)]).uCentres(qualityUnitInd))
                probes{iCondPlusAll}{area}{numel(probes{iCondPlusAll}{area})+1} = probe;
              end
            end
            assert(size(areaPhaseFOIindividual{iCondPlusAll}{area},1) == numel(areaProbeChFOIindividual{iCondPlusAll}{area}));
            assert(size(areaCohFOIindividual{iCondPlusAll}{area},1) == numel(areaProbeChFOIindividual{iCondPlusAll}{area}));
          end
        end
      end
    end
  end
  areas = areasReverse;
end

% Determine the file name and either save or load the data
if reverse
  if qualityCheck %#ok<*UNRCH>
    filename = [mainFolder filesep 'globalUnitsLocs_ca_reverse_quality.mat'];
  else
    filename = [mainFolder filesep 'globalUnitsLocs_ca_reverse.mat'];
  end
else
  if qualityCheck
    filename = [mainFolder filesep 'globalUnitsLocs_ca_quality.mat'];
  else
    filename = [mainFolder filesep 'globalUnitsLocs_ca.mat'];
  end
end
if fullRun
  if ~exist(mainFolder, 'file')
    mkdir(mainFolder);
  end
  save(filename, 'conditions','areas','FOI','areaPhaseFOIindividual','areaCohFOIindividual',...
    'areaProbeChFOIindividual','areaChFOIindividual','areaProbeLocFOIindividual','areaLocFOIindividual','probes');
else
  load(filename);
end


%% DISPLAY PHASE LOCATION PROFILES
figPath = fileparts(filename);
nSlices = [2 3 4 5 6];
rFOI = {};
pvalFOI = {};
nFOI = {};
distributionStats = {};
for iCond = 1:numel(conditions)
  for iArea = 1:numel(areas)
    if ~isempty(areaPhaseFOIindividual{iCond}{iArea})
      disp(['Processing unit location data for ' conditions{iCond} ' ' areas{iArea}...
        ' (comparison # ' num2str((iCond-1)*numel(areas) + iArea) '/' num2str(numel(conditions)*numel(areas)) ')']);
      
      % Phase location profile graphs for individual frequencies and location profile slice histograms
      [slicePhaseFig, rFOI{iCond}{iArea}, pvalFOI{iCond}{iArea}, nFOI{iCond}{iArea}, distributionStats{iCond}{iArea}] = phaseLocationProfile(...
        areaPhaseFOIindividual{iCond}{iArea}, edges, areaLocFOIindividual{iCond}{iArea}, FOI, nSlices, areas{iArea}, conditions{iCond},...
        individualGraphs, figPath, true);
      close all;
      
      % Summary graphs
      if summaryGraphs
        sbPhase = phaseLocationProfileSummary(areaPhaseFOIindividual{iCond}{iArea}, areaLocFOIindividual{iCond}{iArea},...
          FOI, rFOI{iCond}{iArea}, pvalFOI{iCond}{iArea}, nFOI{iCond}{iArea}, areas{iArea}, conditions{iCond}, figPath);
        close(sbPhase);
      end
    end
  end
end

% Save statistical test results
save(filename, 'rFOI','pvalFOI','nFOI','distributionStats', '-append');