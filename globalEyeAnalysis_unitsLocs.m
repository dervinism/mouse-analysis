% Run this script to produce unit phase location profiles and histograms
% for area-to-pupil comparisons.
%
% The following data files are produced:
% dataDir\area2pupil\unitsFolder\phaseLocationProfilesSubfolder\globalUnitsLocs_area2pupil.mat or
%   dataDir\area2pupil\qualityUnitsFolder\phaseLocationProfilesSubfolder\globalUnitsLocs_area2pupil_quality.mat.
%   These files contain phase and coherence frequency profiles, unit
%   channel, probe location, and area location information. Correlation
%   analysis and statistical test information regarding phase location
%   distributions is also saved in these files.
%
% Phase location profile figures and summary subplots are saved either in
%   dataDir\area2pupil\unitsFolder\phaseLocationProfilesSubfolder or
%   dataDir\area2pupil\qualityUnitsFolder\phaseLocationProfilesSubfolder.

cleanUp


%% INITIALISE PARAMETERS
params
lists

animals = animalsOI;
areas = areasOI;

fullRun = true;
qualityCheck = false;
individualGraphs = true;
summaryGraphs = true;

if qualityCheck
  mainFolder = [dataDir filesep area2pupilDir filesep qualityUnitsFolder filesep phaseLocationProfilesSubfolder];
else
  mainFolder = [dataDir filesep area2pupilDir filesep unitsFolder filesep phaseLocationProfilesSubfolder];
end


%% COMPUTE VARIABLES NEEDED FOR DISPLAYING UNIT PHASE LOCATION PROFILES
if fullRun
  for animal = 1:numel(animals) % Loop through animals
    disp(['Loading data for animal ' animals{animal}]);
    load([dataDir filesep animals{animal} filesep animals{animal} '.mat'])
    fnsData = fieldnames(dataStruct.seriesData);
    
    % Initialise storage variables
    FOI = dataStruct.seriesData.(fnsData{1}).FOI;
    FOI(FOI > 2) = [];
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
    
    for dbCount = 1:numel(fnsData) % Loop through database entries
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
      
      % Test for exceptions
      if exceptionTest(except, seriesName)
        continue
      end
      
      % Determine if population rate > 0
      if firingRateTest(sum(dbStruct.popData.MUAsAll,1), dbStruct.conf.params.srData)
        continue
      end
      
      % Determine recording area
      area = determineArea(seriesName);
      
      % Determine recording condition (i.e., awake or anaesthesia)
      [breakClause, iCond] = series2condition(awake, anaesthesia, seriesName);
      if breakClause
        continue
      end
      
      % Report an error for multishank recordings and remove S1 units obtained with Neuronexus probes
      [probeConfFile, probeConfFileShort] = findProbeConfFile(dbStruct.io.dataDir);
      load(probeConfFile);
      probe = probeConfFileShort(8:end-4);
      [~, ~, areaName] = determineArea(seriesName);
      if numel(fieldnames(dbStruct.shankData)) > 1
        error('The analysis script is currently not set up to work with multi-shank recordings.')
      elseif ~strcmpi(probe, 'Neuropixels') && strcmpi(areaName, 'S1')
        continue
      end
      
      % Disqualify low quality units if needed
      units = [];
      for sh = 1:numel(fieldnames(dbStruct.shankData))
        units = [units; dbStruct.shankData.(['shank' num2str(sh)]).units]; %#ok<*AGROW>
      end
      if qualityCheck
        fullSeriesName = [animals{animal} '_s' seriesName];
        qualitySeries = [dataStruct.seriesData.(fullSeriesName).io.baseFilename '.qua.1.mat'];
        [~, qualityUnitInd] = qualityTest(dbStruct, qualitySeries, units, cluDist, refractCont);
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
          
          % Load and store phase of units
          phase = dbStruct.popData.pupil.unitData.phaseFOI(qualityUnitInd,:);
          if isempty(areaPhaseFOIindividual{iCondPlusAll}{iAreaPlusAll})
            areaPhaseFOIindividual{iCondPlusAll}{iAreaPlusAll} = phase;
          else
            areaPhaseFOIindividual{iCondPlusAll}{iAreaPlusAll} = [areaPhaseFOIindividual{iCondPlusAll}{iAreaPlusAll}; phase];
          end
          
          % Load and store coherence
          coh = dbStruct.popData.pupil.unitData.cohFOI(qualityUnitInd,:);
          if isempty(areaCohFOIindividual{iCondPlusAll}{iAreaPlusAll})
            areaCohFOIindividual{iCondPlusAll}{iAreaPlusAll} = coh;
          else
            areaCohFOIindividual{iCondPlusAll}{iAreaPlusAll} = [areaCohFOIindividual{iCondPlusAll}{iAreaPlusAll}; coh];
          end
          
          % Compute and store unit locations
          for sh = 1:numel(fieldnames(dbStruct.shankData))
            if isempty(areaProbeChFOIindividual{iCondPlusAll}{iAreaPlusAll})
                areaProbeChFOIindividual{iCondPlusAll}{iAreaPlusAll} = dbStruct.shankData.(['shank' num2str(sh)]).uCentres(qualityUnitInd);
                areaChFOIindividual{iCondPlusAll}{iAreaPlusAll} = areaProbeChFOIindividual{iCondPlusAll}{iAreaPlusAll} - min(dbStruct.dbSeries.chOI) + 1;
                areaProbeLocFOIindividual{iCondPlusAll}{iAreaPlusAll} = ycoords(areaProbeChFOIindividual{iCondPlusAll}{iAreaPlusAll})';
                areaLocFOIindividual{iCondPlusAll}{iAreaPlusAll} = ...
                  areaProbeLocFOIindividual{iCondPlusAll}{iAreaPlusAll} - ycoords(min(dbStruct.dbSeries.chOI))' + ycoords(1)/2;
              else
                uCentres = dbStruct.shankData.(['shank' num2str(sh)]).uCentres(qualityUnitInd);
                areaProbeChFOIindividual{iCondPlusAll}{iAreaPlusAll} = [areaProbeChFOIindividual{iCondPlusAll}{iAreaPlusAll}; uCentres]; %#ok<*AGROW>
                areaChFOIindividual{iCondPlusAll}{iAreaPlusAll} = [areaChFOIindividual{iCondPlusAll}{iAreaPlusAll}; uCentres - min(dbStruct.dbSeries.chOI) + 1];
                areaProbeLocFOIindividual{iCondPlusAll}{iAreaPlusAll} = [areaProbeLocFOIindividual{iCondPlusAll}{iAreaPlusAll}; ycoords(uCentres)'];
                areaLocFOIindividual{iCondPlusAll}{iAreaPlusAll} = [areaLocFOIindividual{iCondPlusAll}{iAreaPlusAll};...
                  ycoords(uCentres)' - ycoords(min(dbStruct.dbSeries.chOI))' + ycoords(1)/2];
            end
            for u = 1:numel(dbStruct.shankData.(['shank' num2str(sh)]).uCentres(qualityUnitInd))
              probes{iCondPlusAll}{iAreaPlusAll}{numel(probes{iCondPlusAll}{iAreaPlusAll})+1} = probe;
            end
          end
          assert(size(areaPhaseFOIindividual{iCondPlusAll}{iAreaPlusAll},1) == numel(areaProbeChFOIindividual{iCondPlusAll}{iAreaPlusAll}));
          assert(size(areaCohFOIindividual{iCondPlusAll}{iAreaPlusAll},1) == numel(areaProbeChFOIindividual{iCondPlusAll}{iAreaPlusAll}));
        end
      end
    end
  end
end

% Determine the file name and either save or load the data
if qualityCheck
  filename = [mainFolder filesep 'globalUnitsLocs_area2pupil_quality.mat']; %#ok<*UNRCH>
else
  filename = [mainFolder filesep 'globalUnitsLocs_area2pupil.mat'];
end
if fullRun
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
        individualGraphs, figPath, true, 'PUPIL', 'VsPupil');
      close all;
      
      % Summary graphs
      if summaryGraphs
        sbPhase = phaseLocationProfileSummary(areaPhaseFOIindividual{iCond}{iArea}, areaLocFOIindividual{iCond}{iArea},...
          FOI, rFOI{iCond}{iArea}, pvalFOI{iCond}{iArea}, nFOI{iCond}{iArea}, areas{iArea}, conditions{iCond}, figPath, 'PUPIL', 'VsPupil');
        close(sbPhase);
      end
    end
  end
end

% Save statistical test results
save(filename, 'rFOI','pvalFOI','nFOI','distributionStats', '-append');