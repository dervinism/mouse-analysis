function concatenatedData = concatenateSeries(dataDir, series, opt)
% Concatenates data series into a single data file. A helper function of
% globalFigs.
% Input: topDir - a folder with all your data stored.
%        analysisDir - a folder where loaded data is stored for each
%                      animal.
%        series - a cell array of string identifying series to concatenate.
%                 Each string takes the following form
%                 '<animal_id>_s<series>'.
%        qualitySeries - a cell array of string identifying cluster quality
%                        files corresponding to series data files.
% Output: concatenatedData structure.


% Initialise storage variables
concatenatedData.units = [];
concatenatedData.uCentres = [];
concatenatedData.uSeries = [];
concatenatedData.previousUnits = [];
concatenatedData.series = [];
concatenatedData.pr.beta = [];


% Loop over series
lastUnit = 0;
for s = 1:numel(series)
  disp(['processing series number ' num2str(s) '/' num2str(numel(series))]);
  
  % Load the file
  separator = find(series{s} == 's');
  animal = series{s}(1:separator-2);
  filename = [dataDir filesep animal filesep animal '.mat'];
  load(filename, 'dataStruct');
  
  % Load variables and check unit quality
  dbStruct = dataStruct.seriesData.(series{s});
  FOI = dbStruct.conf.FOI;
  fnsData = fieldnames(dbStruct.shankData);
  units = [];
  uCentres = [];
  uSeries = [];
  for sh = 1:numel(fnsData)
    if ~isempty(dbStruct.shankData.(['shank' num2str(sh)]).units)
      unitMetadata = dbStruct.shankData.(['shank' num2str(sh)]).unitMetadata;
      for u = torow(dbStruct.shankData.(['shank' num2str(sh)]).units)
        ind = find(u == dbStruct.shankData.(['shank' num2str(sh)]).units);
        if ~opt.qualityCheck || (unitMetadata(ind,7) >= cluDist && unitMetadata(ind,6) <= refractCont) % quality check
          units = [units; u]; %#ok<*AGROW>
          uCentres = [uCentres; unitMetadata(ind,3)];
          uSeries = [uSeries; s];
        end
      end
    end
  end
  units = sort(units);
  if isfield(dbStruct.popData, 'phaseCoh')
    pr.phaseCohFOI = dbStruct.popData.phaseCoh;
    pr.phaseCohHalvesFOI = dbStruct.popData.phaseCohHalves;
  else
    pr = [];
  end
  if (opt.pupilPlot || opt.mixedPlot) && isfield(dbStruct.popData, 'pupil')
    pupil.unitData = dbStruct.popData.pupil.unitData;
    pupil.popData = dbStruct.popData.pupil.popData;
  else
    pupil = [];
  end
  
  % Concatenate data
  if  s == 1
    concatenatedData.FOI = FOI;
  else
    assert(~sum(~(concatenatedData.FOI == FOI)));
  end
  concatenatedData.units = [concatenatedData.units; units];
  concatenatedData.uCentres = [concatenatedData.uCentres; uCentres];
  concatenatedData.uSeries = [concatenatedData.uSeries; uSeries];
  concatenatedData.previousUnits = [concatenatedData.previousUnits; ones(size(units))*lastUnit];
  concatenatedData.series = [concatenatedData.series; strings(size(units))+series{s}];
  unitCount = 0;
  for sh = 1:numel(fnsData)
    for u = 1:numel(units)
      iU = find(dbStruct.shankData.(['shank' num2str(sh)]).units == units(u));
      if isempty(iU)
        continue
      end
      if ~isempty(pr)
        if ~isempty(pr.phaseCohFOI.beta(unitCount+iU))
          concatenatedData.pr.beta = [concatenatedData.pr.beta; pr.phaseCohFOI.beta(unitCount+iU)];
        else
          concatenatedData.pr.beta = [concatenatedData.pr.beta; NaN];
        end
        concatenatedData.pr.phaseCohFOI{lastUnit+u}.phaseFOI = pr.phaseCohFOI.phaseFOI(unitCount+iU,:);
        concatenatedData.pr.phaseCohFOI{lastUnit+u}.cohFOI = pr.phaseCohFOI.cohFOI(unitCount+iU,:);
        concatenatedData.pr.phaseCohFOI{lastUnit+u}.coh_confFOI = pr.phaseCohFOI.coh_confFOI(unitCount+iU,:);
        actualFOI = dbStruct.shankData.(['shank' num2str(sh)]).phaseCoh{iU}.actualFOI;
        for iF = 1:numel(FOI)
          f = FOI(iF);
          if (sum(isnan(actualFOI)) == numel(actualFOI)) || (~isnan(f) && ((f >= 100 && abs(actualFOI(iF)-f) > 10) ||...
              (f < 100 && f >= 10 && abs(actualFOI(iF)-f) > 5) || (f < 10 && f >= 1 && abs(actualFOI(iF)-f) > 1) ||...
              (f < 1 && f >= 0.1 && abs(actualFOI(iF)-f) > 0.1) || (f < 0.1 && f >= 0.01 && abs(actualFOI(iF)-f) > 0.01)))
            concatenatedData.pr.phaseCohFOI{lastUnit+u}.phaseFOI(iF) = nan;
            concatenatedData.pr.phaseCohFOI{lastUnit+u}.cohFOI(iF) = nan;
            concatenatedData.pr.phaseCohFOI{lastUnit+u}.coh_confFOI(iF) = nan;
            actualFOI(iF) = nan;
          end
        end
        concatenatedData.pr.phaseCohFOI{lastUnit+u}.FOI = actualFOI;
        concatenatedData.pr.mfrHalves{lastUnit+u}.half1 = pr.phaseCohFOI.mfr_1sthalf(unitCount+iU,:);
        concatenatedData.pr.mfrHalves{lastUnit+u}.half2 = pr.phaseCohFOI.mfr_2ndhalf(unitCount+iU,:);
        if isempty(pr.phaseCohHalvesFOI) || iU > numel(pr.phaseCohHalvesFOI) || isempty(pr.phaseCohHalvesFOI{iU})
          uPhaseCohHalvesFOI.half1.phase = nan(1,numel(FOI));
          uPhaseCohHalvesFOI.half1.coh = nan(1,numel(FOI));
          uPhaseCohHalvesFOI.half1.coh_conf = nan(1,numel(FOI));
          uPhaseCohHalvesFOI.half1.FOI = nan(1,numel(FOI));
          uPhaseCohHalvesFOI.half2.phase = nan(1,numel(FOI));
          uPhaseCohHalvesFOI.half2.coh = nan(1,numel(FOI));
          uPhaseCohHalvesFOI.half2.coh_conf = nan(1,numel(FOI));
          uPhaseCohHalvesFOI.half2.FOI = nan(1,numel(FOI));
        else
          uPhaseCohHalvesFOI = pr.phaseCohHalvesFOI{iU};
          uPhaseCohHalvesFOI.half1 = rmfield(uPhaseCohHalvesFOI.half1,{'fInds','unit','series'});
          uPhaseCohHalvesFOI.half2 = rmfield(uPhaseCohHalvesFOI.half2,{'fInds','unit','series'});
        end
        concatenatedData.pr.phaseCohHalvesFOI{lastUnit+u} = uPhaseCohHalvesFOI;
        actualFOI = concatenatedData.pr.phaseCohHalvesFOI{lastUnit+u}.half1.FOI;
        for iF = 1:numel(FOI)
          f = FOI(iF);
          if (sum(isnan(actualFOI)) == numel(actualFOI)) || (~isnan(f) && ((f >= 100 && abs(actualFOI(iF)-f) > 10) ||...
              (f < 100 && f >= 10 && abs(actualFOI(iF)-f) > 5) || (f < 10 && f >= 1 && abs(actualFOI(iF)-f) > 1) ||...
              (f < 1 && f >= 0.1 && abs(actualFOI(iF)-f) > 0.1) || (f < 0.1 && f >= 0.01 && abs(actualFOI(iF)-f) > 0.01)))
            concatenatedData.pr.phaseCohFOI{lastUnit+u}.half1.phase(iF) = nan;
            concatenatedData.pr.phaseCohFOI{lastUnit+u}.half1.coh(iF) = nan;
            concatenatedData.pr.phaseCohFOI{lastUnit+u}.half1.coh_conf(iF) = nan;
            concatenatedData.pr.phaseCohFOI{lastUnit+u}.half2.phase(iF) = nan;
            concatenatedData.pr.phaseCohFOI{lastUnit+u}.half2.coh(iF) = nan;
            concatenatedData.pr.phaseCohFOI{lastUnit+u}.half2.coh_conf(iF) = nan;
            actualFOI(iF) = nan;
          end
        end
        concatenatedData.pr.phaseCohFOI{lastUnit+u}.half1.FOI = actualFOI;
        concatenatedData.pr.phaseCohFOI{lastUnit+u}.half2.FOI = actualFOI;
      end
      if (opt.pupilPlot || opt.mixedPlot) && ~isempty(pupil)
        concatenatedData.pupil.unitData.phaseCohFOI{lastUnit+u}.phaseFOI = pupil.unitData.phaseFOI(unitCount+iU,:);
        concatenatedData.pupil.unitData.phaseCohFOI{lastUnit+u}.cohFOI = pupil.unitData.cohFOI(unitCount+iU,:);
        concatenatedData.pupil.unitData.phaseCohFOI{lastUnit+u}.coh_confFOI = pupil.unitData.coh_confFOI(unitCount+iU,:);
        actualFOI = dbStruct.shankData.(['shank' num2str(sh)]).pupil.unitData{iU}.actualFOI;
        for iF = 1:numel(FOI)
          f = FOI(iF);
          if (sum(isnan(actualFOI)) == numel(actualFOI)) || (~isnan(f) && ((f >= 100 && abs(actualFOI(iF)-f) > 10) ||...
              (f < 100 && f >= 10 && abs(actualFOI(iF)-f) > 5) || (f < 10 && f >= 1 && abs(actualFOI(iF)-f) > 1) ||...
              (f < 1 && f >= 0.1 && abs(actualFOI(iF)-f) > 0.1) || (f < 0.1 && f >= 0.01 && abs(actualFOI(iF)-f) > 0.01)))
            concatenatedData.pupil.unitData.phaseCohFOI{lastUnit+u}.phaseFOI(iF) = nan;
            concatenatedData.pupil.unitData.phaseCohFOI{lastUnit+u}.cohFOI(iF) = nan;
            concatenatedData.pupil.unitData.phaseCohFOI{lastUnit+u}.coh_confFOI(iF) = nan;
            actualFOI(iF) = nan;
          end
        end
        concatenatedData.pupil.unitData.phaseCohFOI{lastUnit+u}.FOI = actualFOI;
        concatenatedData.pupil.unitData.phaseCohHalvesFOI{lastUnit+u}.half1.phase = pupil.unitData.phase_1sthalfFOI(unitCount+iU,:);
        concatenatedData.pupil.unitData.phaseCohHalvesFOI{lastUnit+u}.half1.coh = pupil.unitData.coh_1sthalfFOI(unitCount+iU,:);
        concatenatedData.pupil.unitData.phaseCohHalvesFOI{lastUnit+u}.half1.coh_conf = pupil.unitData.coh_1sthalf_confFOI(unitCount+iU,:);
        concatenatedData.pupil.unitData.phaseCohHalvesFOI{lastUnit+u}.half2.phase = pupil.unitData.phase_2ndhalfFOI(unitCount+iU,:);
        concatenatedData.pupil.unitData.phaseCohHalvesFOI{lastUnit+u}.half2.coh = pupil.unitData.coh_2ndhalfFOI(unitCount+iU,:);
        concatenatedData.pupil.unitData.phaseCohHalvesFOI{lastUnit+u}.half2.coh_conf = pupil.unitData.coh_2ndhalf_confFOI(unitCount+iU,:);
        if isfield(dbStruct.shankData.(['shank' num2str(sh)]).pupil.unitData{iU}, 'actualFOI_1sthalf')
          actualFOI = dbStruct.shankData.(['shank' num2str(sh)]).pupil.unitData{iU}.actualFOI_1sthalf;
        else
          actualFOI = concatenatedData.pr.phaseCohHalvesFOI{lastUnit+u}.half1.FOI;
        end
        for iF = 1:numel(FOI)
          f = FOI(iF);
          if (sum(isnan(actualFOI)) == numel(actualFOI)) || (~isnan(f) && ((f >= 100 && abs(actualFOI(iF)-f) > 10) ||...
              (f < 100 && f >= 10 && abs(actualFOI(iF)-f) > 5) || (f < 10 && f >= 1 && abs(actualFOI(iF)-f) > 1) ||...
              (f < 1 && f >= 0.1 && abs(actualFOI(iF)-f) > 0.1) || (f < 0.1 && f >= 0.01 && abs(actualFOI(iF)-f) > 0.01)))
            concatenatedData.pupil.unitData.phaseCohHalvesFOI{lastUnit+u}.half1.phase(iF) = nan;
            concatenatedData.pupil.unitData.phaseCohHalvesFOI{lastUnit+u}.half1.coh(iF) = nan;
            concatenatedData.pupil.unitData.phaseCohHalvesFOI{lastUnit+u}.half1.coh_conf(iF) = nan;
            concatenatedData.pupil.unitData.phaseCohHalvesFOI{lastUnit+u}.half2.phase(iF) = nan;
            concatenatedData.pupil.unitData.phaseCohHalvesFOI{lastUnit+u}.half2.coh(iF) = nan;
            concatenatedData.pupil.unitData.phaseCohHalvesFOI{lastUnit+u}.half2.coh_conf(iF) = nan;
            actualFOI(iF) = nan;
          end
        end
        concatenatedData.pupil.unitData.phaseCohHalvesFOI{lastUnit+u}.half1.FOI = actualFOI;
        concatenatedData.pupil.unitData.phaseCohHalvesFOI{lastUnit+u}.half2.FOI = actualFOI;
        if ~isempty(pr)
          concatenatedData.mixed{lastUnit+u}.spiking.phase = pr.phaseCohFOI.phaseFOI(unitCount+iU,:);
          concatenatedData.mixed{lastUnit+u}.spiking.coh = pr.phaseCohFOI.cohFOI(unitCount+iU,:);
          concatenatedData.mixed{lastUnit+u}.spiking.coh_conf = pr.phaseCohFOI.coh_confFOI(unitCount+iU,:);
          concatenatedData.mixed{lastUnit+u}.spiking.FOI = concatenatedData.pr.phaseCohFOI{lastUnit+u}.FOI;
        else
          concatenatedData.mixed{lastUnit+u}.spiking.phase = [];
          concatenatedData.mixed{lastUnit+u}.spiking.coh = [];
          concatenatedData.mixed{lastUnit+u}.spiking.coh_conf = [];
          concatenatedData.mixed{lastUnit+u}.spiking.FOI = [];
        end
        concatenatedData.mixed{lastUnit+u}.pupil.phase = pupil.unitData.phaseFOI(unitCount+iU,:);
        concatenatedData.mixed{lastUnit+u}.pupil.coh = pupil.unitData.cohFOI(unitCount+iU,:);
        concatenatedData.mixed{lastUnit+u}.pupil.coh_conf = pupil.unitData.coh_confFOI(unitCount+iU,:);
        concatenatedData.mixed{lastUnit+u}.pupil.FOI = concatenatedData.pupil.unitData.phaseCohFOI{lastUnit+u}.FOI;
      else
        concatenatedData.pupil.unitData.phaseCohFOI{lastUnit+u}.phaseFOI = [];
        concatenatedData.pupil.unitData.phaseCohFOI{lastUnit+u}.cohFOI = [];
        concatenatedData.pupil.unitData.phaseCohFOI{lastUnit+u}.coh_confFOI = [];
        concatenatedData.pupil.unitData.phaseCohHalvesFOI{lastUnit+u}.half1.phase = [];
        concatenatedData.pupil.unitData.phaseCohHalvesFOI{lastUnit+u}.half1.coh = [];
        concatenatedData.pupil.unitData.phaseCohHalvesFOI{lastUnit+u}.half1.coh_conf = [];
        concatenatedData.pupil.unitData.phaseCohHalvesFOI{lastUnit+u}.half1.FOI = [];
        concatenatedData.pupil.unitData.phaseCohHalvesFOI{lastUnit+u}.half2.phase = [];
        concatenatedData.pupil.unitData.phaseCohHalvesFOI{lastUnit+u}.half2.coh = [];
        concatenatedData.pupil.unitData.phaseCohHalvesFOI{lastUnit+u}.half2.coh_conf = [];
        concatenatedData.pupil.unitData.phaseCohHalvesFOI{lastUnit+u}.half2.FOI = [];
        if ~isempty(pr)
          concatenatedData.mixed{lastUnit+u}.spiking.phase = pr.phaseCohFOI.phaseFOI(unitCount+iU,:);
          concatenatedData.mixed{lastUnit+u}.spiking.coh = pr.phaseCohFOI.cohFOI(unitCount+iU,:);
          concatenatedData.mixed{lastUnit+u}.spiking.coh_conf = pr.phaseCohFOI.coh_confFOI(unitCount+iU,:);
          concatenatedData.mixed{lastUnit+u}.spiking.FOI = FOI;
        else
          concatenatedData.mixed{lastUnit+u}.spiking.phase = [];
          concatenatedData.mixed{lastUnit+u}.spiking.coh = [];
          concatenatedData.mixed{lastUnit+u}.spiking.coh_conf = [];
          concatenatedData.mixed{lastUnit+u}.spiking.FOI = [];
        end
        concatenatedData.mixed{lastUnit+u}.pupil.phase = [];
        concatenatedData.mixed{lastUnit+u}.pupil.coh = [];
        concatenatedData.mixed{lastUnit+u}.pupil.coh_conf = [];
        concatenatedData.mixed{lastUnit+u}.pupil.FOI = [];
      end
    end
    unitCount = unitCount + numel(dbStruct.shankData.(['shank' num2str(sh)]).units);
  end
  if (opt.pupilPlot || opt.mixedPlot) && ~isempty(pupil)
    concatenatedData.pupil.popData{s}.freq = pupil.popData.freq;
    concatenatedData.pupil.popData{s}.phase = pupil.popData.phase;
    concatenatedData.pupil.popData{s}.coh = pupil.popData.coh;
    concatenatedData.pupil.popData{s}.coh_conf = pupil.popData.coh_conf;
    concatenatedData.pupil.popData{s}.rateadjust_kappa = pupil.popData.rateadjust_kappa;
  else
    concatenatedData.pupil.popData{s}.freq = [];
    concatenatedData.pupil.popData{s}.phase = [];
    concatenatedData.pupil.popData{s}.coh = [];
    concatenatedData.pupil.popData{s}.coh_conf = [];
    concatenatedData.pupil.popData{s}.rateadjust_kappa = [];
  end
  lastUnit = numel(concatenatedData.units);
end