% Run this script to perform analyses comparing LFP and MUA of different brain areas.


%% LOAD PRE-PROCESSED DATA
load(dataFile);


%% INITIALISE PARAMETERS
visibility = 'on';
overwrite = true;
intermediateSaving = false;


%% DETERMINE SERIES TO COMPARE
series2compare


%% CORRELATE LFP WITH CROSS_AREA POPULATION RATE
% LOOP THROUGH DB ENTRIES
figs = {};
animal = dataStruct.seriesData.(fnsData{1}).db(1).animal;
if ~isempty(dbEntries_c) && dbEntries_c(1) == inf
  dbEntries_c = 1:numel(series_c);
end
for dbCount = dbEntries_c
  for sca = 1:numel(series_ca{dbCount})
    
    % LOAD THE CONTENTS OF THE DB STRUCTURE VARIABLE
    dbStruct = dataStruct.seriesData.([animal '_s' series_c{dbCount}]);
    srData = dbStruct.conf.samplingParams.srData;
    entryName = [animal '_s' series_c{dbCount}];
    sca_name = [entryName '__' animal '_s' series_ca{dbCount}{sca}];
    ind = strfind(sca_name, '__');
    if ~strcmpi(entryName, sca_name(1:ind-1))
      continue
    end
    dbStruct_ca = dataStruct.seriesData.([animal '_s' series_ca{dbCount}{sca}]);
    
    % TEST IF THE ANALYSIS DATA ALREADY EXISTS
    if ~overwrite && isfield(dataStruct.seriesData_ca, dbStruct_ca) && isfield(dataStruct.seriesData_ca.(sca_name), 'lfpphaseCohDataPR')
      continue
    end
    
    strSep = strfind(sca_name,'__');
    sca_name1 = sca_name(1:strSep-1);
    strSep1 = strfind(sca_name1,'s');
    sca_name1 = sca_name1(strSep1+1:end);
    sca_name2 = sca_name(strSep+2:end);
    strSep2 = strfind(sca_name2,'s');
    sca_name2 = sca_name2(strSep2+1:end);
    figsubdirname = [sca_name1 'v' sca_name2];
    if ~exist(figsubdirname,'dir')
      mkdir(figsubdirname)
    end
    
    PR = sum(dbStruct.popData.MUAsAll,1);
    interpTimes = 1/srData:1/srData:numel(PR)/srData;
    if ~isempty(PR)
      [inds1, PR] = determineInds(dbStruct.dbSeries.period, srData, PR);
    end
    if isempty(PR)
      continue
    end
    
    % LOOP THROUGH LFP CHANNELS
    if isfield(dbStruct_ca, 'lfpPowerData')
      chOIDB = dbStruct_ca.lfpPowerData.chOIDB;
    else
      continue
    end
    for iCh = 1:numel(chOIDB)
      
      % LOAD LFP MEASURES
      lfpTimes = dbStruct_ca.lfpPowerData.lfpTimes;
      if iscell(dbStruct.dbSeries.period)
        assert(iscell(dbStruct_ca.dbSeries.period))
        for iPeriod = 1:numel(dbStruct.dbSeries.period)
          assert(dbStruct.dbSeries.period{iPeriod}(1) == dbStruct_ca.dbSeries.period{iPeriod}(1));
          assert(dbStruct.dbSeries.period{iPeriod}(2) == dbStruct_ca.dbSeries.period{iPeriod}(2));
        end
      else
        assert(dbStruct.dbSeries.period(1) == dbStruct_ca.dbSeries.period(1));
        assert(dbStruct.dbSeries.period(2) == dbStruct_ca.dbSeries.period(2));
      end
      [inds2, rippleRate] = determineInds(dbStruct_ca.dbSeries.period, 1/(lfpTimes(2)-lfpTimes(1)), dbStruct_ca.lfpPowerData.rippleRate{iCh});
      if isempty(rippleRate)
        continue
      end
      meanRippleRate = dbStruct_ca.lfpPowerData.meanRippleRate{iCh};
      theta2deltaRatio = dbStruct_ca.lfpPowerData.theta2deltaRatio{iCh}(inds2);
      slowPower = dbStruct_ca.lfpPowerData.slowPower{iCh}(inds2);
      fastPower = dbStruct_ca.lfpPowerData.fastPower{iCh}(inds2);
      ultraFastPower = dbStruct_ca.lfpPowerData.ultraFastPower{iCh}(inds2);
      
      % INTERPOLATE LFP MEASURES
      rippleRate = interp1(lfpTimes(inds2), rippleRate, interpTimes(inds1))';
      theta2deltaRatio = interp1(lfpTimes(inds2), theta2deltaRatio, interpTimes(inds1))';
      slowPower = interp1(lfpTimes(inds2), slowPower, interpTimes(inds1))';
      fastPower = interp1(lfpTimes(inds2), fastPower, interpTimes(inds1))';
      ultraFastPower = interp1(lfpTimes(inds2), ultraFastPower, interpTimes(inds1))';
      artefact = 10*srData;
      
      % CLEAN NAN's
      rippleRate(isnan(rippleRate)) = 0;
      theta2deltaRatio(isnan(theta2deltaRatio)) = 0;
      slowPower(isnan(slowPower)) = 0;
      fastPower(isnan(fastPower)) = 0;
      ultraFastPower(isnan(ultraFastPower)) = 0;
      
      % CORRELATION ANALYSES
      %     lfpPhaseCohDB.rippleRate{iCh}.corr = zeros(2,2);
      %     [lfpPhaseCohDB.rippleRate{iCh}.corr(1,1), lfpPhaseCohDB.rippleRate{iCh}.corr(1,2)] = corrSimple(PR(artefact:end-artefact), rippleRate(artefact:end-artefact), 'Pearson');
      %     [lfpPhaseCohDB.rippleRate{iCh}.corr(2,1), lfpPhaseCohDB.rippleRate{iCh}.corr(2,2)] = corrSimple(PR(artefact:end-artefact), rippleRate(artefact:end-artefact), 'Spearman');
      %     lfpPhaseCohDB.t2dRatio{iCh}.corr = zeros(2,2);
      %     [lfpPhaseCohDB.t2dRatio{iCh}.corr(1,1), lfpPhaseCohDB.t2dRatio{iCh}.corr(1,2)] = corrSimple(PR(artefact:end-artefact), theta2deltaRatio(artefact:end-artefact), 'Pearson');
      %     [lfpPhaseCohDB.t2dRatio{iCh}.corr(2,1), lfpPhaseCohDB.t2dRatio{iCh}.corr(2,2)] = corrSimple(PR(artefact:end-artefact), theta2deltaRatio(artefact:end-artefact), 'Spearman');
      %     lfpPhaseCohDB.slowPower{iCh}.corr = zeros(2,2);
      %     [lfpPhaseCohDB.slowPower{iCh}.corr(1,1), lfpPhaseCohDB.slowPower{iCh}.corr(1,2)] = corrSimple(PR(artefact:end-artefact), slowPower(artefact:end-artefact), 'Pearson');
      %     [lfpPhaseCohDB.slowPower{iCh}.corr(2,1), lfpPhaseCohDB.slowPower{iCh}.corr(2,2)] = corrSimple(PR(artefact:end-artefact), slowPower(artefact:end-artefact), 'Spearman');
      %     lfpPhaseCohDB.fastPower{iCh}.corr = zeros(2,2);
      %     [lfpPhaseCohDB.fastPower{iCh}.corr(1,1), lfpPhaseCohDB.fastPower{iCh}.corr(1,2)] = corrSimple(PR(artefact:end-artefact), fastPower(artefact:end-artefact), 'Pearson');
      %     [lfpPhaseCohDB.fastPower{iCh}.corr(2,1), lfpPhaseCohDB.fastPower{iCh}.corr(2,2)] = corrSimple(PR(artefact:end-artefact), fastPower(artefact:end-artefact), 'Spearman');
      %     lfpPhaseCohDB.ultraFastPower{iCh}.corr = zeros(2,2);
      %     [lfpPhaseCohDB.ultraFastPower{iCh}.corr(1,1), lfpPhaseCohDB.ultraFastPower{iCh}.corr(1,2)] = corrSimple(PR(artefact:end-artefact), ultraFastPower(artefact:end-artefact), 'Pearson');
      %     [lfpPhaseCohDB.ultraFastPower{iCh}.corr(2,1), lfpPhaseCohDB.ultraFastPower{iCh}.corr(2,2)] = corrSimple(PR(artefact:end-artefact), ultraFastPower(artefact:end-artefact), 'Spearman');
      
      % Plot graphs
      %figure; plot(PR,rippleRate, 'r.', 'MarkerSize',1)
      %figure; plot(PR,theta2deltaRatio, 'r.', 'MarkerSize',1)
      %figure; plot(PR,slowPower, 'r.', 'MarkerSize',1)
      %figure; plot(PR,fastPower, 'r.', 'MarkerSize',1)
      %figure; plot(PR,ultraFastPower, 'r.', 'MarkerSize',1)
%       figs{iCh} = figure('Visible', visibility); %#ok<*SAGROW>
%       plot(interpTimes,PR./mean(PR, 'omitnan'), 'LineWidth',2); hold on;
%       plot(interpTimes,rippleRate./mean(rippleRate, 'omitnan'))
%       plot(interpTimes,theta2deltaRatio./mean(theta2deltaRatio, 'omitnan'))
%       plot(interpTimes,slowPower./mean(slowPower, 'omitnan'))
%       plot(interpTimes,fastPower./mean(fastPower, 'omitnan'))
%       plot(interpTimes,ultraFastPower./mean(ultraFastPower, 'omitnan'))
%       legend('PR','Ripple rate','Theta2delta ratio','Slow power','Fast power','Ultra fast power')
%       xlabel('Time (s)')
%       ylabel('Normalised signal')
%       figName = [figsubdirname filesep entryName '_LFP_variousSignals_v_MUA_' 'ch' num2str(chOIDB(iCh))];
%       set(figs{iCh}, 'Name',figName);
%       hgsave(figs{iCh}, figName);
%       close(figs{iCh});
      
      % CALCULATE PHASE AND COHERENCE MEASURES
      optLFP = dbStruct.conf.optCoh;
      optLFP.winfactor = winfactor;
      optLFP.freqfactor = freqfactor;
      optLFP.monotoneFreq = true;
      optLFP.typespk1 = 'pb';
      optLFP.typespk2 = 'c';
      % Calculations
      [lfpPhaseCohDB.rippleRate{iCh}.freq, lfpPhaseCohDB.rippleRate{iCh}.coh, lfpPhaseCohDB.rippleRate{iCh}.phase,...
        lfpPhaseCohDB.rippleRate{iCh}.coh_conf, lfpPhaseCohDB.rippleRate{iCh}.phase_confU, lfpPhaseCohDB.rippleRate{iCh}.phase_confL,...
        lfpPhaseCohDB.rippleRate{iCh}.coh_halves_freq, lfpPhaseCohDB.rippleRate{iCh}.coh_halves, lfpPhaseCohDB.rippleRate{iCh}.coh_conf_halves,...
        lfpPhaseCohDB.rippleRate{iCh}.phase_halves, lfpPhaseCohDB.rippleRate{iCh}.phase_conf_halves, lfpPhaseCohDB.rippleRate{iCh}.stPRgr0,...
        lfpPhaseCohDB.rippleRate{iCh}.stPR, lfpPhaseCohDB.rippleRate{iCh}.m, lfpPhaseCohDB.rippleRate{iCh}.mv1, lfpPhaseCohDB.rippleRate{iCh}.mv2,...
        lfpPhaseCohDB.rippleRate{iCh}.mean_var, lfpPhaseCohDB.rippleRate{iCh}.mean_var_timeBin] = phaseCohLFP(PR(artefact:end-artefact),...
        rippleRate(artefact:end-artefact), srData, optLFP, 1e2);
      [lfpPhaseCohDB.t2dRatio{iCh}.freq, lfpPhaseCohDB.t2dRatio{iCh}.coh, lfpPhaseCohDB.t2dRatio{iCh}.phase,...
        lfpPhaseCohDB.t2dRatio{iCh}.coh_conf, lfpPhaseCohDB.t2dRatio{iCh}.phase_confU, lfpPhaseCohDB.t2dRatio{iCh}.phase_confL,...
        lfpPhaseCohDB.t2dRatio{iCh}.coh_halves_freq, lfpPhaseCohDB.t2dRatio{iCh}.coh_halves, lfpPhaseCohDB.t2dRatio{iCh}.coh_conf_halves,...
        lfpPhaseCohDB.t2dRatio{iCh}.phase_halves, lfpPhaseCohDB.t2dRatio{iCh}.phase_conf_halves, lfpPhaseCohDB.t2dRatio{iCh}.stPRgr0,...
        lfpPhaseCohDB.t2dRatio{iCh}.stPR, lfpPhaseCohDB.t2dRatio{iCh}.m, lfpPhaseCohDB.t2dRatio{iCh}.mv1, lfpPhaseCohDB.t2dRatio{iCh}.mv2,...
        lfpPhaseCohDB.t2dRatio{iCh}.mean_var, lfpPhaseCohDB.t2dRatio{iCh}.mean_var_timeBin] = phaseCohLFP(PR(artefact:end-artefact),...
        theta2deltaRatio(artefact:end-artefact), srData, optLFP, 1e2);
      [lfpPhaseCohDB.slowPower{iCh}.freq, lfpPhaseCohDB.slowPower{iCh}.coh, lfpPhaseCohDB.slowPower{iCh}.phase,...
        lfpPhaseCohDB.slowPower{iCh}.coh_conf, lfpPhaseCohDB.slowPower{iCh}.phase_confU, lfpPhaseCohDB.slowPower{iCh}.phase_confL,...
        lfpPhaseCohDB.slowPower{iCh}.coh_halves_freq, lfpPhaseCohDB.slowPower{iCh}.coh_halves, lfpPhaseCohDB.slowPower{iCh}.coh_conf_halves,...
        lfpPhaseCohDB.slowPower{iCh}.phase_halves, lfpPhaseCohDB.slowPower{iCh}.phase_conf_halves, lfpPhaseCohDB.slowPower{iCh}.stPRgr0,...
        lfpPhaseCohDB.slowPower{iCh}.stPR, lfpPhaseCohDB.slowPower{iCh}.m, lfpPhaseCohDB.slowPower{iCh}.mv1, lfpPhaseCohDB.slowPower{iCh}.mv2,...
        lfpPhaseCohDB.slowPower{iCh}.mean_var, lfpPhaseCohDB.slowPower{iCh}.mean_var_timeBin] = phaseCohLFP(PR(artefact:end-artefact),...
        slowPower(artefact:end-artefact), srData, optLFP, 1e2);
      [lfpPhaseCohDB.fastPower{iCh}.freq, lfpPhaseCohDB.fastPower{iCh}.coh, lfpPhaseCohDB.fastPower{iCh}.phase,...
        lfpPhaseCohDB.fastPower{iCh}.coh_conf, lfpPhaseCohDB.fastPower{iCh}.phase_confU, lfpPhaseCohDB.fastPower{iCh}.phase_confL,...
        lfpPhaseCohDB.fastPower{iCh}.coh_halves_freq, lfpPhaseCohDB.fastPower{iCh}.coh_halves, lfpPhaseCohDB.fastPower{iCh}.coh_conf_halves,...
        lfpPhaseCohDB.fastPower{iCh}.phase_halves, lfpPhaseCohDB.fastPower{iCh}.phase_conf_halves, lfpPhaseCohDB.fastPower{iCh}.stPRgr0,...
        lfpPhaseCohDB.fastPower{iCh}.stPR, lfpPhaseCohDB.fastPower{iCh}.m, lfpPhaseCohDB.fastPower{iCh}.mv1, lfpPhaseCohDB.fastPower{iCh}.mv2,...
        lfpPhaseCohDB.fastPower{iCh}.mean_var, lfpPhaseCohDB.fastPower{iCh}.mean_var_timeBin] = phaseCohLFP(PR(artefact:end-artefact),...
        fastPower(artefact:end-artefact), srData, optLFP, 1e2);
      [lfpPhaseCohDB.ultraFastPower{iCh}.freq, lfpPhaseCohDB.ultraFastPower{iCh}.coh, lfpPhaseCohDB.ultraFastPower{iCh}.phase,...
        lfpPhaseCohDB.ultraFastPower{iCh}.coh_conf, lfpPhaseCohDB.ultraFastPower{iCh}.phase_confU, lfpPhaseCohDB.ultraFastPower{iCh}.phase_confL,...
        lfpPhaseCohDB.ultraFastPower{iCh}.coh_halves_freq, lfpPhaseCohDB.ultraFastPower{iCh}.coh_halves, lfpPhaseCohDB.ultraFastPower{iCh}.coh_conf_halves,...
        lfpPhaseCohDB.ultraFastPower{iCh}.phase_halves, lfpPhaseCohDB.ultraFastPower{iCh}.phase_conf_halves, lfpPhaseCohDB.ultraFastPower{iCh}.stPRgr0,...
        lfpPhaseCohDB.ultraFastPower{iCh}.stPR, lfpPhaseCohDB.ultraFastPower{iCh}.m, lfpPhaseCohDB.ultraFastPower{iCh}.mv1, lfpPhaseCohDB.ultraFastPower{iCh}.mv2,...
        lfpPhaseCohDB.ultraFastPower{iCh}.mean_var, lfpPhaseCohDB.ultraFastPower{iCh}.mean_var_timeBin] = phaseCohLFP(PR(artefact:end-artefact),...
        ultraFastPower(artefact:end-artefact), srData, optLFP, 1e2);
      
      % Adjustments
      [mfr, mfr_1sthalf, mfr_2ndhalf] = rateCalc(PR(artefact:end-artefact), srData);
      dbStruct.conf.optCoh.winfactor = winfactor;
      dbStruct.conf.optCoh.freqfactor = freqfactor;
      dbStruct.conf.optCoh.monotoneFreq = true;
      [~, psd_halves, ~, psd] = psdCalc(PR(artefact:end-artefact), srData, dbStruct.conf.optCoh);
      [lfpPhaseCohDB.rippleRate{iCh}.rateadjust_kappa, lfpPhaseCohDB.rippleRate{iCh}.rateadjust_kappa_halves] = kappaCalc(...
        mfr, mfr_1sthalf, mfr_2ndhalf, psd, psd_halves, lfpPhaseCohDB.rippleRate{iCh}.coh, lfpPhaseCohDB.rippleRate{iCh}.coh_halves);
      [lfpPhaseCohDB.t2dRatio{iCh}.rateadjust_kappa, lfpPhaseCohDB.t2dRatio{iCh}.rateadjust_kappa_halves] = kappaCalc(...
        mfr, mfr_1sthalf, mfr_2ndhalf, psd, psd_halves, lfpPhaseCohDB.t2dRatio{iCh}.coh, lfpPhaseCohDB.t2dRatio{iCh}.coh_halves);
      [lfpPhaseCohDB.slowPower{iCh}.rateadjust_kappa, lfpPhaseCohDB.slowPower{iCh}.rateadjust_kappa_halves] = kappaCalc(...
        mfr, mfr_1sthalf, mfr_2ndhalf, psd, psd_halves, lfpPhaseCohDB.slowPower{iCh}.coh, lfpPhaseCohDB.slowPower{iCh}.coh_halves);
      [lfpPhaseCohDB.fastPower{iCh}.rateadjust_kappa, lfpPhaseCohDB.fastPower{iCh}.rateadjust_kappa_halves] = kappaCalc(...
        mfr, mfr_1sthalf, mfr_2ndhalf, psd, psd_halves, lfpPhaseCohDB.fastPower{iCh}.coh, lfpPhaseCohDB.fastPower{iCh}.coh_halves);
      [lfpPhaseCohDB.ultraFastPower{iCh}.rateadjust_kappa, lfpPhaseCohDB.ultraFastPower{iCh}.rateadjust_kappa_halves] = kappaCalc(...
        mfr, mfr_1sthalf, mfr_2ndhalf, psd, psd_halves, lfpPhaseCohDB.ultraFastPower{iCh}.coh, lfpPhaseCohDB.ultraFastPower{iCh}.coh_halves);
      
      % Graphs
      if sum(~isnan(lfpPhaseCohDB.rippleRate{iCh}.freq))
        figTmp = plotPhaseSpectrum(lfpPhaseCohDB.rippleRate{iCh}.freq', lfpPhaseCohDB.rippleRate{iCh}.coh, lfpPhaseCohDB.rippleRate{iCh}.phase,...
          lfpPhaseCohDB.rippleRate{iCh}.phase_confU, lfpPhaseCohDB.rippleRate{iCh}.phase_confL, lfpPhaseCohDB.rippleRate{iCh}.coh_conf, visibility);
        subplot(4,2,2) % Rate-adjusted coherence
        hold on, semilogx(lfpPhaseCohDB.rippleRate{iCh}.freq', lfpPhaseCohDB.rippleRate{iCh}.coh .* lfpPhaseCohDB.rippleRate{iCh}.rateadjust_kappa, 'b.-')
        hold on, semilogx(lfpPhaseCohDB.rippleRate{iCh}.freq',...
          (lfpPhaseCohDB.rippleRate{iCh}.coh + lfpPhaseCohDB.rippleRate{iCh}.coh_conf) .* lfpPhaseCohDB.rippleRate{iCh}.rateadjust_kappa, 'c--')
        hold on, semilogx(lfpPhaseCohDB.rippleRate{iCh}.freq',...
          (lfpPhaseCohDB.rippleRate{iCh}.coh - lfpPhaseCohDB.rippleRate{iCh}.coh_conf) .* lfpPhaseCohDB.rippleRate{iCh}.rateadjust_kappa, 'c--')
        figName = [figsubdirname filesep entryName '_LFP_rippleRate_v_MUA' '_PSDpolar_' 'ch' num2str(chOIDB(iCh))];
        set(figTmp, 'Name',figName);
        hgsave(figTmp, figName);
        close(figTmp);
        halfPhaseCohFig(lfpPhaseCohDB.rippleRate{iCh},...
          [figsubdirname filesep entryName '_LFP_rippleRate_v_MUA' '_phaseCohHalves_' 'ch' num2str(chOIDB(iCh))], 'Phase', 'Coherence', visibility);
      end
      
      if sum(~isnan(lfpPhaseCohDB.t2dRatio{iCh}.freq))
        figTmp = plotPhaseSpectrum(lfpPhaseCohDB.t2dRatio{iCh}.freq', lfpPhaseCohDB.t2dRatio{iCh}.coh, lfpPhaseCohDB.t2dRatio{iCh}.phase,...
          lfpPhaseCohDB.t2dRatio{iCh}.phase_confU, lfpPhaseCohDB.t2dRatio{iCh}.phase_confL, lfpPhaseCohDB.t2dRatio{iCh}.coh_conf, visibility);
        subplot(4,2,2) % Rate-adjusted coherence
        hold on, semilogx(lfpPhaseCohDB.t2dRatio{iCh}.freq', lfpPhaseCohDB.t2dRatio{iCh}.coh .* lfpPhaseCohDB.t2dRatio{iCh}.rateadjust_kappa, 'b.-')
        hold on, semilogx(lfpPhaseCohDB.t2dRatio{iCh}.freq',...
          (lfpPhaseCohDB.t2dRatio{iCh}.coh + lfpPhaseCohDB.t2dRatio{iCh}.coh_conf) .* lfpPhaseCohDB.t2dRatio{iCh}.rateadjust_kappa, 'c--')
        hold on, semilogx(lfpPhaseCohDB.t2dRatio{iCh}.freq',...
          (lfpPhaseCohDB.t2dRatio{iCh}.coh - lfpPhaseCohDB.t2dRatio{iCh}.coh_conf) .* lfpPhaseCohDB.t2dRatio{iCh}.rateadjust_kappa, 'c--')
        figName = [figsubdirname filesep entryName '_LFP_t2dRatio_v_MUA' '_PSDpolar_' 'ch' num2str(chOIDB(iCh))];
        set(figTmp, 'Name',figName);
        hgsave(figTmp, figName);
        close(figTmp);
        halfPhaseCohFig(lfpPhaseCohDB.t2dRatio{iCh},...
          [figsubdirname filesep entryName '_LFP_t2dRatio_v_MUA' '_phaseCohHalves_' 'ch' num2str(chOIDB(iCh))], 'Phase', 'Coherence', visibility);
      end
      
      if sum(~isnan(lfpPhaseCohDB.slowPower{iCh}.freq))
        figTmp = plotPhaseSpectrum(lfpPhaseCohDB.slowPower{iCh}.freq', lfpPhaseCohDB.slowPower{iCh}.coh, lfpPhaseCohDB.slowPower{iCh}.phase,...
          lfpPhaseCohDB.slowPower{iCh}.phase_confU, lfpPhaseCohDB.slowPower{iCh}.phase_confL, lfpPhaseCohDB.slowPower{iCh}.coh_conf, visibility);
        subplot(4,2,2) % Rate-adjusted coherence
        hold on, semilogx(lfpPhaseCohDB.slowPower{iCh}.freq', lfpPhaseCohDB.slowPower{iCh}.coh .* lfpPhaseCohDB.slowPower{iCh}.rateadjust_kappa, 'b.-')
        hold on, semilogx(lfpPhaseCohDB.slowPower{iCh}.freq',...
          (lfpPhaseCohDB.slowPower{iCh}.coh + lfpPhaseCohDB.slowPower{iCh}.coh_conf) .* lfpPhaseCohDB.slowPower{iCh}.rateadjust_kappa, 'c--')
        hold on, semilogx(lfpPhaseCohDB.slowPower{iCh}.freq',...
          (lfpPhaseCohDB.slowPower{iCh}.coh - lfpPhaseCohDB.slowPower{iCh}.coh_conf) .* lfpPhaseCohDB.slowPower{iCh}.rateadjust_kappa, 'c--')
        figName = [figsubdirname filesep entryName '_LFP_slowPower_v_MUA' '_PSDpolar_' 'ch' num2str(chOIDB(iCh))];
        set(figTmp, 'Name',figName);
        hgsave(figTmp, figName);
        close(figTmp);
        halfPhaseCohFig(lfpPhaseCohDB.slowPower{iCh},...
          [figsubdirname filesep entryName '_LFP_slowPower_v_MUA' '_phaseCohHalves_' 'ch' num2str(chOIDB(iCh))], 'Phase', 'Coherence', visibility);
      end
      
      if sum(~isnan(lfpPhaseCohDB.fastPower{iCh}.freq))
        figTmp = plotPhaseSpectrum(lfpPhaseCohDB.fastPower{iCh}.freq', lfpPhaseCohDB.fastPower{iCh}.coh, lfpPhaseCohDB.fastPower{iCh}.phase,...
          lfpPhaseCohDB.fastPower{iCh}.phase_confU, lfpPhaseCohDB.fastPower{iCh}.phase_confL, lfpPhaseCohDB.fastPower{iCh}.coh_conf, visibility);
        subplot(4,2,2) % Rate-adjusted coherence
        hold on, semilogx(lfpPhaseCohDB.fastPower{iCh}.freq', lfpPhaseCohDB.fastPower{iCh}.coh .* lfpPhaseCohDB.fastPower{iCh}.rateadjust_kappa, 'b.-')
        hold on, semilogx(lfpPhaseCohDB.fastPower{iCh}.freq',...
          (lfpPhaseCohDB.fastPower{iCh}.coh + lfpPhaseCohDB.fastPower{iCh}.coh_conf) .* lfpPhaseCohDB.fastPower{iCh}.rateadjust_kappa, 'c--')
        hold on, semilogx(lfpPhaseCohDB.fastPower{iCh}.freq',...
          (lfpPhaseCohDB.fastPower{iCh}.coh - lfpPhaseCohDB.fastPower{iCh}.coh_conf) .* lfpPhaseCohDB.fastPower{iCh}.rateadjust_kappa, 'c--')
        figName = [figsubdirname filesep entryName '_LFP_fastPower_v_MUA' '_PSDpolar_' 'ch' num2str(chOIDB(iCh))];
        set(figTmp, 'Name',figName);
        hgsave(figTmp, figName);
        close(figTmp);
        halfPhaseCohFig(lfpPhaseCohDB.fastPower{iCh},...
          [figsubdirname filesep entryName '_LFP_fastPower_v_MUA' '_phaseCohHalves_' 'ch' num2str(chOIDB(iCh))], 'Phase', 'Coherence', visibility);
      end
      
      if sum(~isnan(lfpPhaseCohDB.ultraFastPower{iCh}.freq))
        figTmp = plotPhaseSpectrum(lfpPhaseCohDB.ultraFastPower{iCh}.freq', lfpPhaseCohDB.ultraFastPower{iCh}.coh, lfpPhaseCohDB.ultraFastPower{iCh}.phase,...
          lfpPhaseCohDB.ultraFastPower{iCh}.phase_confU, lfpPhaseCohDB.ultraFastPower{iCh}.phase_confL, lfpPhaseCohDB.ultraFastPower{iCh}.coh_conf, visibility);
        subplot(4,2,2) % Rate-adjusted coherence
        hold on, semilogx(lfpPhaseCohDB.ultraFastPower{iCh}.freq', lfpPhaseCohDB.ultraFastPower{iCh}.coh .* lfpPhaseCohDB.ultraFastPower{iCh}.rateadjust_kappa, 'b.-')
        hold on, semilogx(lfpPhaseCohDB.ultraFastPower{iCh}.freq',...
          (lfpPhaseCohDB.ultraFastPower{iCh}.coh + lfpPhaseCohDB.ultraFastPower{iCh}.coh_conf) .* lfpPhaseCohDB.ultraFastPower{iCh}.rateadjust_kappa, 'c--')
        hold on, semilogx(lfpPhaseCohDB.ultraFastPower{iCh}.freq',...
          (lfpPhaseCohDB.ultraFastPower{iCh}.coh - lfpPhaseCohDB.ultraFastPower{iCh}.coh_conf) .* lfpPhaseCohDB.ultraFastPower{iCh}.rateadjust_kappa, 'c--')
        figName = [figsubdirname filesep entryName '_LFP_ultraFastPower_v_MUA' '_PSDpolar_' 'ch' num2str(chOIDB(iCh))];
        set(figTmp, 'Name',figName);
        hgsave(figTmp, figName);
        close(figTmp);
        halfPhaseCohFig(lfpPhaseCohDB.ultraFastPower{iCh},...
          [figsubdirname filesep entryName '_LFP_ultraFastPower_v_MUA' '_phaseCohHalves_' 'ch' num2str(chOIDB(iCh))], 'Phase', 'Coherence', visibility);
      end
    end
    
    % SAVE DATA
    dataString = ['dataStruct.seriesData_ca.' animal '_s' num2str(series_c{dbCount}) '__' animal '_s' series_ca{dbCount}{sca} '.lfpphaseCohDataPR = lfpPhaseCohDB;'];
    eval(dataString);
    if intermediateSaving
      save(dataFile,'dataStruct','-v7.3'); %#ok<*UNRCH>
    end
  end
end

if ~intermediateSaving
  save(dataFile,'dataStruct','-v7.3');
end
clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca



function [freq, coh, phase, coh_conf, phase_confU, phase_confL, coh_halves_freq, coh_halves, coh_conf_halves, phase_halves,...
  phase_conf_halves, stPRgr0, stPR, m, mv1, mv2, mean_var, mean_var_timeBin] = phaseCohLFP(PR, LFPsignal,...
  srData, optLFP, stPRwindow)

% PHASE AND COHERENCE ANALYSIS OF POPULATION SPIKING RELATIVE TO PUPIL AREA
[freq, coh, phase, coh_conf, phase_confU, phase_confL, coh_halves_freq, coh_halves, coh_conf_halves, phase_halves,...
  phase_conf_halves] = phaseCohCalc(PR, LFPsignal-mean(LFPsignal), srData, optLFP);

% CROSSCORRELATE POPULATION SPIKING AND PUPIL AREA
stPRgr0 = stprCalc(LFPsignal, PR, stPRwindow);
stPR = stprHalfCalc(LFPsignal, PR, stPRwindow);

% COMPARE MEAN VS VARIANCE USING DIFFERENT BIN SIZES
[m, mv1, mv2, mean_var, mean_var_timeBin] = meanVarCalc(LFPsignal, srData);
end

function halfPhaseCohFig(phaseCohData, figName, titleStr1, titleStr2, visibility)

% Half phase and coherence comparisons
figTmp = figure('units', 'normalized', 'position', [0.002, .04, 1/2, .88/1.04], 'Visible', visibility);
subplot(2, 1, 1);
phaseCohData.phase_halves(1,isnan(phaseCohData.phase_conf_halves(1,:))) = NaN;
phaseCohData.phase_halves(1,isnan(phaseCohData.phase_conf_halves(2,:))) = NaN;
phaseCohData.phase_halves(2,isnan(phaseCohData.phase_conf_halves(3,:))) = NaN;
phaseCohData.phase_halves(2,isnan(phaseCohData.phase_conf_halves(4,:))) = NaN;
semilogx(phaseCohData.coh_halves_freq, phaseCohData.phase_halves(1,:), 'r');
hold on, semilogx(phaseCohData.coh_halves_freq, phaseCohData.phase_halves(2,:), 'm');
hold on, semilogx(phaseCohData.coh_halves_freq, phaseCohData.phase_conf_halves(1,:), '--r');
hold on, semilogx(phaseCohData.coh_halves_freq, phaseCohData.phase_conf_halves(2,:), '--r');
hold on, semilogx(phaseCohData.coh_halves_freq, phaseCohData.phase_conf_halves(3,:), '--m');
hold on, semilogx(phaseCohData.coh_halves_freq, phaseCohData.phase_conf_halves(4,:), '--m');
title(titleStr1, 'Interpreter', 'none')
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');
legend('1st half', '2nd half');

subplot(2, 1, 2);
if isfield(phaseCohData, 'rateadjust_kappa_halves')
  factor = phaseCohData.rateadjust_kappa_halves;
else
  factor = ones(1,numel(phaseCohData.coh_halves(1,:)));
end
semilogx(phaseCohData.coh_halves_freq, phaseCohData.coh_halves(1,:) .* factor, 'r');
hold on, semilogx(phaseCohData.coh_halves_freq, phaseCohData.coh_halves(2,:), 'm');
hold on, semilogx(phaseCohData.coh_halves_freq, (phaseCohData.coh_halves(1,:) + phaseCohData.coh_conf_halves(1,:)) .* factor, '--r');
hold on, semilogx(phaseCohData.coh_halves_freq, (phaseCohData.coh_halves(1,:) - phaseCohData.coh_conf_halves(1,:)) .* factor, '--r');
hold on, semilogx(phaseCohData.coh_halves_freq, (phaseCohData.coh_halves(2,:) + phaseCohData.coh_conf_halves(2,:)) .* factor, '--m');
hold on, semilogx(phaseCohData.coh_halves_freq, (phaseCohData.coh_halves(2,:) - phaseCohData.coh_conf_halves(2,:)) .* factor, '--m');
title(titleStr2, 'Interpreter', 'none')
xlabel('Frequency (Hz)');
ylabel('Coherence');
legend('1st half', '2nd half');

set(figTmp, 'Name',figName);
hgsave(figTmp, figName);
close(figTmp);
end