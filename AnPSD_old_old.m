% Run this script to perform coherence and phase analyses on mouse data (even older; definitely don't use).

% close all
% clear all %#ok<CLALL>
% clc

% CREATE DB AND INITIALISE BASIC PARAMETERS
if ispc
  % addpath \\basket.cortexlab.net\data\mush\circStat
  addpath R:\CSN\Shared\Dynamics\Code\circStat
else
  addpath /data/mush/circStat
end

fileList = dir('spk*.mat');
if ~isempty(fileList)
  shankDataFile = fileList(end);
  load(shankDataFile.name,'db','figsubdirname','UPF','q','x_lim','params','opt','shankIDs');
end

if ~exist('db', 'var')
  error('Please run some makedb script');
end

if ~exist('params', 'var')
  params.Fs = 400;
end
if ~exist('opt', 'var')
  opt.maxFreq = 300; %110;
  opt.winfactor = 10;
  opt.freqfactor = 1.333;
  opt.tapers = 5;
end

if ~exist('figsubdirname', 'var')
  figsubdirname = 'AnPSD';
end
if ~exist('UPF', 'var')
  UPF = 10; % Units Per Figure
end
if ~exist('x_lim', 'var')
  x_lim = [0.02 50];
end

if ~exist('q', 'var')
  q = []; % for saving results
end

% LOOP THROUGH DB ENTRIES
for dbCount = length(db):-1:1
  if isempty(fileList)
    baseFilename = [topDir filesep db(dbCount).animal filesep num2str(db(dbCount).series) filesep ...
      db(dbCount).basefilename];
    spkAllShanks = {};
    
% SPIKES ON ALL SHANKS
    shankIDs = cell2mat(db(dbCount).shankGroups);
    shankMUAs = loadAsMUA(baseFilename, shankIDs, SR, 1/params.Fs);
    
    % There will be more than one shank group if there was more than one
    % probe used in the recording (or we decided to subdivide the shanks/tetrodes for
    % any other reason).
% SPIKES ON ALL SHANKS GROUPED BY A SHANK GROUP
    for i = 1:length(db(dbCount).shankGroups)
      spkAllShanks{i} = [];
      spkAllShanks{i} = zeros(numel(db(dbCount).shankGroups{i}), size(shankMUAs, 2));
      for shCount = 1:numel(db(dbCount).shankGroups{i})
        spkAllShanks{i}(shCount, :) = shankMUAs(cell2mat(db(dbCount).shankGroups) == db(dbCount).shankGroups{i}(shCount), :);
      end
      if ~isempty(db(dbCount).evtFilename) % cut to duration of the eventfile
        error('no longer supported?')
        evt = load([topDir filesep db(dbCount).animal filesep num2str(db(dbCount).series) filesep ...
          db(dbCount).evtFilename]); %#ok<*UNRCH>
        
        evt = evt(:, 1) * params.Fs/1000; % convert to the time resolution we use
        spkAllShanks{i} = spkAllShanks{i}(:, ceil(evt(1)):floor(evt(end)));
      end
    end
  end

% LOOP THROUGH SHANKS
  for sh = shankIDs
    fprintf('%s shank %d -------------------\n', db(dbCount).entryName, sh);
    
    for i = 1:length(db(dbCount).shankGroups) % to which group of shanks the current one belongs
      if any(db(dbCount).shankGroups{i} == sh)
        currentGrp  = i;
      end
    end
    
    shankDataFile = ['spk' num2str(dbCount) '_' num2str(sh)];
    if exist([shankDataFile '.mat'],'file')
      load([shankDataFile '.mat']);
    else
      
% SPIKES FOR SHANK OF INTEREST
      spkAll = shankMUAs(shankIDs(sh), :);  % population rate on this shank/tetrode

% TRUNCATE TO THE EVENT PERIOD
      if ~isempty(db(dbCount).evtFilename) % cut to duration of the eventfile
        error('no longer supported?')
        evt = load([topDir filesep db(dbCount).animal filesep num2str(db(dbCount).series) filesep ...
          db(dbCount).evtFilename]);
        evt = evt(:, 1) * params.Fs/1000; % convert to the time resolution we use
        raster = raster(:, ceil(evt(1)):floor(evt(end)));
        clear evt
      elseif size(spkAll, 2) < size(spkAllShanks{currentGrp}, 2)
        spkAll(size(spkAllShanks{currentGrp}, 2)) = 0; % make sure it's at least size(spkAllShanks{currentGrp}, 2) long, i.e., expand it to be the same duration...
      end

% SPIKE CLUSTER (UNIT) IDs FOR SHANK OF INTEREST
      clu = load([baseFilename '.clu.' num2str(sh)]); clu = clu(2:end); % first entry is the total number of units (remove it)

% UNITS FOR SHANK OF INTEREST
      units = unique(clu); 
      units = units(units > 1); % remove noise and MUA spikes

% REMOVE UNITS WITH FEWER THAN 100 SPIKES
      badUnits = [];
      for u = 1:numel(units)
        %       [uQ, acg_violation] = getUnitQuality(db(dbCount).animal, db(dbCount).series, sh, units(u));
        %       if isempty(uQ)
        %         warning('no unit quality data?!')
        %       elseif acg_violation > 0.004 || uQ < 15 || (uQ < 20 && acg_violation > 0.002)
        %         badUnits(end+1) = units(u);
        %       end

        if sum(clu == units(u)) < 100 % do not even have 100 spikes
          badUnits(end+1) = units(u); %#ok<*SAGROW>
        end
        if isempty(strfind(dbID, 'acute')) %#ok<*STREMP>
          Tstart = 1; % we presume it's a chronic recording ==> assume stationary throughout
          Tend = numel(spkAll);
        else
          error('not handled right now')
          % [Tstart, Tend] = findStationaryInterval(T, SR);
          %         if isnan(Tstart) || (Tend-Tstart)/SR < 2400
          %           badUnits(end+1) = units(u); % Don't have 40min. of reasonably stationary activity (with reasonable firing rate)
          %         end
        end
      end
      [units, ix] = setdiff(units, badUnits);    
      clear badUnits uQ acg_violation ix clu i shCount Tstart Tend u

% CREATE RASTER MATRIX
      load_opt.selectedUnits = units;
      spk = loadAsRasterSparse(baseFilename, sh, SR, 1/params.Fs, load_opt);
      save([shankDataFile '.mat'],'spk','db','SR','dbID','metadb', 'units',...
        'opt','params','topDir','UPF','x_lim','shankMUAs','spkAll',...
        'spkAllShanks','q','figsubdirname','shankIDs','baseFilename',...
        'currentGrp','load_opt');
    end
    
% LOOP THROUGH UNITS
    figH = figure('units', 'normalized', 'position', [0.002, .04, 1, .88]);
    for u = 1:numel(units)
      clear freq1 P1 freq2 P2 freq1112 P1112 P1112c C2 phi1 phi2 confC2 C11 C12 C1 phi1 C1c phi1u phi1l q
      
% IDENTIFY FIGURE AND SUBFIGURE
      if numel(units) < UPF
        whichSubplot = u;
      else
        whichSubplot = mod(u, UPF); if ~whichSubplot; whichSubplot=UPF; end
      end
      
% CREATE SUMMARY DATA STRUCTURE
      q = cell2struct({db(dbCount).animal, db(dbCount).series, db(dbCount).entryName, sh, units(u)}, {'animal', 'series', 'dbentryName', 'shank', 'unit'}, 2);
      
      if isempty(strfind(dbID, 'acute'))
        Tstart = 1;
        Tend = numel(spkAll);
        % Tstart = T(1); % we presume it's a chronic recording ==> assume stationary throughout
        % Tend = T(end);
      else
        error('not handled anymore')
        % [Tstart, Tend] = findStationaryInterval(T, SR);
      end

% SPIKES FOR UNIT OF INTEREST
      spkOI = full(spk(u,:));
      
% SPIKES FOR THE VICINITY OF UNIT OF INTEREST (in this case the whole shank)
      if isempty(strfind(dbID, 'IMEC'))
        localMUA = []; % this is not an IMEC probe
      else
%         localMUA = loadMUAaroundcontact([topDir filesep db(dbCount).animal filesep num2str(db(dbCount).series)], units(u), SR, 1/params.Fs);
%         assert(all(spk <= localMUA) && all(localMUA <= spkAll))
        localMUA = spkAll;
      end
      assert(all(spk <= spkAll))

% SUMMARY DATA
      % mean firing rate (mfr)
      q.mfr = mean(spk(Tstart:Tend))*params.Fs;
      q.mfr_1sthalf = mean(spk(Tstart:round((Tend+Tstart)/2)))*params.Fs;
      q.mfr_2ndhalf = mean(spk(1+round((Tend+Tstart)/2):Tend))*params.Fs;
      
      % PSD
      [freq_1sthalf, psd_1sthalf] = freqDependentWindowCoherence(spk(Tstart:round((Tend+Tstart)/2)), [], 1/params.Fs, [], opt); % psd
      [freq_2ndhalf, psd_2ndhalf] = freqDependentWindowCoherence(spk(1+round((Tend+Tstart)/2):Tend), [], 1/params.Fs, [], opt); % psd
      [freq, psd, ~, psd_conf] = freqDependentWindowCoherence(spk(Tstart:Tend), [], 1/params.Fs, [], opt); % psd
% PLOT PSD
      subplot(5, min(UPF,numel(units)), whichSubplot);
      hold on, plot(freq_1sthalf, psd_1sthalf, 'r');
      hold on, plot(freq_2ndhalf, psd_2ndhalf, 'm'); %hold on, plot(freq2, P2/P2(end), 'm');
      
      q.psd_halves = [psd_1sthalf; psd_2ndhalf];
      q.psd_halves_freq = freq_1sthalf; assert(numel(freq_1sthalf) == numel(freq_2ndhalf) && max(abs(freq_1sthalf - freq_2ndhalf)) < 1e-9)
      q.psd_freq = freq; q.psd = psd; q.psd_conf = psd_conf; q.psd_numelSignal = Tend - Tstart + 1;
      clear psd psd_conf freq freq_1sthalf freq_2ndhalf
      
      % Coherence with MUA
      subplot(5, min(UPF,numel(units)), whichSubplot + min(UPF,numel(units)));
      if length(spkAllShanks) == 1 && length(db(dbCount).shank) == 1 && ~isempty(localMUA) % all we've got is one shank, which we handle only for IMEC recordings
        [freq_1sthalf, coh_1sthalf, phase_1sthalf] = ...
          freqDependentWindowCoherence(spkAll(Tstart:round((Tend+Tstart)/2))'-spk(Tstart:round((Tend+Tstart)/2))', ...
          spk(Tstart:round((Tend+Tstart)/2))', 1/params.Fs, [], opt);
        [freq_2ndhalf, coh_2ndhalf, phase_2ndhalf] = ...
          freqDependentWindowCoherence(spkAll(1+round((Tend+Tstart)/2):Tend)'-spk(1+round((Tend+Tstart)/2):Tend)', ...
          spk(1+round((Tend+Tstart)/2):Tend)', 1/params.Fs, [], opt);
% PLOT COHERENCE
        hold on, plot(freq_1sthalf, coh_1sthalf, '.-r');
        hold on, plot(freq_2ndhalf, coh_2ndhalf, '.-m');
% PLOT PHASE
        subplot(5, min(UPF,numel(units)), whichSubplot + 2*min(UPF,numel(units)));
        hold on, plot(freq_1sthalf, unwrap(phase_1sthalf), '.-r');
        hold on, plot(freq_2ndhalf, unwrap(phase_2ndhalf), '.-m');
        
        [freq, coh, phase, coh_conf, phase_confU, phase_confL] = ...
              freqDependentWindowCoherence(spkAll(Tstart:Tend)'-spk(Tstart:Tend)', ...
              spk(Tstart:Tend)', 1/params.Fs, [], opt);
        q = catstruct(q, ...
              cell2struct({freq, coh, phase, coh_conf, phase_confU, phase_confL}, ...
              {'freq','coh', 'phase', 'coh_conf', 'phase_confU', 'phase_confL'}, 2));            
        q.coh_halves_freq = freq_1sthalf; assert(numel(freq_1sthalf) == numel(freq_2ndhalf) && max(abs(freq_1sthalf - freq_2ndhalf)) < 1e-9)
        q.coh_halves = [torow(coh_1sthalf); torow(coh_2ndhalf)];
        q.phase_halves = [torow(phase_1sthalf); torow(phase_2ndhalf)];    
            
% PLOT CIRCULAR PHASE
        figTmp = plotPhaseSpectrum(freq', coh, phase, phase_confU, phase_confL, coh_conf);
        % Compute the rate adjustment factors for coherence report, per Aoi et al.
        % we adjust as if the neuron's rate is 1 spk/s
        q.rateadjust_kappa = (1 + (q.mfr - 1)*q.mfr ./ q.psd);
        q.rateadjust_kappa(q.rateadjust_kappa < 0) = NaN; % the values are not supposed to be negative, it sometimes happens for neurons with low firing rate,
        % presumably because of the effect refractory period has on q.psd
        q.rateadjust_kappa = q.rateadjust_kappa.^-0.5;
        assert(isreal(q.psd) && isreal(coh) && isreal(q.rateadjust_kappa))

% PLOT RATE-ADJUSTED COHERENCE
        subplot(4,2,2)
        hold on, semilogx(q.freq, coh .* q.rateadjust_kappa, 'b.-')
        hold on, semilogx(q.freq, (coh + coh_conf) .* q.rateadjust_kappa, 'c--')
        hold on, semilogx(q.freq, max(0,(coh - coh_conf) .* q.rateadjust_kappa), 'c--')
        ylabel('coh, rate-adjusted coh')
        title([db(dbCount).entryName '_' num2str(sh) '_u' num2str(units(u))], ...
          'Interpreter', 'none')
%f        hgsave(figTmp, [figsubdirname filesep db(dbCount).entryName '_' num2str(sh) ...
%f          '_PSDpolar_u' num2str(units(u))]);
        close(figTmp);
        
        q.mfr_popRate = mean(spkAll(Tstart:Tend))*params.Fs;
        q.mfr_localMUARate = mean(localMUA(Tstart:Tend))*params.Fs;
      else % we have more than 1 shank
        [freq2, C2, phi2, confC2] = freqDependentWindowCoherence(spkAll(Tstart:Tend)' - spk(Tstart:Tend)', ...
          spk(Tstart:Tend)', 1/params.Fs, [], opt);
        hold on, plot(freq2, C2, '.-b');
        hold on, fill([freq2; freq2(end:-1:1)], [C2+confC2; C2(end:-1:1)-confC2(end:-1:1)], 'b', 'FaceAlpha', 0.25, 'LineStyle', 'none');
        for grp = 1:length(spkAllShanks)
          if grp == currentGrp
            
            q.mfr_popRate = mean(sum(spkAllShanks{grp}(:, Tstart:Tend), 1)-spkAll(Tstart:Tend))*params.Fs; %f GROUP OF SHANKS MINUS LOCAL SHANK
            q.mfr_localPopRate = mean(spkAll(Tstart:Tend) - spk(Tstart:Tend))*params.Fs; %f LOCAL SHANK MINUS UNIT
            assert(~any(spkAll > sum(spkAllShanks{grp}, 1)))
            
            % compute coherency with the local (i.e. on the same shank/tetrode) population rate. If spike sorting excludes the
            % possibility of more than 1 spike in (say) 1ms detection window, a spurious correlation might be introduced on fast timescale
            
            [~, coh, phase, coh_conf, phase_confU, phase_confL] = ... we skip the first output (freq) because it's the same as that below
              freqDependentWindowCoherence(spkAll(Tstart:Tend)' - spk(Tstart:Tend)', spk(Tstart:Tend)', 1/params.Fs, [], opt);
            
            q = catstruct(q, ...
              cell2struct({coh, phase, coh_conf, phase_confU, phase_confL}, ...
              {'coh_local', 'phase_local', 'coh_conf_local', 'phase_confU_local', 'phase_confL_local'}, 2));
            
            % compute coherency where population rate excludes the local shank/tetrode
            [freq_1sthalf, coh_1sthalf, phase_1sthalf] = ...
              freqDependentWindowCoherence(sum(spkAllShanks{grp}(:, Tstart:round((Tend+Tstart)/2)), 1)'-spkAll(Tstart:round((Tend+Tstart)/2))', ...
              spk(Tstart:round((Tend+Tstart)/2))', 1/params.Fs, [], opt);
            [freq_2ndhalf, coh_2ndhalf, phase_2ndhalf] = ...
              freqDependentWindowCoherence(sum(spkAllShanks{grp}(:, 1+round((Tend+Tstart)/2):Tend), 1)'-spkAll(1+round((Tend+Tstart)/2):Tend)', ...
              spk(1+round((Tend+Tstart)/2):Tend)', 1/params.Fs, [], opt);
            
            [freq, coh, phase, coh_conf, phase_confU, phase_confL] = ...
              freqDependentWindowCoherence(sum(spkAllShanks{grp}(:, Tstart:Tend), 1)'-spkAll(Tstart:Tend)', ...
              spk(Tstart:Tend)', 1/params.Fs, [], opt);
            [freq_s, coh_s, phase_s, coh_conf_s, phase_confU_s, phase_confL_s] = ... shuffled (spk is taken from end to start)
              freqDependentWindowCoherence(sum(spkAllShanks{grp}(:, Tstart:Tend), 1)'-spkAll(Tstart:Tend)', ...
              spk(Tend:-1:Tstart)', 1/params.Fs, [], opt);
            
            %             assert(max(abs(freq1112(numel(freq1112)-numel(freq11)+1:end) - freq11)) == 0)
            %             freq1112 = freq1112(numel(freq1112)-numel(freq11)+1:end);
            %             C1112 = C1112(numel(C1112)-numel(C11)+1:end);
            %             phi1112 = phi1112(numel(phi1112)-numel(C11)+1:end);
            %             phi1112Cu = phi1112Cu(numel(phi1112Cu)-numel(C11)+1:end);
            %             phi1112Cl = phi1112Cl(numel(phi1112Cl)-numel(C11)+1:end);
            
            hold on, plot(freq_1sthalf, coh_1sthalf, '.-r');
            hold on, plot(freq_2ndhalf, coh_2ndhalf, '.-m');
            %q.freq = freq11;
            %q.coh =  C11;
            q = catstruct(q, ...
              cell2struct({freq, coh, phase, coh_conf, phase_confU, phase_confL}, ...
              {'freq','coh', 'phase', 'coh_conf', 'phase_confU', 'phase_confL'}, 2));
            
            q.coh_halves_freq = freq_1sthalf; assert(numel(freq_1sthalf) == numel(freq_2ndhalf) && max(abs(freq_1sthalf - freq_2ndhalf)) < 1e-9)
            q.coh_halves = [torow(coh_1sthalf); torow(coh_2ndhalf)];
            q.phase_halves = [torow(phase_1sthalf); torow(phase_2ndhalf)];
            
          else
            % NOTE THAT if there are more than 2 groups (e.g. 3 probes), only the last one is saved, i.e. effectively more than 2
            % groups are NOT supported now.
            [freq1, C1, phi1, C1c, phi1u, phi1l] = ...
              freqDependentWindowCoherence(sum(spkAllShanks{grp}(:, Tstart:Tend), 1)', spk(Tstart:Tend)', 1/params.Fs, [], opt);
            hold on, plot(freq1, C1, '.-k');
            q = catstruct(q, ...
              cell2struct({freq1, C1, phi1, C1c, phi1u, phi1l}, {'freq_oGrp','C_oGrp', 'phi_oGrp', 'Cc_oGrp', 'phiCu_oGrp', 'phiCl_oGrp'}, 2)); % coherence and phase with population rate of distal shank group
          end
        end % looping over shank/tetrode groups
        
        assert(numel(q.psd_freq) == numel(q.freq) && max(abs(q.psd_freq - q.freq)) < 1e-9)
        
        subplot(5, min(UPF,numel(units)), whichSubplot + 2*min(UPF,numel(units)));
        hold on, plot(freq_1sthalf, bestUnwrap(phase_1sthalf), '.-r')
        hold on, plot(freq_2ndhalf, bestUnwrap(phase_2ndhalf), '.-m')
        if length(spkAllShanks) == 2
          hold on, plot(freq1, bestUnwrap(phi1), '.-k')
        end
        %         figTmp = plotPhaseSpectrum(freq11', [C11; C12; C1112], [phi11; phi12; phi1112], ...
        %           [NaN(size(phi11)); NaN(size(phi12)); phi1112Cu], [NaN(size(phi11)); NaN(size(phi12)); phi1112Cl]);
        figTmp = plotPhaseSpectrum(freq', coh, phase, phase_confU, phase_confL, coh_conf);
        
        % Compute the rate adjustment factors for coherence report, per Aoi et al.
        % we adjust as if the neuron's rate is 1 spk/s
        q.rateadjust_kappa = (1 + (q.mfr - 1)*q.mfr ./ q.psd);
        q.rateadjust_kappa(q.rateadjust_kappa < 0) = NaN; % the values are not supposed to be negative, it sometimes happens for neurons with low firing rate,
        % presumably because of the effect refractory period has on q.psd
        q.rateadjust_kappa = q.rateadjust_kappa.^-0.5;
        assert(isreal(q.psd) && isreal(coh) && isreal(q.rateadjust_kappa))
        
        subplot(4,2,2)
        hold on, semilogx(q.freq, coh .* q.rateadjust_kappa, 'b.-')
        hold on, semilogx(q.freq, (coh + coh_conf) .* q.rateadjust_kappa, 'c--')
        hold on, semilogx(q.freq, max(0,(coh - coh_conf) .* q.rateadjust_kappa), 'c--')
        ylabel('coh, rate-adjusted coh')
        title([db(dbCount).entryName '_' num2str(sh) '_u' num2str(units(u))], ...
          'Interpreter', 'none')
%f        hgsave(figTmp, [figsubdirname filesep db(dbCount).entryName '_' num2str(sh) ...
%f          '_PSDpolar_u' num2str(units(u))]);
        close(figTmp);
        figTmp = plotPhaseSpectrum(freq_s', coh_s, phase_s, phase_confU_s,  phase_confL_s, coh_conf_s);
%f        hgsave(figTmp, [figsubdirname filesep db(dbCount).entryName '_' num2str(sh) ...
%f          '_PSDpolar_u' num2str(units(u)) '_shuf']);
        close(figTmp);
        if length(spkAllShanks) == 2
          figTmp = plotPhaseSpectrum(q.freq_oGrp', q.C_oGrp, q.phi_oGrp, q.phiCu_oGrp, q.phiCl_oGrp, q.Cc_oGrp);
%f          hgsave(figTmp, [figsubdirname filesep db(dbCount).entryName '_' num2str(sh) ...
%f            '_PSDpolar_u' num2str(units(u)) '_oGrp']);
          close(figTmp);
        end
        figure(figH);
      end % different handling on single-shank and multiple-shank recordings
      
      %opt.freq = [0.005 25];
      %[freq, C, phi] = dumbCoherency(spk', spkAll'-spk', 100, -0.01, opt);
      %hold on, plot(freq, C, 'm', 'LineWidth', 2);
      
      %       subplot(4, numel(units), u + 2*numel(units)); % Coherence phases
      %       hold on, plot(freq1, phi1, '.-k');
      %       hold on, plot(freq2, phi2, '.-b');

% ADJUST GRAPHS
      subplot(5, min(UPF,numel(units)), whichSubplot);
      set(gca, 'XScale', 'log');
      set(gca, 'XTick', [0.03 0.1 0.3 1 3 10 30]);
      xlim(x_lim);
      
      title(sprintf('u %d, mfrs: %2.1f %2.1f spk/s', units(u), ...
        q.mfr_1sthalf, q.mfr_2ndhalf), 'Interpreter', 'none');
      
      subplot(5, min(UPF,numel(units)), whichSubplot + min(UPF,numel(units)));
      set(gca, 'XScale', 'log');
      set(gca, 'XTick', [0.03 0.1 0.3 1 3 10 30]);
      xlim(x_lim);
      ylim([0 1]);
      
      subplot(5, min(UPF,numel(units)), whichSubplot + 2*min(UPF,numel(units)));
      set(gca, 'XScale', 'log');
      set(gca, 'XTick', [0.03 0.1 0.3 1 3 10 30]);
      xlim(x_lim);

% CALCULATE AND PLOT SPIKE-TRIGERRED POPULATION RATE
      % the good old stPR (spike-triggered Population rate)
      %subplot(5, numel(units), u + 3*numel(units));
      subplot(5, min(UPF,numel(units)), whichSubplot + 4*min(UPF,numel(units)));
      for grp = 1:length(spkAllShanks)
        if grp == currentGrp % stPR with the population rate on the same probe
          stPR = xcorr(spk(Tstart:Tend)', sum(spkAllShanks{grp}(:, Tstart:Tend), 1)'-spkAll(Tstart:Tend)', 1e2)/sum(spk(Tstart:Tend));
          stPR = stPR / mean(stPR([1:10 190:201]));
          hold on, plot(-100:100, stPR, 'r');
          q.stPR = stPR;
        else % stPR with the population rate on the other probe
          stPR = xcorr(spk(Tstart:Tend)', sum(spkAllShanks{grp}(:, Tstart:Tend), 1)', 1e2)/sum(spk(Tstart:Tend));
          stPR = stPR / mean(stPR([1:10 190:201]));
          hold on, plot(-100:100, stPR, 'g');
        end
      end
      % stPR with the local (same shank/tetrode) population rate
      stPR = xcorr(spk(Tstart:round((Tend+Tstart)/2))', ...
        spkAll(Tstart:round((Tend+Tstart)/2))'-spk(Tstart:round((Tend+Tstart)/2))', 1e2)/sum(spk(Tstart:round((Tend+Tstart)/2)));
      stPR = stPR / mean(stPR([1:10 190:201]));
      %stPR = stPR - mean(stPR([1:10 190:201]));
      hold on, plot(-100:100, stPR, 'b');
      q.stPRl(1,:) = stPR; % stPR "local"
      stPR = xcorr(spk(1+round((Tend+Tstart)/2):Tend)', ...
        spkAll(1+round((Tend+Tstart)/2):Tend)'-spk(1+round((Tend+Tstart)/2):Tend)', 1e2)/sum(spk(1+round((Tend+Tstart)/2):Tend));
      stPR = stPR / mean(stPR([1:10 190:201]));
      hold on, plot(-100:100, stPR, 'b');
      q.stPRl(2,:) = stPR;
      xlim(100*[-1 1]);
      set(gca, 'XTick', -100:50:100);
      set(gca, 'XTickLabel', {num2str(-100/params.Fs), num2str(-50/params.Fs), '0', num2str(50/params.Fs), num2str(100/params.Fs)})
      
      % Mean vs var using different bin sizes (computed but not plotted)
      subplot(5, min(UPF,numel(units)), whichSubplot + 3*min(UPF,numel(units)));
      hold on;
      mv1 = zeros(14, 2);
      mv2 = mv1;
      m = min(mean(spk(Tstart:round((Tend+Tstart)/2))), mean(spk(1+round((Tend+Tstart)/2):Tend)));
      
      spk = spk(Tstart:Tend); % spk is anyhow changing in each iteration of the for-loop below:
      for s = 1:14 % starting from 10ms bins multiplying by 2 each step
        mv1(s, :) = [mean(spk(1:round(numel(spk)/2))), var(spk(1:round(numel(spk)/2)))];
        mv2(s, :) = [mean(spk(1+round(numel(spk)/2):end)), var(spk(1+round(numel(spk)/2):end))];
        spk = spk(1:floor(length(spk)/2)*2); spk = sum(reshape(spk, 2, []))';
      end
      q.mean_var = cat(3, mv1, mv2);
      q.mean_var_timeBin = params.Fs;
      hold on, plot(log(mv1(:,1)), log(mv1(:,2)) , '*-r');
      hold on, plot(log(mv2(:,1)), log(mv2(:,2)) , '*-m');
      hold on, plot([-2 10], [-2 10], 'g--');
      if whichSubplot == 1
        xlabel('log(mean)')
        ylabel('log(var)')
      end
      if mv1(end, 1)/mv2(end,1) > 2 || mv1(end, 1)/mv2(end,1) < 0.5 || ... means differ more than 100%
          mv1(end, 2)/(mv2(end,2) * mv1(end, 1)^2/mv2(end,1)^2) > 2 || mv1(end, 2)/(mv2(end,2) * mv1(end, 1)^2/mv2(end,1)^2) < 0.5 %variance differ > 100%
        title(['UNSTABLE ' num2str(mv1(end, 2)/(mv2(end,2) * mv1(end, 1)^2/mv2(end,1)^2))]);
      else
        b = regress([log(mv1(end-4:end,2)); log(mv2(end-4:end,2))], [[log(mv1(end-4:end,1)); log(mv2(end-4:end,1))] ones(10, 1)]);
        title(num2str(b(1)));
      end
      
      
      %       hold on, plot(mv1(:,1), sqrt(mv1(:,2))./mv1(:,1), '*-r');
      %       hold on, plot(mv2(:,1), sqrt(mv2(:,2))./mv2(:,1), '*-m');
      %       xlabel('mean'); ylabel('std/mean');
      
      
      if u == numel(units) || ~mod(u, UPF) 
        % have to save a figure (one figure per UPF units)
%f        hgsave(figH, [ figsubdirname filesep db(dbCount).entryName '_' num2str(sh) ...
%f          '_PSD_p' num2str(ceil(u/UPF))])
        close all
        figH = figure;
      end
      if ~exist('Q', 'var')
        Q = q;
      else
        Q(end+1) = q;
      end
      clear q
      
    end % loop on units
%f    save([figsubdirname filesep db(dbCount).entryName '.mat'], 'Q', 'opt'); % to be able to look at the data before all shanks are done...
  end % loop on shanks
%f  save([figsubdirname filesep db(dbCount).entryName '.mat'], 'Q', 'opt');
  clear Q
end % loop on db entries

