% Run this script to perform coherence and phase analyses on mouse data (old; don't use).


pool = gcp('nocreate');
delete(pool);
fclose all;
close all
% clear all
clc

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
  load(shankDataFile.name,'db','figsubdirname','UPF','x_lim','params','opt','shankIDs');
end

if ~exist('db', 'var')
  error('Please run some makedb script');
end

if ~exist('params', 'var')
  params.Fs = 400;
end
if ~exist('opt', 'var')
  opt.maxFreq = 200; %110;
  opt.winfactor = 10;
  opt.freqfactor = 1.333;
  opt.tapers = 5;
end
if opt.maxFreq > params.Fs/2
  error('opt.maxFreq > Nyquist?!')
end
  
if ~exist('figsubdirname', 'var')
  figsubdirname = '';
end
if ~exist('UPF', 'var')
  UPF = 5; % Units Per Figure
end
if ~exist('x_lim', 'var')
  x_lim = [0.02 50];
end

plotFigs1 = true;
plotFigs2 = true;

% INITIALISE PARALLEL POOL
pool = parpool;
feature('numcores');

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
          badUnits(end+1) = units(u); %#ok<SAGROW>
        end
        if isempty(strfind(dbID, 'acute')) %#ok<STREMP>
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
      selectedUnits = units;
      load_opt.selectedUnits = [];
      for u = 1:numel(selectedUnits)
        if any(selectedUnits(u) == units)
          load_opt.selectedUnits = [load_opt.selectedUnits selectedUnits(u)];
        end
      end
      spk = loadAsRasterSparse(baseFilename, sh, SR, 1/params.Fs, load_opt);
      uCentres = zeros(numel(load_opt.selectedUnits),1);
      load([baseFilename '.chanMap.mat']);
      %load([baseFilename '.imec.ap_CAR' '.chanMap.mat']);
      for u = 1:numel(load_opt.selectedUnits)
        i = find(chanMap(:,1) == load_opt.selectedUnits(u));
        if ~isempty(i)
          uCentres(u) = chanMap(i,2);
        end
      end
      setup = load([topDir filesep db(dbCount).animal filesep num2str(db(dbCount).series) filesep...
        'forPRBimecP3opt3.mat']);
      sites = setup.chanMap;
      connectedSites = setup.connected;
      xcoords = setup.xcoords;
      ycoords = setup.ycoords;
      save([shankDataFile '.mat'],'spk','db','SR','dbID','metadb','units',...
        'opt','params','topDir','UPF','x_lim','shankMUAs','spkAll',...
        'spkAllShanks','figsubdirname','shankIDs','baseFilename',...
        'currentGrp','load_opt','uCentres','shankDataFile','sites',...
        'connectedSites','xcoords','ycoords','chanMap');
    end
    
% LOOP THROUGH UNITS
    parfor u = 1:numel(load_opt.selectedUnits)
      if ~any(units == load_opt.selectedUnits(u)) %#ok<PFBNS>
        continue
      end
      fprintf('Started processing unit %i\n',load_opt.selectedUnits(u));
      
% CREATE SUMMARY DATA STRUCTURE
      q = cell2struct({db(dbCount).animal, db(dbCount).series, db(dbCount).entryName, sh, load_opt.selectedUnits(u)}, ...
        {'animal', 'series', 'dbentryName', 'shank', 'unit'}, 2);

% SPIKES FOR UNIT OF INTEREST
      if isempty(strfind(dbID, 'acute')) %#ok<STREMP>
        Tstart = 1;
        Tend = numel(spkAll);
        % Tstart = T(1); % we presume it's a chronic recording ==> assume stationary throughout
        % Tend = T(end);
      else
        error('not handled anymore')
        % [Tstart, Tend] = findStationaryInterval(T, SR);
      end
      Tmid = round((Tend+Tstart)/2);
      spkOI = full(spk(u,Tstart:Tend)); %#ok<PFBNS>
      
% SPIKES FOR THE VICINITY OF UNIT OF INTEREST (in this case the whole shank)
      if isempty(strfind(dbID, 'IMEC')) %#ok<STREMP>
        localMUA = []; % this is not an IMEC probe
      else
%         localMUA = loadMUAaroundcontact([topDir filesep db(dbCount).animal filesep num2str(db(dbCount).series)], units(u), SR, 1/params.Fs);
%         assert(all(spkOI <= localMUA) && all(localMUA <= spkAll))
        localMUA = spkAll;
      end
      assert(all(spkOI <= spkAll))

% SUMMARY DATA
      % Mean firing rate (mfr)
      q.mfr = mean(spkOI)*params.Fs; %#ok<PFBNS>
      q.mfr_1sthalf = mean(spkOI(Tstart:Tmid))*params.Fs;
      q.mfr_2ndhalf = mean(spkOI(Tmid+1:Tend))*params.Fs;
      
      % PSD
      [freq_1sthalf, psd_1sthalf] = freqDependentWindowCoherence(spkOI(Tstart:Tmid), [], 1/params.Fs, [], opt);
      [freq_2ndhalf, psd_2ndhalf] = freqDependentWindowCoherence(spkOI(Tmid+1:Tend), [], 1/params.Fs, [], opt);
      [freq, psd, ~, psd_conf] = freqDependentWindowCoherence(spkOI(Tstart:Tend), [], 1/params.Fs, [], opt);
      q.psd_halves = [psd_1sthalf; psd_2ndhalf];
      q.psd_halves_freq = freq_1sthalf; assert(numel(freq_1sthalf) == numel(freq_2ndhalf) && max(abs(freq_1sthalf - freq_2ndhalf)) < 1e-9)
      q.psd_freq = freq; q.psd = psd; q.psd_conf = psd_conf; q.psd_numelSignal = Tend - Tstart + 1;
      
      % Coherence with MUA
      if length(spkAllShanks) == 1 && length(db(dbCount).shank) == 1 && ~isempty(localMUA) % all we've got is one shank, which we handle only for IMEC recordings
        [freq_1sthalf, coh_1sthalf, phase_1sthalf] = ...
          freqDependentWindowCoherence(spkAll(Tstart:Tmid)'-spkOI(Tstart:Tmid)', spkOI(Tstart:Tmid)', 1/params.Fs, [], opt);
        [freq_2ndhalf, coh_2ndhalf, phase_2ndhalf] = ...
          freqDependentWindowCoherence(spkAll(Tmid+1:Tend)'-spkOI(Tmid+1:Tend)', spkOI(Tmid+1:Tend)', 1/params.Fs, [], opt);
        [freq, coh, phase, coh_conf, phase_confU, phase_confL] = freqDependentWindowCoherence(spkAll(Tstart:Tend)'-spkOI', spkOI', 1/params.Fs, [], opt);
        q = catstruct(q, cell2struct({freq, coh, phase, coh_conf, phase_confU, phase_confL}, ...
              {'freq', 'coh', 'phase', 'coh_conf', 'phase_confU', 'phase_confL'}, 2));            
        q.coh_halves_freq = freq_1sthalf; assert(numel(freq_1sthalf) == numel(freq_2ndhalf) && max(abs(freq_1sthalf - freq_2ndhalf)) < 1e-9)
        q.coh_halves = [torow(coh_1sthalf); torow(coh_2ndhalf)];
        q.phase_halves = [torow(phase_1sthalf); torow(phase_2ndhalf)];    
            
        % Compute the rate adjustment factors for coherence report, per Aoi et al.
        % we adjust as if the neuron's rate is 1 spk/s
        q.rateadjust_kappa = coherenceAdjustment(q.mfr, 1, q.psd, 0, 1);
        assert(isreal(q.psd) && isreal(coh) && isreal(q.rateadjust_kappa))

        % Population mean firing rate
        q.mfr_popRate = mean(spkAll(Tstart:Tend))*params.Fs;
        q.mfr_localMUARate = mean(localMUA(Tstart:Tend))*params.Fs;
        
      else % we have more than 1 shank
        [freq2, C2, phi2, confC2] = freqDependentWindowCoherence(spkAll(Tstart:Tend)' - spkOI', spkOI', 1/params.Fs, [], opt);
% LOOP OVER SHANK GROUPS
        for grp = 1:length(spkAllShanks)
          if grp == currentGrp
            
            % Local mean firing rate
            q.mfr_popRate = mean(sum(spkAllShanks{grp}(:,Tstart:Tend),1)-spkAll(Tstart:Tend))*params.Fs; % group of shanks minus local shank
            q.mfr_localPopRate = mean(spkAll(Tstart:Tend) - spkOI)*params.Fs; % local shank minus unit
            assert(~any(spkAll > sum(spkAllShanks{grp}, 1)))
            
            % Coherence
            % compute coherence with the local (i.e. on the same shank/tetrode) population rate. If spike sorting excludes the
            % possibility of more than 1 spike in (say) 1ms detection window, a spurious correlation might be introduced on fast timescale
            [~, coh, phase, coh_conf, phase_confU, phase_confL] = ... we skip the first output (freq) because it's the same as that below
              freqDependentWindowCoherence(spkAll(Tstart:Tend)' - spkOI', spkOI', 1/params.Fs, [], opt);
            q = catstruct(q, cell2struct({coh, phase, coh_conf, phase_confU, phase_confL}, ...
              {'coh_local', 'phase_local', 'coh_conf_local', 'phase_confU_local', 'phase_confL_local'}, 2));
            
            % compute coherence where population rate excludes the local shank/tetrode
            [freq_1sthalf, coh_1sthalf, phase_1sthalf] = ...
              freqDependentWindowCoherence(sum(spkAllShanks{grp}(:,Tstart:Tmid),1)'-spkAll(Tstart:Tmid)', spkOI(Tstart:Tmid)', 1/params.Fs, [], opt);
            [freq_2ndhalf, coh_2ndhalf, phase_2ndhalf] = ...
              freqDependentWindowCoherence(sum(spkAllShanks{grp}(:,Tmid+1:Tend),1)'-spkAll(Tmid+1:Tend)', spkOI(Tmid+1:Tend)', 1/params.Fs, [], opt);
            
            [freq, coh, phase, coh_conf, phase_confU, phase_confL] = ...
              freqDependentWindowCoherence(sum(spkAllShanks{grp}(:,Tstart:Tend),1)'-spkAll(Tstart:Tend)', spkOI', 1/params.Fs, [], opt);
%             [freq_s, coh_s, phase_s, coh_conf_s, phase_confU_s, phase_confL_s] = ... shuffled (spk is taken from end to start)
%               freqDependentWindowCoherence(sum(spkAllShanks{grp}(:,Tstart:Tend),1)'-spkAll(Tstart:Tend)', spkOI(Tend:-1:Tstart)', 1/params.Fs, [], opt);
            q = catstruct(q, cell2struct({freq, coh, phase, coh_conf, phase_confU, phase_confL}, ...
              {'freq','coh', 'phase', 'coh_conf', 'phase_confU', 'phase_confL'}, 2));
            q.coh_halves_freq = freq_1sthalf; assert(numel(freq_1sthalf) == numel(freq_2ndhalf) && max(abs(freq_1sthalf - freq_2ndhalf)) < 1e-9)
            q.coh_halves = [torow(coh_1sthalf); torow(coh_2ndhalf)];
            q.phase_halves = [torow(phase_1sthalf); torow(phase_2ndhalf)];
            
          else
            % NOTE THAT if there are more than 2 groups (e.g. 3 probes), only the last one is saved,
            % i.e. effectively more than 2 groups are NOT supported now.
            % coherence and phase with population rate of the last shank group
            [freq1, C1, phi1, C1c, phi1u, phi1l] = ...
              freqDependentWindowCoherence(sum(spkAllShanks{grp}(:,Tstart:Tend),1)', spkOI', 1/params.Fs, [], opt);
            q = catstruct(q, cell2struct({freq1, C1, phi1, C1c, phi1u, phi1l}, ...
              {'freq_oGrp','C_oGrp', 'phi_oGrp', 'Cc_oGrp', 'phiCu_oGrp', 'phiCl_oGrp'}, 2));
          end
        end % loop over shank/tetrode groups
        
        assert(numel(q.psd_freq) == numel(q.freq) && max(abs(q.psd_freq - q.freq)) < 1e-9)
        
        % Compute the rate adjustment factors for coherence report, per Aoi et al.
        % we adjust as if the neuron's rate is 1 spk/s
        q.rateadjust_kappa = coherenceAdjustment(q.mfr, 1, q.psd, 0, 1);
        assert(isreal(q.psd) && isreal(coh) && isreal(q.rateadjust_kappa))
      end % different handling on single-shank and multiple-shank recordings

% CALCULATE SPIKE-TRIGERRED POPULATION RATE (stPR)
      for grp = 1:length(spkAllShanks)
        if grp == currentGrp % stPR with the population rate on the same probe
          stPR = xcorr(spkOI', sum(spkAllShanks{grp}(:,Tstart:Tend),1)'-spkAll(Tstart:Tend)', 1e2)/sum(spkOI);
          stPR = stPR / mean(stPR([1:10 190:201]));
          q.stPRgr0 = stPR;
        else % stPR with the population rate on the other probe
          % NOTE THAT if there are more than 2 groups (e.g. 3 probes), only the last one is saved,
          % i.e. effectively more than 2 groups are NOT supported now.
          stPR = xcorr(spkOI', sum(spkAllShanks{grp}(:,Tstart:Tend),1)', 1e2)/sum(spkOI);
          stPR = stPR / mean(stPR([1:10 190:201]));
          q.stPRgr1 = stPR;
        end
      end
      % stPR with the local (same shank/tetrode) population rate
      stPR = xcorr(spkOI(Tstart:Tmid)', spkAll(Tstart:Tmid)'-spkOI(Tstart:Tmid)', 1e2)/sum(spkOI(Tstart:Tmid));
      stPR = stPR / mean(stPR([1:10 190:201]));
      q.stPR(1,:) = stPR;
      stPR = xcorr(spkOI(Tmid+1:Tend)', spkAll(Tmid+1:Tend)'-spkOI(Tmid+1:Tend)', 1e2)/sum(spkOI(Tmid+1:Tend));
      stPR = stPR / mean(stPR([1:10 190:201]));
      q.stPR(2,:) = stPR;
      
      % Mean vs var using different bin sizes
      mv1 = zeros(14, 2);
      mv2 = mv1;
      q.m = min(mean(spkOI(Tstart:Tmid)), mean(spkOI(Tmid+1:Tend)));
      
      for s = 1:14 % starting from 10ms bins multiplying by 2 each step
        mv1(s, :) = [mean(spkOI(1:round(numel(spkOI)/2))), var(spkOI(1:round(numel(spkOI)/2)))];
        mv2(s, :) = [mean(spkOI(1+round(numel(spkOI)/2):end)), var(spkOI(1+round(numel(spkOI)/2):end))];
        spkOI = spkOI(1:floor(length(spkOI)/2)*2); spkOI = sum(reshape(spkOI, 2, []))';
      end
      q.mv1 = mv1;
      q.mv2 = mv2;
      q.mean_var = cat(3, mv1, mv2);
      q.mean_var_timeBin = params.Fs;
%       if mv1(end, 1)/mv2(end,1) > 2 || mv1(end, 1)/mv2(end,1) < 0.5 || ... means differ more than 100%
%           mv1(end, 2)/(mv2(end,2) * mv1(end, 1)^2/mv2(end,1)^2) > 2 || mv1(end, 2)/(mv2(end,2) * mv1(end, 1)^2/mv2(end,1)^2) < 0.5 %variance differ > 100%
%       else
%         b = regress([log(mv1(end-4:end,2)); log(mv2(end-4:end,2))], [[log(mv1(end-4:end,1)); log(mv2(end-4:end,1))] ones(10, 1)]);
%       end
      
% SAVE UNIT SUMMARY DATA
      saveParfor([shankDataFile '_' num2str(load_opt.selectedUnits(u)) '.mat'], q);
      fprintf('Finished processing unit %i\n',load_opt.selectedUnits(u));
      
    end % loop over units
    
% SAVE SHANK SUMMARY DATA
    fileList = dir([shankDataFile '_*.mat']);
    if ~isempty(fileList)
      Q = struct([]);
      for i = 1:size(fileList,1)
        q = load(fileList(i).name);
        Q{end+1} = q.q; %#ok<SAGROW>
      end
      save(shankDataFile,'Q','-append')
      delete([shankDataFile '_*.mat']);
    end
    
    fprintf('Finished processing shank %i\n',sh);
  end % loop over shanks
  
  fprintf('Finished processing db entry %i\n',dbCount);
end % loop over db entries
%delete(pool);

% PLOT SUMMARY DATA
if plotFigs1
  for dbCount = length(db):-1:1
    for sh = shankIDs
      shankDataFile = ['spk' num2str(dbCount) '_' num2str(sh)];
      if exist([shankDataFile '.mat'],'file')
        load([shankDataFile '.mat']);
      else
        continue
      end
      nUnits = numel(Q);
      for u = 1:nUnits
        q = Q{u};
        
        if mod(u, UPF) == 1
          figH = figure('units', 'normalized', 'position', [0.002, .04, 1, .88]);
        end
        if numel(Q) < UPF
          whichSubplot = u;
        else
          whichSubplot = mod(u, UPF); if ~whichSubplot; whichSubplot=UPF; end
        end
        
        % PSD
        subplot(5, min(UPF,numel(Q)), whichSubplot);
        hold on, plot(q.psd_halves_freq, q.psd_halves(1,:), 'r');
        hold on, plot(q.psd_halves_freq, q.psd_halves(2,:), 'm');
        titleStr = sprintf('sh%i_u%i_mfr1_%3.1f_mfr2_%3.1f', sh, q.unit, q.mfr_1sthalf, q.mfr_2ndhalf);
        title(titleStr, 'Interpreter', 'none')
        if whichSubplot == 1
          uFig = zeros(1,UPF);
          ylabel('Power (APs^2/Hz)')
        end
        
        % Coherence with MUA
        subplot(5, min(UPF,numel(Q)), whichSubplot + min(UPF,numel(Q)));
        if length(spkAllShanks) == 1 && length(db(dbCount).shank) == 1 && ~isempty(spkAll) % all we've got is one shank, which we handle only for IMEC recordings
          hold on, plot(q.coh_halves_freq, bestUnwrap(q.coh_halves(1,:)), '.-r');
          hold on, plot(q.coh_halves_freq, bestUnwrap(q.coh_halves(2,:)), '.-m');
          if whichSubplot == 1
            ylabel('Coherence')
          end
          
          % Phase
          subplot(5, min(UPF,numel(Q)), whichSubplot + 2*min(UPF,numel(Q)));
          hold on, plot(q.coh_halves_freq, q.phase_halves(1,:), '.-r');
          hold on, plot(q.coh_halves_freq, q.phase_halves(2,:), '.-m');
          xlabel('Frequency (Hz)')
          if whichSubplot == 1
            ylabel('Phase (rad)')
          end
          
          % Circular phase
          figTmp = plotPhaseSpectrum(q.freq', q.coh, q.phase, q.phase_confU, q.phase_confL, q.coh_conf);
          
          % Rate-adjusted coherence
          subplot(4,2,2)
          hold on, semilogx(q.freq, q.coh .* q.rateadjust_kappa, 'b.-')
          hold on, semilogx(q.freq, (q.coh + q.coh_conf) .* q.rateadjust_kappa, 'c--')
          hold on, semilogx(q.freq, max(0,(q.coh - q.coh_conf) .* q.rateadjust_kappa), 'c--')
          ylabel('coh, rate-adjusted coh')
          title([db(dbCount).animal '_' num2str(db(dbCount).series) '_' num2str(sh) '_u' num2str(units(u))], 'Interpreter', 'none')
          if strcmpi(figsubdirname,'')
            hgsave(figTmp, [db(dbCount).animal '_' num2str(db(dbCount).series) '_' 'PSDpolar_u' num2str(q.unit)]);
          else
            hgsave(figTmp, [figsubdirname filesep db(dbCount).animal '_' num2str(db(dbCount).series) '_' ...
              'PSDpolar_u' num2str(q.unit)]);
          end
          close(figTmp);
        else
          error('Data plotting for more than one shank is not supported by the current version of the program.')
        end
        
        subplot(5, min(UPF,numel(Q)), whichSubplot);
        set(gca, 'XScale', 'log');
        set(gca, 'XTick', [0.03 0.1 0.3 1 3 10 30]);
        xlim(x_lim);
        
        subplot(5, min(UPF,numel(Q)), whichSubplot + min(UPF,numel(Q)));
        set(gca, 'XScale', 'log');
        set(gca, 'XTick', [0.03 0.1 0.3 1 3 10 30]);
        xlim(x_lim);
        ylim([0 1]);
        
        subplot(5, min(UPF,numel(Q)), whichSubplot + 2*min(UPF,numel(Q)));
        set(gca, 'XScale', 'log');
        set(gca, 'XTick', [0.03 0.1 0.3 1 3 10 30]);
        xlim(x_lim);
        
        % Mean vs var using different bin sizes
        subplot(5, min(UPF,numel(Q)), whichSubplot + 3*min(UPF,numel(Q)));
        hold on, plot(log(q.mv1(:,1)), log(q.mv1(:,2)) , '*-r');
        hold on, plot(log(q.mv2(:,1)), log(q.mv2(:,2)) , '*-m');
        hold on, plot([-2 10], [-2 10], 'g--');
        xlabel('log(mean)')
        if whichSubplot == 1
          ylabel('log(var)')
        end
        
        % stPR
        subplot(5, min(UPF,numel(Q)), whichSubplot + 4*min(UPF,numel(Q)));
        hold on, plot(-100:100, q.stPRgr0, 'r');
        hold on, plot(-100:100, q.stPR(1,:), 'b');
        hold on, plot(-100:100, q.stPR(2,:), 'b');
        xlim(100*[-1 1]);
        xlabel('Time (ms)')
        if whichSubplot == 1
          ylabel('XCorr (1/mean)')
        end
        set(gca, 'XTick', -100:50:100);
        set(gca, 'XTickLabel', {num2str(-100/params.Fs), num2str(-50/params.Fs), '0', num2str(50/params.Fs), num2str(100/params.Fs)})
        
        uFig(whichSubplot) = q.unit;
        if u == nUnits || ~mod(u, UPF)
          suffix = [];
          for i = 1:numel(uFig)
            suffix = [suffix sprintf('_%i',uFig(i))]; %#ok<AGROW>
          end
          if strcmpi(figsubdirname,'')
            hgsave(figH, [db(dbCount).animal '_' num2str(db(dbCount).series) '_' 'PSD_p' num2str(ceil(u/UPF)) suffix])
          else
            hgsave(figH, [figsubdirname filesep db(dbCount).animal '_' num2str(db(dbCount).series) '_' ...
              'PSD_p' num2str(ceil(u/UPF)) suffix])
          end
          close(figH);
        end
      end % loop over units
    end % loop over shanks
  end % loop over db entries
end



%% Compare electrode segments

% CREATE ELECTRODE MAP
vDim = round(numel(sites)/2);
hDim = 4;
basicMotive = [1 0 1 0 0 1 0 1];
eMapExt = repmat(basicMotive,1,vDim/2);
siteIDs = sites;
siteIDs(~connectedSites) = NaN;
eMap1(logical(eMapExt)) = siteIDs';
eMap1(~logical(eMapExt)) = NaN;
eMap1 = reshape(eMap1,hDim,vDim);
eMap1 = rot90(eMap1',2);
eMap2 = eMap1(:,[2 1 4 3]);
eMap3 = fliplr(eMap1);
eMap4 = fliplr(eMap2);
eMap = eMap2;
save(shankDataFile,'eMap','eMap1','eMap2','eMap3','eMap4','-append')

% ASSIGN ELECTRODES TO SEGMENTS AND OBTAIN COMPARISON PAIRS
parts = [2 3 4 6 8 12];
segCount = 0;
pairs = [];
for nP = 1:numel(parts)
  vDimP = round(vDim/parts(nP));
  pairs = [pairs; nchoosek(segCount+1:segCount+parts(nP),2)]; %#ok<AGROW>
  for s = parts(nP):-1:1
    segCount = segCount + 1;
    segs{segCount} = sort(reshape(eMap(vDimP*s:-1:vDimP*(s-1)+1,:), [1 vDimP*hDim]));
    segs{segCount} = segs{segCount}(~isnan(segs{segCount}));
  end
end

% LOAD CHANNEL MAP AND LOOP OVER SEGMENT PAIRS
assert(numel(load_opt.selectedUnits) == numel(units))
assert(max(chanMap(:,2)) <= max(max(eMap)))
for dbCount = length(db):-1:1
  for sh = shankIDs
    shankDataFile = ['spk' num2str(dbCount) '_' num2str(sh)];
    load([shankDataFile '.mat']);
    parfor p = 1:size(pairs,1)
      fprintf('Started processing pair %i out of %i\n',p,size(pairs,1));
      
% CREATE SUMMARY DATA STRUCTURE
      qp = cell2struct({db(dbCount).animal, db(dbCount).series, db(dbCount).entryName, sh, p}, ...
        {'animal', 'series', 'dbentryName', 'shank', 'pair'}, 2);
      
% SPIKES FOR THE PAIR SEGMENTS
      if isempty(strfind(dbID, 'acute')) %#ok<STREMP>
        Tstart = 1;
        Tend = numel(spkAll);
        % Tstart = T(1); % we presume it's a chronic recording ==> assume stationary throughout
        % Tend = T(end);
      else
        error('not handled anymore')
        % [Tstart, Tend] = findStationaryInterval(T, SR);
      end
      Tmid = round((Tend+Tstart)/2);
      segOI = segs{pairs(p,1)}; %#ok<PFBNS>
      qp.seg1 = segOI;
      segOI1 = [];
      for i = 1:numel(segOI)
        u = find(uCentres == segOI(i));
        if ~isempty(u)
          segOI1 = [segOI1; u];
        end
      end
      qp.units1 = segOI1;
      segOI1 = full(sum(spk(segOI1,Tstart:Tend),1)); %#ok<PFBNS>
      assert(all(segOI1 <= spkAll))
      segOI = segs{pairs(p,2)};
      qp.seg2 = segOI;
      segOI2 = [];
      for i = 1:numel(segOI)
        u = find(uCentres == segOI(i));
        if ~isempty(u)
          segOI2 = [segOI2; u];
        end
      end
      qp.units2 = segOI2;
      segOI2 = full(sum(spk(segOI2,Tstart:Tend),1));
      assert(all(segOI1 <= spkAll))
      
% SUMMARY DATA
      % Mean firing rate (mfr)
      qp.mfr1 = mean(segOI1)*params.Fs; %#ok<PFBNS>
      qp.mfr1_1sthalf = mean(segOI1(Tstart:Tmid))*params.Fs;
      qp.mfr1_2ndhalf = mean(segOI1(Tmid+1:Tend))*params.Fs;
      qp.mfr2 = mean(segOI2)*params.Fs;
      qp.mfr2_1sthalf = mean(segOI2(Tstart:Tmid))*params.Fs;
      qp.mfr2_2ndhalf = mean(segOI2(Tmid+1:Tend))*params.Fs;
      
      % PSD
      [freq1_1sthalf, psd1_1sthalf] = freqDependentWindowCoherence(segOI1(Tstart:Tmid), [], 1/params.Fs, [], opt);
      [freq1_2ndhalf, psd1_2ndhalf] = freqDependentWindowCoherence(segOI1(Tmid+1:Tend), [], 1/params.Fs, [], opt);
      [freq1, psd1, ~, psd1_conf] = freqDependentWindowCoherence(segOI1(Tstart:Tend), [], 1/params.Fs, [], opt);
      qp.psd1_halves = [psd1_1sthalf; psd1_2ndhalf];
      qp.psd1_halves_freq = freq1_1sthalf; assert(numel(freq1_1sthalf) == numel(freq1_2ndhalf) && max(abs(freq1_1sthalf - freq1_2ndhalf)) < 1e-9)
      qp.psd1_freq = freq1; qp.psd1 = psd1; qp.psd1_conf = psd1_conf; qp.psd1_numelSignal = Tend - Tstart + 1;
      [freq2_1sthalf, psd2_1sthalf] = freqDependentWindowCoherence(segOI2(Tstart:Tmid), [], 1/params.Fs, [], opt);
      [freq2_2ndhalf, psd2_2ndhalf] = freqDependentWindowCoherence(segOI2(Tmid+1:Tend), [], 1/params.Fs, [], opt);
      [freq2, psd2, ~, psd2_conf] = freqDependentWindowCoherence(segOI2(Tstart:Tend), [], 1/params.Fs, [], opt);
      qp.psd2_halves = [psd2_1sthalf; psd2_2ndhalf];
      qp.psd2_halves_freq = freq2_1sthalf; assert(numel(freq2_1sthalf) == numel(freq2_2ndhalf) && max(abs(freq2_1sthalf - freq2_2ndhalf)) < 1e-9)
      qp.psd2_freq = freq2; qp.psd2 = psd2; qp.psd2_conf = psd2_conf; qp.psd2_numelSignal = Tend - Tstart + 1;
      
      % Coherence
      [freq_1sthalf, coh_1sthalf, phase_1sthalf] = freqDependentWindowCoherence(segOI2(Tstart:Tmid)', segOI1(Tstart:Tmid)', 1/params.Fs, [], opt);
      [freq_2ndhalf, coh_2ndhalf, phase_2ndhalf] = freqDependentWindowCoherence(segOI2(Tmid+1:Tend)', segOI1(Tmid+1:Tend)', 1/params.Fs, [], opt);
      [freq, coh, phase, coh_conf, phase_confU, phase_confL] = freqDependentWindowCoherence(segOI2', segOI1', 1/params.Fs, [], opt);
      qp = catstruct(qp, cell2struct({freq, coh, phase, coh_conf, phase_confU, phase_confL}, ...
        {'freq', 'coh', 'phase', 'coh_conf', 'phase_confU', 'phase_confL'}, 2));
      qp.coh_halves_freq = freq_1sthalf; assert(numel(freq_1sthalf) == numel(freq_2ndhalf) && max(abs(freq_1sthalf - freq_2ndhalf)) < 1e-9)
      qp.coh_halves = [torow(coh_1sthalf); torow(coh_2ndhalf)];
      qp.phase_halves = [torow(phase_1sthalf); torow(phase_2ndhalf)];
      
      % Compute the rate adjustment factors for coherence report, per Aoi et al.
      qp.rateadjust_kappa1 = coherenceAdjustment(qp.mfr1, qp.mfr2, qp.psd1, 0, 1);
      assert(isreal(qp.psd1) && isreal(coh) && isreal(qp.rateadjust_kappa1))
      qp.rateadjust_kappa2 = coherenceAdjustment(qp.mfr2, qp.mfr1, qp.psd2, 0, 1);
      assert(isreal(qp.psd2) && isreal(coh) && isreal(qp.rateadjust_kappa2))
      
      % stPR
      stPR = xcorr(segOI1(Tstart:Tmid)', segOI2(Tstart:Tmid)', 1e2)/sum(segOI1(Tstart:Tmid));
      stPR = stPR / mean(stPR([1:10 190:201]));
      qp.stPR(1,:) = stPR;
      stPR = xcorr(segOI1(Tmid+1:Tend)', segOI2(Tmid+1:Tend)', 1e2)/sum(segOI1(Tmid+1:Tend));
      stPR = stPR / mean(stPR([1:10 190:201]));
      qp.stPR(2,:) = stPR;
      
      % Mean vs var using different bin sizes
      mv1 = zeros(14, 2);
      mv2 = mv1;
      qp.m_1 = min(mean(segOI1(Tstart:Tmid)), mean(segOI1(Tmid+1:Tend)));
      
      for s = 1:14 % starting from 10ms bins multiplying by 2 each step
        mv1(s, :) = [mean(segOI1(1:round(numel(segOI1)/2))), var(segOI1(1:round(numel(segOI1)/2)))];
        mv2(s, :) = [mean(segOI1(1+round(numel(segOI1)/2):end)), var(segOI1(1+round(numel(segOI1)/2):end))];
        segOI1 = segOI1(1:floor(length(segOI1)/2)*2); segOI1 = sum(reshape(segOI1, 2, []))';
      end
      qp.mv1_1 = mv1;
      qp.mv2_1 = mv2;
      qp.mean_var_1 = cat(3, mv1, mv2);
      qp.mean_var_timeBin_1 = params.Fs;
      %       if mv1(end, 1)/mv2(end,1) > 2 || mv1(end, 1)/mv2(end,1) < 0.5 || ... means differ more than 100%
      %           mv1(end, 2)/(mv2(end,2) * mv1(end, 1)^2/mv2(end,1)^2) > 2 || mv1(end, 2)/(mv2(end,2) * mv1(end, 1)^2/mv2(end,1)^2) < 0.5 %variance differ > 100%
      %       else
      %         b = regress([log(mv1(end-4:end,2)); log(mv2(end-4:end,2))], [[log(mv1(end-4:end,1)); log(mv2(end-4:end,1))] ones(10, 1)]);
      %       end
      
      mv1 = zeros(14, 2);
      mv2 = mv1;
      qp.m_2 = min(mean(segOI2(Tstart:Tmid)), mean(segOI2(Tmid+1:Tend)));
      
      for s = 1:14 % starting from 10ms bins multiplying by 2 each step
        mv1(s, :) = [mean(segOI2(1:round(numel(segOI2)/2))), var(segOI2(1:round(numel(segOI2)/2)))];
        mv2(s, :) = [mean(segOI2(1+round(numel(segOI2)/2):end)), var(segOI2(1+round(numel(segOI2)/2):end))];
        segOI2 = segOI2(1:floor(length(segOI2)/2)*2); segOI2 = sum(reshape(segOI2, 2, []))';
      end
      qp.mv1_2 = mv1;
      qp.mv2_2 = mv2;
      qp.mean_var_2 = cat(3, mv1, mv2);
      qp.mean_var_timeBin_2 = params.Fs;
      %       if mv1(end, 1)/mv2(end,1) > 2 || mv1(end, 1)/mv2(end,1) < 0.5 || ... means differ more than 100%
      %           mv1(end, 2)/(mv2(end,2) * mv1(end, 1)^2/mv2(end,1)^2) > 2 || mv1(end, 2)/(mv2(end,2) * mv1(end, 1)^2/mv2(end,1)^2) < 0.5 %variance differ > 100%
      %       else
      %         b = regress([log(mv1(end-4:end,2)); log(mv2(end-4:end,2))], [[log(mv1(end-4:end,1)); log(mv2(end-4:end,1))] ones(10, 1)]);
      %       end
      
% SAVE UNIT SUMMARY DATA
      saveParfor(['s_' shankDataFile '_' num2str(p) '.mat'], qp);
      fprintf('Finished processing pair %i out of %i\n',p,size(pairs,1));
      
    end % loop over segment pairs
    
% SAVE SHANK SUMMARY DATA
    fileList = dir(['s_' shankDataFile '_*.mat']);
    if ~isempty(fileList)
      QP = struct([]);
      for i = 1:size(fileList,1)
        qp = load(fileList(i).name);
        QP{end+1} = qp.q; %#ok<SAGROW>
      end
      save(shankDataFile,'QP','-append')
      delete(['s_' shankDataFile '_*.mat']);
    end
  end % loop over shanks
end % loop over db entries

% PLOT PHASE PREFERENCE FOR SEGMENT COMPARISONS
colours = [[0 0 255];    %blue
           [255 0 0];    %red
           [0 255 0];    %green
           [0 255 255];  %cyan
           [205 205 0];  %yellow
           [255 0 255]]; %magenta
colours = colours./255;

if plotFigs2
  for dbCount = length(db):-1:1
    for sh = shankIDs
      shankDataFile = ['spk' num2str(dbCount) '_' num2str(sh)];
      load([shankDataFile '.mat']);
      for p = 1:numel(QP)
        phase = QP{p}.phase;
        countSignificant = length(phase(~isnan(QP{p}.phase_confU)));
        countTotal = length(phase);
        titleStr = ['_' num2str(countSignificant) '_' num2str(countTotal) sprintf('mfr1_%4.1f_mfr2_%4.1f', QP{p}.mfr1, QP{p}.mfr2)];
        figStr = ['_' num2str(countSignificant) '_' num2str(countTotal)];

        prefix = [db(dbCount).animal '_' num2str(db(dbCount).series) '_pair' num2str(QP{p}.pair)];
        figH = phaseFigPairMouse(QP{p}, eMap, colours, titleStr, figStr, prefix);
        
        close(figH)
      end
    end
  end
end

movefile([shankDataFile '.mat'],[db(dbCount).animal '_' num2str(db(dbCount).series) '.mat'],'f');
