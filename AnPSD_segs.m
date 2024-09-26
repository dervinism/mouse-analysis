% A part of the old AnPSD script adapted to perform coherence and phase
% analyses comparing segment spiking rate data.


%% LOAD PRE-PROCESSED DATA
if ~exist('dataStruct','var')
  load(dataFile);
end


% INITIALISE PARALLEL POOL
delete(gcp('nocreate'))
pool = parpool;
feature('numcores');


% LOOP THROUGH DB ENTRIES
fnsData = fieldnames(dataStruct);
for dbCount = numel(fnsData):-1:1
  dbStruct = dataStruct.(fnsData{dbCount});
  shankIDs = fieldnames(dbStruct);

% LOOP THROUGH SHANKS
  for sh = 1:numel(shankIDs)
    shankStruct = dbStruct.(shankIDs{sh});
    fprintf('%s shank %d -------------------\n', fnsData{dbCount}, sh);
    
% LOAD THE CONTENTS OF THE SHANK STRUCTURE VARIABLE
    spk = shankStruct.spk;
    dbID = shankStruct.dbID;
    units = shankStruct.units;
    opt = shankStruct.opt;
    params = shankStruct.params;
    spkAll = shankStruct.spkAll;
    load_opt = shankStruct.load_opt;
    for dbCountShank = 1:numel(shankStruct.db)
      if shankStruct.db(dbCountShank).entryName == fnsData{dbCount}
        animal = shankStruct.db(dbCountShank).animal;
        series = shankStruct.db(dbCountShank).series;
        entryName = shankStruct.db(dbCountShank).entryName;
        break
      end
    end
    sites = shankStruct.sites;
    connectedSites = shankStruct.connectedSites;
    chanMap = shankStruct.chanMap;
    uCentres = shankStruct.uCentres;

% CREATE ELECTRODE MAP
    [eMap, hDim, vDim] = electrodeMap(probe, sites, connectedSites);

% ASSIGN ELECTRODES TO SEGMENTS AND OBTAIN COMPARISON PAIRS
    parts = [2 3 4 6 8 12];
    segCount = 0;
    pairs = [];
    for nP = 1:numel(parts)
      vDimP = ceil(vDim/parts(nP));
      pairs = [pairs; nchoosek(segCount+1:segCount+parts(nP),2)]; %#ok<AGROW>
      for s = parts(nP):-1:1
        segCount = segCount + 1;
        if vDimP*s > size(eMap,1)
          lowEnd = size(eMap,1);
        else
          lowEnd = vDimP*s;
        end
        if vDimP*(s-1)+1 > size(eMap,1)
          upEnd = size(eMap,1);
        else
          upEnd = vDimP*(s-1)+1;
        end
        segs{segCount} = sort(reshape(eMap(lowEnd:-1:upEnd,:), [1 (lowEnd-upEnd+1)*hDim]));
        segs{segCount} = segs{segCount}(~isnan(segs{segCount}));
      end
    end

% SAVE ELECTRODE MAPS AND LOOP OVER SEGMENT PAIRS
    assert(numel(load_opt.selectedUnits) == numel(units))
    assert(max(chanMap(:,2)) <= max(max(eMap)))
    dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.eMap = eMap;']; eval(dataString);
%     dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.eMap1 = eMap1;']; eval(dataString);
%     dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.eMap2 = eMap2;']; eval(dataString);
%     dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.eMap3 = eMap3;']; eval(dataString);
%     dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.eMap4 = eMap4;']; eval(dataString);
    save(dataFile,'dataStruct');
    parfor p = 1:size(pairs,1)
      fprintf('Started processing pair %i out of %i\n',p,size(pairs,1));
      
% CREATE SUMMARY DATA STRUCTURE
      qp = cell2struct({animal, series, entryName, sh, p}, ...
        {'animal', 'series', 'dbentryName', 'shank', 'pair'}, 2);
      
% SPIKES FOR THE PAIR SEGMENTS
      Tstart = 1;
      Tend = numel(spkAll);
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
      saveParfor(['s_' entryName '_' num2str(p) '.mat'], qp);
      fprintf('Finished processing pair %i out of %i\n',p,size(pairs,1));
      
    end % loop over segment pairs
    
% SAVE SHANK SUMMARY DATA
    fileList = dir(['s_' entryName '_*.mat']);
    if ~isempty(fileList)
      QP = struct([]);
      for i = 1:size(fileList,1)
        qp = load(fileList(i).name);
        QP{end+1} = qp.q; %#ok<SAGROW>
      end
      dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.QP = QP;'];
      eval(dataString);
      save(dataFile,'dataStruct');
      delete(['s_' entryName '_*.mat']);
    end
  end % loop over shanks
end % loop over db entries
delete(pool);
