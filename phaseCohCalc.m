function [freq, coh, phase, coh_conf, phase_confU, phase_confL, coh_halves_freq, coh_halves, coh_conf_halves, phase_halves,...
  phase_conf_halves] = phaseCohCalc(data1, data2, sr, opt)
% Phase and coherence calculator. A helper function to AnPSD_units and
% eyeAnalysis.

data1 = double(data1);
data2 = double(data2);
data1(isnan(data1)) = 0;
data2(isnan(data2)) = 0;
Tstart = 1;
Tend = floor(numel(data1)/2)*2;
Tmid = round(Tend/2);
if isfield(opt,'typespk1') && strcmpi(opt.typespk1,'c')
  data1_1sthalf = data1(Tstart:Tmid) - mean(data1(Tstart:Tmid), 'omitnan');
  data1_2ndhalf = data1(Tmid+1:Tend) - mean(data1(Tmid+1:Tend), 'omitnan');
  data1 = data1 - mean(data1, 'omitnan');
else
  data1_1sthalf = data1(Tstart:Tmid);
  data1_2ndhalf = data1(Tmid+1:Tend);
end
if isfield(opt,'typespk2') && strcmpi(opt.typespk2,'c')
  data2_1sthalf = data2(Tstart:Tmid) - mean(data2(Tstart:Tmid), 'omitnan');
  data2_2ndhalf = data2(Tmid+1:Tend) - mean(data2(Tmid+1:Tend), 'omitnan');
  data2 = data2 - mean(data2, 'omitnan');
else
  data2_1sthalf = data2(Tstart:Tmid);
  data2_2ndhalf = data2(Tmid+1:Tend);
end

% Half phase and coherence
if sum(data2_1sthalf) && sum(data1_1sthalf)
  [freq_1sthalf, coh_1sthalf, phase_1sthalf, coh_1sthalf_conf, phase_1sthalf_confU, phase_1sthalf_confL] = ...
    freqDependentWindowCoherenceMD(data1_1sthalf', data2_1sthalf', 1/sr, [], opt);
end
if sum(data2_2ndhalf) && sum(data1_2ndhalf)
  [freq_2ndhalf, coh_2ndhalf, phase_2ndhalf, coh_2ndhalf_conf, phase_2ndhalf_confU, phase_2ndhalf_confL] = ...
    freqDependentWindowCoherenceMD(data1_2ndhalf', data2_2ndhalf', 1/sr, [], opt);
  if ~sum(data1_1sthalf) || ~exist('freq_1sthalf','var')
    freq_1sthalf = NaN(size(freq_2ndhalf)); coh_1sthalf = NaN(size(coh_2ndhalf)); phase_1sthalf = NaN(size(phase_2ndhalf));
  end
else
  if sum(data1_1sthalf)
    if ~exist('freq_1sthalf','var')
      freq_1sthalf = NaN; coh_1sthalf = NaN; phase_1sthalf = NaN;
    end
    freq_2ndhalf = NaN(size(freq_1sthalf)); coh_2ndhalf = NaN(size(coh_1sthalf)); phase_2ndhalf = NaN(size(phase_1sthalf));
  else
    freq_1sthalf = NaN; coh_1sthalf = NaN; phase_1sthalf = NaN;
    freq_2ndhalf = NaN; coh_2ndhalf = NaN; phase_2ndhalf = NaN;
  end
end

% Full phase and coherence
if sum(data1) && sum(data2)
  [freq, coh, phase, coh_conf, phase_confU, phase_confL] = freqDependentWindowCoherenceMD(data1', data2', 1/sr, [], opt);
else
  freq = NaN; coh = NaN; phase = NaN; coh_conf = NaN; phase_confU = NaN; phase_confL = NaN;
end

% Remove phase values for which no confidence interval exists
if exist('phase_1sthalf','var') && exist('phase_1sthalf_confU','var') && exist('phase_1sthalf_confL','var')
  %phase_1sthalf(isnan(phase_1sthalf_confU)) = NaN;
  %phase_1sthalf(isnan(phase_1sthalf_confL)) = NaN;
else
  phase_1sthalf = NaN(size(freq_1sthalf));
  phase_1sthalf_confU = NaN(size(freq_1sthalf)); phase_1sthalf_confL = NaN(size(freq_1sthalf));
  coh_1sthalf = NaN(size(freq_1sthalf)); coh_1sthalf_conf = NaN(size(freq_1sthalf));
end
if exist('phase_2ndhalf','var') && exist('phase_2ndhalf_confU','var') && exist('phase_2ndhalf_confL','var')
  %phase_2ndhalf(isnan(phase_2ndhalf_confU)) = NaN;
  %phase_2ndhalf(isnan(phase_2ndhalf_confL)) = NaN;
else
  phase_2ndhalf = NaN(size(freq_2ndhalf));
  phase_2ndhalf_confU = NaN(size(freq_2ndhalf)); phase_2ndhalf_confL = NaN(size(freq_2ndhalf));
  coh_2ndhalf = NaN(size(freq_2ndhalf)); coh_2ndhalf_conf = NaN(size(freq_2ndhalf));
end
%phase(isnan(phase_confU)) = NaN;
%phase(isnan(phase_confL)) = NaN;

% Correct for misamtching frequencies (Chronux sometimes messes those up)
inds = ismember(freq_1sthalf, freq_2ndhalf);
freq_1sthalf = freq_1sthalf(:,inds);
coh_1sthalf = coh_1sthalf(:,inds);
phase_1sthalf = phase_1sthalf(:,inds);
coh_1sthalf_conf = coh_1sthalf_conf(:,inds);
phase_1sthalf_confU = phase_1sthalf_confU(:,inds);
phase_1sthalf_confL = phase_1sthalf_confL(:,inds);
  
inds = ismember(freq_2ndhalf, freq_1sthalf);
freq_2ndhalf = freq_2ndhalf(:,inds);
coh_2ndhalf = coh_2ndhalf(:,inds);
phase_2ndhalf = phase_2ndhalf(:,inds);
coh_2ndhalf_conf = coh_2ndhalf_conf(:,inds);
phase_2ndhalf_confU = phase_2ndhalf_confU(:,inds);
phase_2ndhalf_confL = phase_2ndhalf_confL(:,inds);

% Output
if exist('freq_1sthalf','var') && ~sum(isnan(freq_1sthalf))
  coh_halves_freq = freq_1sthalf;
else
  coh_halves_freq = freq_2ndhalf;
end
coh_halves = [torow(coh_1sthalf); torow(coh_2ndhalf)];
coh_conf_halves = [torow(coh_1sthalf_conf); torow(coh_2ndhalf_conf)];
phase_halves = [torow(phase_1sthalf); torow(phase_2ndhalf)];
phase_conf_halves = [torow(phase_1sthalf_confU); torow(phase_1sthalf_confL); torow(phase_2ndhalf_confU); torow(phase_2ndhalf_confL)];