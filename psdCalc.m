function [freq_1sthalf, psd_halves, freq, psd, psd_conf, psd_numelSignal, beta_halves, beta] = psdCalc(data, sr, opt)
% A power spectral density calculator. A helper function to AnPSD_units and
% eyeAnalysis.

Tstart = 1;
Tend = numel(data);
Tmid = round((Tend+Tstart)/2);

% Half PSD
[freq_1sthalf, psd_1sthalf] = freqDependentWindowCoherenceMD(data(Tstart:Tmid), [], 1/sr, [], opt);
[freq_2ndhalf, psd_2ndhalf] = freqDependentWindowCoherenceMD(data(Tmid+1:Tend), [], 1/sr, [], opt);
beta_1sthalf = psdBeta(freq_1sthalf, psd_1sthalf);
beta_2ndhalf = psdBeta(freq_2ndhalf, psd_2ndhalf);

% Full PSD
[freq, psd, ~, psd_conf] = freqDependentWindowCoherenceMD(data(Tstart:Tend), [], 1/sr, [], opt);
beta = psdBeta(freq, psd);

% Output
psd_halves = [psd_1sthalf; psd_2ndhalf];
beta_halves = [beta_1sthalf; beta_2ndhalf];
assert(numel(freq_1sthalf) == numel(freq_2ndhalf) && max(abs(freq_1sthalf - freq_2ndhalf)) < 1e-9)
psd_numelSignal = Tend - Tstart + 1;