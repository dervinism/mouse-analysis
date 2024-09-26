function [psd_halves_freq, psd_halves, freq, psd, psd_conf, coh, phase, coh_conf, phase_confU, phase_confL, coh_halves_freq, coh_halves,...
  coh_conf_halves, phase_halves, phase_conf_halves, rateadjust_kappa, rateadjust_kappa_halves] = correctFreq(psd_halves_freq, psd_halves, freqPSD,...
  psd, psd_conf, freq, coh, phase, coh_conf, phase_confU, phase_confL, coh_halves_freq, coh_halves, coh_conf_halves, phase_halves,...
  phase_conf_halves, rateadjust_kappa, rateadjust_kappa_halves)
% Sometimes Chronux outputs different frequencies for two non-identical
% signals that have the same length following coherence analysis. This
% function corrects for missing or non-matching frequencies. 

if numel(coh) == 1 && isnan(coh)
  rateadjust_kappa = 1;
  rateadjust_kappa_halves = [1;1];
else
  
  inds = ismember(freqPSD, freq);
  freqPSD = freqPSD(inds);
  psd = psd(inds);
  psd_conf = psd_conf(inds);
  rateadjust_kappa = rateadjust_kappa(inds);
  
  inds = ismember(freq, freqPSD);
  freq = freq(inds);
  coh = coh(inds);
  phase = phase(inds);
  coh_conf = coh_conf(:,inds);
  phase_confU = phase_confU(:,inds);
  phase_confL = phase_confL(:,inds);
  
  inds = ismember(psd_halves_freq, coh_halves_freq);
  psd_halves_freq = psd_halves_freq(:,inds);
  psd_halves = psd_halves(:,inds);
  rateadjust_kappa_halves = rateadjust_kappa_halves(:,inds);
  
  inds = ismember(coh_halves_freq, psd_halves_freq);
  coh_halves_freq = coh_halves_freq(:,inds);
  coh_halves = coh_halves(:,inds);
  phase_halves = phase_halves(:,inds);
  coh_conf_halves = coh_conf_halves(:,inds);
  phase_conf_halves = phase_conf_halves(:,inds);
end