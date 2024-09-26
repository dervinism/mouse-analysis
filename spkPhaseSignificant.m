function phaseDB = spkPhaseSignificant(phaseCoh, fRef)
% phaseDB = spkPhaseSignificant(phaseCoh, fRef)
%
% Function outputs phase at a particular reference frequency for all units
% or MUAs when it exists.
% Input: phaseCoh - a structure variable holding the coherence analysis
%                   variables.
%        fRef - the frequency of interest.

phaseDB = zeros(1,numel(phaseCoh));
for u = 1:numel(phaseCoh)
  if ~isfield(phaseCoh{u}, 'phase') || sum(isnan(phaseCoh{u}.phase)) == numel(phaseCoh{u}.phase)
    phaseDB(u) = NaN;
    continue
  end
  coh = phaseCoh{u}.coh;
  cohConf = phaseCoh{u}.coh_conf;
  rateadjust_kappa = phaseCoh{u}.rateadjust_kappa;
  phase = phaseCoh{u}.phase;
  phaseConfU = phaseCoh{u}.phase_confU;
  phaseConfL = phaseCoh{u}.phase_confL;
  phase = correctPhaseCoh(phase, phaseConfU, phaseConfL, coh, cohConf, rateadjust_kappa);
  
  freq = phaseCoh{u}.freq;
  freqInterp = unique([freq fRef]);
  
  phaseFull = interp1(freq, phase, freqInterp, 'linear');
  [~, iFreq] = min(abs(freq - fRef));
  phaseDB(u) = recentrePhase(phaseFull(iFreq), 0);
end