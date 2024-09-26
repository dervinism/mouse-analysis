function coh = spkCoh(phaseCoh, fRef)
% coh = spkCoh(phaseCoh, fRef)
%
% Function outputs coherence at a particular reference frequency for all units
% or MUAs.
% Input: phaseCoh - a structure variable holding the coherence analysis
%                   variables.
%        fRef - the frequency of interest.

coh = zeros(1,numel(phaseCoh));
for u = 1:numel(phaseCoh)
  freq = phaseCoh{u}.freq;
  freqInterp = unique([freq fRef]);
  if ~isfield(phaseCoh{u}, 'coh') || sum(isnan(phaseCoh{u}.coh)) == numel(phaseCoh{u}.coh)
    coh(u) = NaN;
    continue
  end
  cohFull = interp1(freq, phaseCoh{u}.coh, freqInterp, 'linear');
  [~, iFreq] = min(abs(freq - fRef));
  coh(u) = recentrePhase(cohFull(iFreq), 0);
end