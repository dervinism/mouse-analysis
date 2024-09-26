function phase = spkPhase(phaseCoh, fRef)
% phase = spkPhase(phaseCoh, fRef)
%
% Function outputs phase at a particular reference frequency for all units
% or MUAs.
% Input: phaseCoh - a structure variable holding the coherence analysis
%                   variables.
%        fRef - the frequency of interest.

phase = zeros(1,numel(phaseCoh));
for u = 1:numel(phaseCoh)
  freq = phaseCoh{u}.freq;
  freqInterp = unique([freq fRef]);
  if ~isfield(phaseCoh{u}, 'phase') || sum(isnan(phaseCoh{u}.phase)) == numel(phaseCoh{u}.phase)
    phase(u) = NaN;
    continue
  end
  phaseFull = interp1(freq, phaseCoh{u}.phase, freqInterp, 'linear');
  [~, iFreq] = min(abs(freq - fRef));
  phase(u) = recentrePhase(phaseFull(iFreq), 0);
end