function [phase, phaseConfU, phaseConfL, coh, cohConfU, cohConfL] = correctPhaseCoh(phase, phaseConfU, phaseConfL, coh, cohConf, rateadjust_kappa)
% [phase, phaseConfU, phaseConfL, coh, cohConfU, cohConfL] = correctPhaseCoh(phase, phaseConfU, phaseConfL, coh, cohConf, rateadjust_kappa)
%
% Function cleans up phase and coherence values so that any non-significant
% values and NaNs are dealt with and phase and coherence values are
% mutually compatible.

% Correct coherence
cohConfU = (coh + cohConf) .* rateadjust_kappa;
cohConfL = (coh - cohConf) .* rateadjust_kappa;
cohConfU(cohConfL <= 0) = NaN;
cohConfL(cohConfL <= 0) = NaN;
coh = coh .* rateadjust_kappa;
coh(isnan(cohConfU) | isnan(cohConfL)) = NaN;

% Correct phase
phase = bestUnwrap(phase);
phase(isnan(phaseConfU) | isnan(phaseConfL) | isnan(coh) | isnan(cohConfU) | isnan(cohConfL)) = NaN;

% Ascertain that significant number of entries is the same for both phase and coherence
coh(isnan(phase) | isnan(phaseConfU) | isnan(phaseConfL)) = NaN;
cohConfU(isnan(phase) | isnan(phaseConfU) | isnan(phaseConfL)) = NaN;
cohConfL(isnan(phase) | isnan(phaseConfU) | isnan(phaseConfL)) = NaN;
assert(sum(isnan(phase)) == sum(isnan(coh)));