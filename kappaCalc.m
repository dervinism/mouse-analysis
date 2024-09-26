function [rateadjust_kappa, rateadjust_kappa_halves] = kappaCalc(mfr, mfr_1sthalf, mfr_2ndhalf, psd, psd_halves, coh, coh_halves)
% It's a helper function to AnPSD_units and eyeAnalysis scripts.
%
% Compute the rate adjustment factors for coherence report, per Aoi et al.
% we adjust as if the neuron's rate is 1 spk/s

rateadjust_kappa = coherenceAdjustment(mfr, 1, psd, 0, 1);

assert(isreal(psd) && isreal(coh) && isreal(rateadjust_kappa))
rateadjust_kappa_1sthalf = coherenceAdjustment(mfr_1sthalf, 1, psd_halves(1,:), 0, 1);

assert(isreal(psd_halves(1,:)) && isreal(coh_halves(1,:)) && isreal(rateadjust_kappa_1sthalf))
rateadjust_kappa_2ndhalf = coherenceAdjustment(mfr_2ndhalf, 1, psd_halves(2,:), 0, 1);

assert(isreal(psd_halves(2,:)) && isreal(coh_halves(2,:)) && isreal(rateadjust_kappa_2ndhalf))
rateadjust_kappa_halves = [rateadjust_kappa_1sthalf; rateadjust_kappa_2ndhalf];