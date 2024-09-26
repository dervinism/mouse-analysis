function [breakClause, mfr] = firingRateTest(raster, sr)
% [breakClause, mfr] = firingRateTest(raster, sr)
%
% Function checks if the supplied firing rate vector has a non-zero firing
% rate. sr is sampling frequency in Hz. If the firing rate is zero a true
% breakClause is issued. Mean firing rate is also outputed.

breakClause = false;
if sum(raster == 0) == numel(raster) || sum(isnan(raster)) || isempty(raster)
  breakClause = true;
  mfr = 0;
  return
end
mfr = rateCalc(raster, sr);
if mfr == 0 || isnan(mfr) || isempty(mfr)
  breakClause = true;
end