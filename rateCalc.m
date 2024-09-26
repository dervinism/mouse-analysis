function [mfr, mfr_1sthalf, mfr_2ndhalf, lfr1, lfr5] = rateCalc(data, sr)
% A spiking rate calculator. A helper function to AnPSD_units and
% eyeAnalysis.

Tstart = 1;
Tend = numel(data);
Tmid = round((Tend+Tstart)/2);

% Mean firing rate
mfr = mean(data)*sr;
mfr_1sthalf = mean(data(Tstart:Tmid))*sr;
mfr_2ndhalf = mean(data(Tmid+1:Tend))*sr;

% Longitudinal firing rate
lfr1 = resampleData(data, sr, 1/60); % per 60 s
lfr5 = resampleData(data, sr, 1/60/5); % per 5 min