function stPR = stprCalc(data1, data2, window)
% Spike triggered population rate (or crosscorrelation) calculator. A
% helper function to AnPSD_units, eyeAnalysis, and stpHalfCalc.

stPR = xcorr(data1', data2', window)/sum(data1);
stPR = stPR' / mean(stPR([1:10 numel(stPR-11):numel(stPR)]));