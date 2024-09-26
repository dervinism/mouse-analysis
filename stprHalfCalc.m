function stPR = stprHalfCalc(data1, data2, window)
% Spike triggered population rate (or crosscorrelation) calculator for data
% halves. A helper function to AnPSD_units, and eyeAnalysis.

Tstart = 1;
Tend = numel(data1);
Tmid = round((Tend+Tstart)/2);

stPR1 = stprCalc(data1(Tstart:Tmid), data2(Tstart:Tmid), window);
stPR2 = stprCalc(data1(Tmid+1:Tend), data2(Tmid+1:Tend), window);

stPR = [stPR1; stPR2];