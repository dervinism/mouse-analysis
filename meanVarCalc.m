function [m, mv1, mv2, mean_var, mean_var_timeBin] = meanVarCalc(data, sr)
% Mean and variance calculator. A helper function to AnPSD_units and
% eyeAnalysis.

Tstart = 1;
Tend = numel(data);
Tmid = round((Tend+Tstart)/2);

mv1 = zeros(14, 2);
mv2 = mv1;
m = min(mean(data(Tstart:Tmid)), mean(data(Tmid+1:Tend)));

for s = 1:14 % starting from 10ms bins multiplying by 2 each step
  mv1(s, :) = [mean(data(1:round(numel(data)/2))), var(data(1:round(numel(data)/2)))];
  mv2(s, :) = [mean(data(1+round(numel(data)/2):end)), var(data(1+round(numel(data)/2):end))];
  data = data(1:floor(length(data)/2)*2); data = sum(reshape(data, 2, []))';
end
mean_var = cat(3, mv1, mv2);
mean_var_timeBin = sr;