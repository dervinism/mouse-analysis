function [results, state] = combineEyeData(ifilename, ofilename)
% Combines two or more pupil data analysis files.

for iFile = 1:numel(ifilename)
  data{iFile} = load(ifilename{iFile}); %#ok<AGROW>
end

results = data{1}.results;
for iFile = 2:numel(ifilename)
  results.x = [results.x; data{iFile}.results.x];
  results.y = [results.y; data{iFile}.results.y];
  results.aAxis = [results.aAxis; data{iFile}.results.aAxis];
  results.bAxis = [results.bAxis; data{iFile}.results.bAxis];
  results.abAxis = [results.abAxis; data{iFile}.results.abAxis];
  results.area = [results.area; data{iFile}.results.area];
  results.goodFit = [results.goodFit; data{iFile}.results.goodFit];
  results.blink = [results.blink; data{iFile}.results.blink];
  results.saturation = [results.saturation; data{iFile}.results.saturation];
  results.threshold = [results.threshold; data{iFile}.results.threshold];
  results.roi = [results.roi; data{iFile}.results.roi];
  results.equation = [results.equation; data{iFile}.results.equation];
  results.xxContour = [results.xxContour; data{iFile}.results.xxContour];
  results.yyContour = [results.yyContour; data{iFile}.results.yyContour];
  results.blinkRho = [results.blinkRho; data{iFile}.results.blinkRho];
end

state = data{1}.state;

save(ofilename, 'results', 'state');