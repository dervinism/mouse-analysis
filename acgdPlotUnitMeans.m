function [fH, output] = acgdPlotUnitMeans(lags, data, areas, condition, srNew, fitEndLag, yLim)


%% Parse user input
if nargin < 7
  yLim = [];
end
if nargin < 6
  fitEndLag = lags(end);
end
srOriginal = 1/mean(lags(2:end)-lags(1:end-1));
if nargin < 5
  srNew = srOriginal;
end
if nargin < 4
  condition = 1;
end
nAreas = numel(areas);


%% Downsample the data
if srNew < srOriginal
  for iArea = 1:nAreas
    areaCode = determineArea(areas{iArea});
    [data{condition}{areaCode}, downsampledLags] = downsampleRasterMatrix(data{condition}{areaCode}, srOriginal, srNew);
    data{condition}{areaCode} = data{condition}{areaCode}./(srOriginal/srNew);
  end
elseif srNew > srOriginal
  error('Data upsampling is not supported');
else
  downsampledLags = (1:numel(data{condition}{areaCode}))./srOriginal;
end
downsampledLags = downsampledLags - (1/srNew)/2;
fitEndInd = find(fitEndLag >= downsampledLags, 1, 'last');


%% Fit the data
decayMean = cell(nAreas,1);
decayCI95 = cell(nAreas,1);
adaptationEnd = cell(nAreas,1);
decayFitMdl = cell(nAreas,1);
decayFit = cell(nAreas,1);
effectiveTau = cell(nAreas,1);
fittedTau = cell(nAreas,1);
modelfun = @(A,downsampledLags) A(1) * (exp(-downsampledLags/A(2)) + A(3));
for iArea = 1:nAreas
  areaCode = determineArea(areas{iArea});
  [decayMean{iArea}, decayCI95{iArea}] = datamean(data{condition}{areaCode});
  [A1Guessed, adaptationEnd{iArea}] = max(decayMean{iArea});
  A3Guessed = decayMean{iArea}(end)/A1Guessed;
  effectiveDecay = (A1Guessed-A3Guessed*A1Guessed)*(1/exp(1));
  effectiveTauInd = find(decayMean{iArea}(adaptationEnd{iArea}:end) <= effectiveDecay+A3Guessed*A1Guessed, 1);
  effectiveTau{iArea} = downsampledLags(effectiveTauInd)-downsampledLags(1);
  A2Guessed = effectiveTau{iArea};
  A0 = [A1Guessed, A2Guessed, A3Guessed];
  decayFitMdl{iArea} = fitnlm(table(downsampledLags(adaptationEnd{iArea}:fitEndInd)',decayMean{iArea}(adaptationEnd{iArea}:fitEndInd)'), modelfun, A0);
  coefficients = decayFitMdl{iArea}.Coefficients{:, 'Estimate'};
  decayFit{iArea} = coefficients(1) .* (exp(-lags./coefficients(2)) + coefficients(3));
  fittedTau{iArea} = coefficients(2);
end


%% Draw the fits
fH = figProperties('Mean unit auto-correlation decay fit', 'normalized', [0, .005, .97, .90], 'w', 'on');
hold on

pH = [];
for iArea = 1:nAreas
  errorbar(downsampledLags, decayMean{iArea}, abs(decayCI95{iArea}(1,:)), abs(decayCI95{iArea}(2,:)), '.',...
    'MarkerSize',10, 'Color',areaColours2(areas{iArea}));
  p = plot(lags,decayFit{iArea}, 'Color',areaColours2(areas{iArea}));
  pH = [pH p];
end
legend(pH, areas);
legend boxoff


%% Graph adjustments
if isempty(yLim)
  yLim = ylim;
end
axesProperties({}, 1, 'normal', 'off', 'w', 'Arial', 25, 4/3, 2, [0.003 0.003], 'out',...
  'on', 'k', {'Lag (s)'}, xlim, xticks,...
  'on', 'k', {'Auto-correlation'}, yLim, yticks);

xLim = xlim;
xAxisLength = xLim(2)-xLim(1);
yAxisLength = yLim(2)-yLim(1);
textStr = [];
for iArea = 1:numel(areas)
  textStr = [textStr '\tau_{' areas{iArea} '}=' num2str(fittedTau{iArea}) '   ']; %#ok<*AGROW>
end
text(xLim(2)-xAxisLength*0.95, yLim(2)-yAxisLength*0.05, textStr, 'FontSize',20);

textStr = [];
for iArea = 1:numel(areas)
  textStr = [textStr '\tau_{eff.' areas{iArea} '}=' num2str(effectiveTau{iArea}) '   ']; %#ok<*AGROW>
end
text(xLim(2)-xAxisLength*0.95, yLim(2)-yAxisLength*0.15, textStr, 'FontSize',20);


%% Assign output
output.decayFitMdl = decayFitMdl;
output.decayFit = decayFit;
output.fittedTau = fittedTau;
output.effectiveTau = effectiveTau;
output.decayMean = decayMean;
output.decayCI95 = decayCI95;
output.lags = downsampledLags;
output.srLag = srNew;
output.fitEndLag = fitEndLag;
output.modelfun = modelfun;
output.areas = areas;
output.condition = condition;