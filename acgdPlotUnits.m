function [fH, output] = acgdPlotUnits(lags, data, areas, condition, srNew, fitEndLag, yLim)


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
end
downsampledLags = downsampledLags - (1/srNew)/2;
fitEndInd = find(fitEndLag >= downsampledLags, 1, 'last');


%% Fit the data
adaptationEnd = cell(nAreas,1);
decayFitMdl = cell(nAreas,1);
%decayFit = cell(nAreas,1);
effectiveTau = cell(nAreas,1);
effectiveTauMean = NaN(1,nAreas);
effectiveTauCI95 = NaN(2,nAreas);
fittedTau = cell(nAreas,1);
fittedTauMean = NaN(1,nAreas);
fittedTauCI95 = NaN(2,nAreas);
A = cell(nAreas,1);
B = cell(nAreas,1);
modelfun = @(A,downsampledLags) A(1) * (exp(-downsampledLags/A(2)) + A(3));
for iArea = 1:nAreas
  areaCode = determineArea(areas{iArea});
  for u = 1:size(data{condition}{areaCode},1)
    if ~sum(~data{condition}{areaCode}(u,:))
      [A1Guessed, adaptationEnd{iArea}(u)] = max(data{condition}{areaCode}(u,:));
      A3Guessed = mean(data{condition}{areaCode}(u,end-5:end))/A1Guessed;
      effectiveDecay = (A1Guessed-A3Guessed*A1Guessed)*(1/exp(1));
      effectiveTauInd = find(data{condition}{areaCode}(u,adaptationEnd{iArea}(u):end) <= effectiveDecay+A3Guessed*A1Guessed, 1);
      if ~isempty(effectiveTauInd)
        effectiveTau{iArea}(u) = downsampledLags(effectiveTauInd)-downsampledLags(1);
        A2Guessed = effectiveTau{iArea}(u);
        A0 = [A1Guessed, A2Guessed, A3Guessed];
        try
          decayFitMdl{iArea}{u} = fitnlm(table(downsampledLags(adaptationEnd{iArea}(u):fitEndInd)',data{condition}{areaCode}(u,adaptationEnd{iArea}(u):fitEndInd)'), modelfun, A0);
          coefficients = decayFitMdl{iArea}{u}.Coefficients{:, 'Estimate'};
          %decayFit{iArea} = [decayFit{iArea}; coefficients(1) .* (exp(-lags./coefficients(2)) + coefficients(3))];
          fittedTau{iArea}(u) = coefficients(2);
          A{iArea}(u) = coefficients(1);
          B{iArea}(u) = coefficients(3);
        catch
          decayFitMdl{iArea}{u} = NaN;
          %decayFit{iArea} = [decayFit{iArea}; NaN(size(lags))];
          fittedTau{iArea}(u) = NaN;
          A{iArea}(u) = NaN;
          B{iArea}(u) = NaN;
        end
%         figure; plot(lags,data{condition}{areaCode}(u,:)); hold on
%         plot(lags(adaptationEnd{iArea}(u)),A1Guessed, 'g.', 'MarkerSize',10);
%         plot(lags(adaptationEnd{iArea}(u)+effectiveTauInd-1),data{condition}{areaCode}(u,adaptationEnd{iArea}(u)+effectiveTauInd-1), 'r.', 'MarkerSize',10);
%         plot(lags([1 end]),[A3Guessed*A1Guessed A3Guessed*A1Guessed]); hold off
      else
        effectiveTau{iArea}(u) = NaN;
        decayFitMdl{iArea}{u} = NaN;
        %decayFit{iArea} = [decayFit{iArea}; NaN(size(lags))];
        fittedTau{iArea}(u) = NaN;
        A{iArea}(u) = NaN;
        B{iArea}(u) = NaN;
      end
    else
      adaptationEnd{iArea}(u) = NaN;
      effectiveTau{iArea}(u) = NaN;
      decayFitMdl{iArea}{u} = NaN;
      %decayFit{iArea} = [decayFit{iArea}; NaN(size(lags))];
      fittedTau{iArea}(u) = NaN;
      A{iArea}(u) = NaN;
      B{iArea}(u) = NaN;
    end
  end
  [effectiveTauMean(iArea), effectiveTauCI95(:,iArea)] = datamean(effectiveTau{iArea}');
  [fittedTauMean(iArea), fittedTauCI95(:,iArea)] = datamean(fittedTau{iArea}');
end


%% Draw the violins
if isempty(yLim)
  options.yLim = [0 1];
else
  options.yLim = yLim;
end
options.yLabel = '\tau_{effective} (s)';
[statsEffective.pval_ttest, statsEffective.area1, statsEffective.area2] = ttestGroup(areas, effectiveTau);
fH1 = multiViolinPlots(effectiveTau, areas, effectiveTauMean, effectiveTauCI95, statsEffective, options);

if isempty(yLim)
  options.yLim = [0 3];
else
  options.yLim = yLim;
end
options.yLabel = '\tau (s)';
[statsFitted.pval_ttest, statsFitted.area1, statsFitted.area2] = ttestGroup(areas, fittedTau);
fH2 = multiViolinPlots(fittedTau, areas, fittedTauMean, fittedTauCI95, statsFitted, options);


%% Assign output
fH = [fH1 fH2];
%output.decayFitMdl = decayFitMdl;
%output.decayFit = decayFit;
output.effectiveTau = effectiveTau;
output.effectiveTauMean = effectiveTauMean;
output.effectiveTauCI95 = effectiveTauCI95;
output.fittedTau = fittedTau;
output.fittedTauMean = fittedTauMean;
output.fittedTauCI95 = fittedTauCI95;
output.A = A;
output.B = B;
output.lags = downsampledLags;
output.srLag = srNew;
output.fitEndLag = fitEndLag;
output.modelfun = modelfun;
output.areas = areas;
output.condition = condition;



%% Local functions
function [p, area1, area2] = ttestGroup(areaNames, data)

areaCombos = nchoosek(1:numel(areaNames), 2);
nCombos = size(areaCombos,1);
area1 = areaNames(areaCombos(:,1));
area2 = areaNames(areaCombos(:,2));
p = zeros(1,nCombos);
for iCombo = 1:nCombos
  areaCode1 = areaCombos(iCombo,1);
  areaCode2 = areaCombos(iCombo,2);
  [~,p(iCombo)] = ttest2(data{areaCode1}, data{areaCode2});
end