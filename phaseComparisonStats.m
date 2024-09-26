function [fPEst, fWTest, strPMethod, pEst, U2, pObs, U2Obs] = phaseComparisonStats(areas, freq, FOI, areaPhaseindividual, areasCritical)
% [fPEst, fWTest, strPMethod, pEst, U2, pObs, U2Obs] = phaseComparisonStats(areas, freq, FOI, areaPhaseindividual, areasCritical)
%
% Function performs Mardia-Watson-Wheeler Uniform-Scores and Watson U2
% tests comparing mean phase values between different conditions of
% vigilance but in the same area.

inds = ismember(freq,FOI);

% Mardia-Watson-Wheeler Uniform-Scores Test
fPEst = {};
fWTest = {};
strPMethod = {};
for iArea = 1:numel(areas)
  for f = 1:numel(freq)
    if inds(f) && ~isempty(areaPhaseindividual{1}{iArea}) && ~isempty(areaPhaseindividual{2}{iArea})
      cvfX{1} = areaPhaseindividual{1}{iArea}(:,f);
      cvfX{2} = areaPhaseindividual{2}{iArea}(:,f);
      [~, fPEst{iArea}(f), fWTest{iArea}(f), strPMethod{iArea}{f}] = mardiatestn_circ_equal(cvfX); %#ok<*AGROW>
    end
  end
end

% Watson U2 test
pEst = {};
U2 = {};
pObs = {};
U2Obs = {};
for iArea = 1:numel(areas)
  disp(['Processing data for ' areas{iArea} ' (area # ' num2str(iArea) '/' num2str(numel(areas)) ')']);
  if sum(ismember(areasCritical, areas{iArea}))
    for f = 1:numel(freq)
      if inds(f) && ~isempty(areaPhaseindividual{1}{iArea}) && ~isempty(areaPhaseindividual{2}{iArea})
        a1 = areaPhaseindividual{1}{iArea}(:,f);
        a2 = areaPhaseindividual{2}{iArea}(:,f);
        [pEst{iArea}(f), U2{iArea}(f)] = watsons_U2_approx_p(a1, a2);
        [pObs{iArea}(f), U2Obs{iArea}(f)] = watsons_U2_perm_test(a1, a2, 1000);
      end
    end
  end
end