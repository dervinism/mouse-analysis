function [fPEst, fWTest, strPMethod, pEst, U2, pObs, U2Obs, comparisons] = phaseComparisonStats_ca(areas, areasLocal, freq, freqLocal, FOI, areaPhaseIndividual, areaPhaseIndividualLocal, areasCritical, areasCriticalLocal)
% [fPEst, fWTest, strPMethod, pEst, U2, pObs, U2Obs] = phaseComparisonStats_ca(areas, freq, freqLocal, FOI, areaPhaseindividual, areaPhaseindividualLocal, areasCritical)
%
% Function performs Mardia-Watson-Wheeler Uniform-Scores and Watson U2
% tests comparing mean phase values between different areas within the same
% condition of vigilance.

assert(numel(freq) == numel(freqLocal));

inds = ismember(freq,FOI);
indsLocal = ismember(freqLocal,FOI);
assert(sum(inds) == sum(logical(inds+indsLocal)));

% Mardia-Watson-Wheeler Uniform-Scores Test: cross-area vs cross-area
fPEst = {};
fWTest = {};
strPMethod = {};
comparison1 = {};
comparison2 = {};
for iCond = 1:min([2 numel(areaPhaseIndividual)])
  for iArea = 1:numel(areas)
    for iArea2 = 1:numel(areas)
      fCount = 0;
      for f = 1:numel(freq)
        disp(['Processing data for condition ' num2str(iCond) ' ' areas{iArea} ' v ' areas{iArea2} ' (area # ' num2str(iArea) '/' num2str(numel(areas)+numel(areasLocal)) ')']);
        if inds(f) && ~isempty(areaPhaseIndividual{iCond}{iArea}) && ~isempty(areaPhaseIndividual{iCond}{iArea2})
          fCount = fCount + 1;
          cvfX{1} = areaPhaseIndividual{iCond}{iArea}(:,f);
          cvfX{2} = areaPhaseIndividual{iCond}{iArea2}(:,f);
          [~, fPEst{iCond}{iArea}{iArea2}(fCount), fWTest{iCond}{iArea}{iArea2}(fCount), strPMethod{iCond}{iArea}{iArea2}{fCount}] = mardiatestn_circ_equal(cvfX); %#ok<*AGROW>
          comparison1{iCond}{iArea} = areas{iArea};
          comparison2{iCond}{iArea}{iArea2} = areas{iArea2};
        end
      end
    end
  end
end

% Mardia-Watson-Wheeler Uniform-Scores Test: local vs cross-area
iAreaInit = iArea;
for iCond = 1:min([2 numel(areaPhaseIndividual)])
  for iArea = 1:numel(areasLocal)
    for iArea2 = 1:numel(areas)
      fCount = 0;
      for f = 1:numel(freq)
        disp(['Processing data for condition ' num2str(iCond) ' ' areasLocal{iArea} ' v ' areas{iArea2} ' (area # ' num2str(iAreaInit+iArea) '/' num2str(numel(areas)+numel(areasLocal)) ')']);
        if inds(f) && ~isempty(areaPhaseIndividualLocal{iCond}{iArea}) && ~isempty(areaPhaseIndividual{iCond}{iArea2})
          fCount = fCount + 1;
          cvfX{1} = areaPhaseIndividualLocal{iCond}{iArea}(:,f);
          cvfX{2} = areaPhaseIndividual{iCond}{iArea2}(:,f);
          [~, fPEst{iCond}{iAreaInit+iArea}{iArea2}(fCount), fWTest{iCond}{iAreaInit+iArea}{iArea2}(fCount), strPMethod{iCond}{iAreaInit+iArea}{iArea2}{fCount}] = mardiatestn_circ_equal(cvfX); %#ok<*AGROW>
          comparison1{iCond}{iAreaInit+iArea} = areasLocal{iArea};
          comparison2{iCond}{iAreaInit+iArea}{iArea2} = areas{iArea2};
        end
      end
    end
  end
end

comparisons.comparison1 = comparison1;
comparisons.comparison2 = comparison2;

% Watson U2 test: cross-area vs cross-area
pEst = {};
U2 = {};
pObs = {};
U2Obs = {};
for iCond = 1:min([2 numel(areaPhaseIndividual)])
  pEstCond = {};
  U2Cond = {};
  pObsCond = {};
  U2ObsCond = {};
  parfor iArea = 1:numel(areas)
    for iArea2 = 1:numel(areas)
      fCount = 0;
      if sum(ismember(areasCritical, areas{iArea})) && sum(ismember(areasCritical, areas{iArea2}))
        for f = 1:numel(freq)
          disp(['Processing data for condition ' num2str(iCond) ' ' areas{iArea} ' v ' areas{iArea2} ' (area # ' num2str(iArea) '/' num2str(numel(areas)+numel(areasLocal)) ')']);
          if inds(f) && ~isempty(areaPhaseIndividual{iCond}{iArea}) && ~isempty(areaPhaseIndividual{iCond}{iArea2})
            fCount = fCount + 1;
            a1 = areaPhaseIndividual{iCond}{iArea}(:,f);
            a2 = areaPhaseIndividual{iCond}{iArea2}(:,f);
            [pEstCond{iArea}{iArea2}(fCount), U2Cond{iArea}{iArea2}(fCount)] = watsons_U2_approx_p(a1, a2);
            [pObsCond{iArea}{iArea2}(fCount), U2ObsCond{iArea}{iArea2}(fCount)] = watsons_U2_perm_test(a1, a2, 1000);
          end
        end
      end
    end
  end
  pEst{iCond} = pEstCond;
  U2{iCond} = U2Cond;
  pObs{iCond} = pObsCond;
  U2Obs{iCond} = U2ObsCond;
end

% Watson U2 test: local vs cross-area
iAreaInit = numel(areas);
for iCond = 1:min([2 numel(areaPhaseIndividual)])
  pEstCond = pEst{iCond};
  U2Cond = U2{iCond};
  pObsCond = pObs{iCond};
  U2ObsCond = U2Obs{iCond};
  parfor iArea = 1:numel(areasLocal)
    for iArea2 = 1:numel(areas)
      fCount = 0;
      if sum(ismember(areasCriticalLocal, areasLocal{iArea})) && sum(ismember(areasCritical, areas{iArea2}))
        for f = 1:numel(freq)
          disp(['Processing data for condition ' num2str(iCond) ' ' areasLocal{iArea} ' v ' areas{iArea2} ' (area # ' num2str(iAreaInit+iArea) '/' num2str(numel(areas)+numel(areasLocal)) ')']);
          if inds(f) && ~isempty(areaPhaseIndividualLocal{iCond}{iArea}) && ~isempty(areaPhaseIndividual{iCond}{iArea2}) %#ok<*PFBNS>
            fCount = fCount + 1;
            a1 = areaPhaseIndividualLocal{iCond}{iArea}(:,f);
            a2 = areaPhaseIndividual{iCond}{iArea2}(:,f);
            [pEstCond{iAreaInit+iArea}{iArea2}(fCount), U2Cond{iAreaInit+iArea}{iArea2}(fCount)] = watsons_U2_approx_p(a1, a2);
            [pObsCond{iAreaInit+iArea}{iArea2}(fCount), U2ObsCond{iAreaInit+iArea}{iArea2}(fCount)] = watsons_U2_perm_test(a1, a2, 1000);
          end
        end
      end
    end
  end
  pEst{iCond} = pEstCond;
  U2{iCond} = U2Cond;
  pObs{iCond} = pObsCond;
  U2Obs{iCond} = U2ObsCond;
end