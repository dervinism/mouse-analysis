function [rPhase, rhoPhase, rhoCircPhase, rCoherence, rhoCoherence, rhoCircCoherence] = corrStates(...
  figPhase, figCoherence, phaseVec, coherenceVec, cohConfVec)
% A helper function to compareHalves_figs, compareSeries_figs and
% globalFigs. It performs correlation analyses between phase values stored
% in a cell array and does the same for coherence. Then updates the
% correlation figures with coefficients and p-values and outputs them.

rPhase = zeros(2,numel(figPhase));
rhoPhase = zeros(2,numel(figPhase));
rhoCircPhase = zeros(2,numel(figPhase));
rCoherence = zeros(2,numel(figCoherence));
rhoCoherence = zeros(2,numel(figCoherence));
rhoCircCoherence = zeros(2,numel(figCoherence));
for i = 1:numel(figPhase)
  inds1 = isnan(phaseVec{i}(1,:));
  inds2 = isnan(phaseVec{i}(2,:));
  inds = ~(inds1 + inds2);
  if sum(inds)
    [R, P] = corr([phaseVec{i}(1,inds)' phaseVec{i}(2,inds)']);
    rPhase(1,i) = R(2);
    rPhase(2,i) = P(2);
    [R, P] = corr([phaseVec{i}(1,inds)' phaseVec{i}(2,inds)'], 'type','Spearman');
    rhoPhase(1,i) = R(2);
    rhoPhase(2,i) = P(2);
    [rhoCircPhase(1,i), rhoCircPhase(2,i)] = circ_corrcc(phaseVec{i}(1,inds), phaseVec{i}(2,inds));
  else
    rPhase(1,i) = NaN;
    rPhase(2,i) = NaN;
    rhoPhase(1,i) = NaN;
    rhoPhase(2,i) = NaN;
    rhoCircPhase(1,i) = NaN;
    rhoCircPhase(2,i) = NaN;
  end
  figure(figPhase(i));
  h1 = get(gca,'title');
  titleStr = get(h1,'string');
  title([titleStr '  r' num2str(rPhase(1,i)) ' p' num2str(rPhase(2,i))...
    ' rho' num2str(rhoPhase(1,i)) ' p' num2str(rhoPhase(2,i))...
    ' rhoCirc' num2str(rhoCircPhase(1,i)) ' p' num2str(rhoCircPhase(2,i))]);
  
  correctedCoh1 = coherenceVec{i}(1,:);
  correctedCoh1(correctedCoh1 - cohConfVec{i}(1,:) <= 0) = NaN;
  correctedCoh2 = coherenceVec{i}(2,:);
  correctedCoh2(correctedCoh2 - cohConfVec{i}(2,:) <= 0) = NaN;
  inds1 = isnan(correctedCoh1);
  inds2 = isnan(correctedCoh2);
  inds = ~(inds1 + inds2);
  if sum(inds)
    [R, P] = corr([correctedCoh1(inds)' correctedCoh2(inds)']);
    rCoherence(1,i) = R(2);
    rCoherence(2,i) = P(2);
    [R, P] = corr([correctedCoh1(inds)' correctedCoh2(inds)'], 'type','Spearman');
    rhoCoherence(1,i) = R(2);
    rhoCoherence(2,i) = P(2);
    [rhoCircCoherence(1,i), rhoCircCoherence(2,i)] = circ_corrcc(correctedCoh1(inds), correctedCoh2(inds));
  else
    rCoherence(1,i) = NaN;
    rCoherence(2,i) = NaN;
    rhoCoherence(1,i) = NaN;
    rhoCoherence(2,i) = NaN;
    rhoCircCoherence(1,i) = NaN;
    rhoCircCoherence(2,i) = NaN;
  end
  figure(figCoherence(i));
  h2 = get(gca,'title');
  titleStr = get(h2,'string');
  title([titleStr '  r' num2str(rCoherence(1,i)) ' p' num2str(rCoherence(2,i))...
    ' rho' num2str(rhoCoherence(1,i)) ' p' num2str(rhoCoherence(2,i))...
    ' rhoCirc' num2str(rhoCircCoherence(1,i)) ' p' num2str(rhoCircCoherence(2,i))]);
end