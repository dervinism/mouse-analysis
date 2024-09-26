function [rMFR, rhoMFR] = corrStatesMFR(figMFR, mfrVEC)
% A helper function to globalFigs. It performs correlation analyses between
% half mean firing rate values stored in a cell array. Then it updates the
% correlation figures with coefficients and p-values and outputs them.

rMFR = zeros(2,numel(figMFR));
rhoMFR = zeros(2,numel(figMFR));
for i = 1:numel(figMFR)
  inds1 = isnan(mfrVEC(1,:));
  inds2 = isnan(mfrVEC(2,:));
  inds = ~(inds1 + inds2);
  if sum(inds)
    [R, P] = corr([mfrVEC(1,inds)' mfrVEC(2,inds)']);
    rMFR(1,i) = R(2);
    rMFR(2,i) = P(2);
    [R, P] = corr([mfrVEC(1,inds)' mfrVEC(2,inds)'], 'type','Spearman');
    rhoMFR(1,i) = R(2);
    rhoMFR(2,i) = P(2);
  else
    rMFR(1,i) = NaN;
    rMFR(2,i) = NaN;
    rhoMFR(1,i) = NaN;
    rhoMFR(2,i) = NaN;
  end
  figure(figMFR(i));
  h1 = get(gca,'title');
  titleStr = get(h1,'string');
  title([titleStr '  r' num2str(rMFR(1,i)) ' p' num2str(rMFR(2,i))...
    ' rho' num2str(rhoMFR(1,i)) ' p' num2str(rhoMFR(2,i))]);
end