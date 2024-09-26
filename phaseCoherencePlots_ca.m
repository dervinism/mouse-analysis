function phaseCoherencePlots_ca(sh, u, UPF, unitData, sr, unit, nUnits, x_lim, figsubdirname, entryName, figType, visibility, isPop, subpop)
% Draw phase and coherence figures. A helper function to
% AnPSD_units_figs_ca.

persistent figH uFig

if (mod(u, UPF) == 1 || isPop) || ~exist('figH', 'var') || ~isgraphics(figH)
  figH = figure('units', 'normalized', 'position', [0.002, .04, 1, .88], 'Visible', visibility);
end
if nUnits < UPF
  whichSubplot = u;
else
  whichSubplot = mod(u, UPF); if ~whichSubplot; whichSubplot=UPF; end
end

% Coherence with MUA on a different probe
subplot(3, min(UPF,nUnits), whichSubplot);
unwrapped = bestUnwrap(unitData.coh_halves(1,:));
if isempty(unwrapped)
  unwrapped = unitData.coh_halves(1,:);
end
hold on, plot(unitData.coh_halves_freq, unwrapped, '.-r');
unwrapped = bestUnwrap(unitData.coh_halves(2,:));
if isempty(unwrapped)
  unwrapped = unitData.coh_halves(2,:);
end
hold on, plot(unitData.coh_halves_freq, unwrapped, '.-m');
if whichSubplot == 1
  ylabel('Coherence')
end

% Phase
subplot(3, min(UPF,nUnits), whichSubplot + min(UPF,nUnits));
hold on, plot(unitData.coh_halves_freq, unitData.phase_halves(1,:), '.-r');
hold on, plot(unitData.coh_halves_freq, unitData.phase_halves(2,:), '.-m');
xlabel('Frequency (Hz)')
if whichSubplot == 1
  ylabel('Phase (rad)')
end

% Circular phase
if sum(~isnan(unitData.coh))
  figTmp = plotPhaseSpectrum(unitData.freq', unitData.coh, unitData.phase, unitData.phase_confU, unitData.phase_confL,...
    unitData.coh_conf, visibility);
end

% Rate-adjusted coherence
subplot(4,2,2)
hold on, semilogx(unitData.freq, unitData.coh .* unitData.rateadjust_kappa, 'b.-')
hold on, semilogx(unitData.freq, (unitData.coh + unitData.coh_conf) .* unitData.rateadjust_kappa, 'c--')
hold on, semilogx(unitData.freq, max(0,(unitData.coh - unitData.coh_conf) .* unitData.rateadjust_kappa), 'c--')
ylabel('coh, rate-adjusted coh')
if isPop
  title([entryName '_' num2str(sh) '_' unit], 'Interpreter', 'none')
else
  title([entryName '_' num2str(sh) '_u' num2str(unit)], 'Interpreter', 'none')
end
if strcmpi(figsubdirname,'')
  if isPop
    figName = [entryName '_' figType '_PSDpolar_sh' num2str(sh) '_' unit];
  else
    figName = [entryName '_' figType '_PSDpolar_sh' num2str(sh) '_u' num2str(unit)];
  end
else
  if isPop
    figName = [figsubdirname filesep entryName '_' figType '_PSDpolar_sh' num2str(sh) '_' unit];
  else
    figName = [figsubdirname filesep entryName '_' figType '_PSDpolar_sh' num2str(sh) '_u' num2str(unit)];
  end
end
if exist('figTmp','var')
  set(figTmp, 'Name',figName);
  hgsave(figTmp, figName);
  close(figTmp);
end

subplot(3, min(UPF,nUnits), whichSubplot);
set(gca, 'XScale', 'log');
set(gca, 'XTick', [0.03 0.1 0.3 1 3 10 30]);
xlim(x_lim);
ylim([0 1]);

subplot(3, min(UPF,nUnits), whichSubplot + min(UPF,nUnits));
set(gca, 'XScale', 'log');
set(gca, 'XTick', [0.03 0.1 0.3 1 3 10 30]);
xlim(x_lim);

% stPR
subplot(3, min(UPF,nUnits), whichSubplot + 2*min(UPF,nUnits));
hold on, plot(-100:100, unitData.stPR, 'k');
hold on, plot(-100:100, unitData.stPR_halves(1,:), 'b');
hold on, plot(-100:100, unitData.stPR_halves(2,:), 'b');
xlim(100*[-1 1]);
xlabel('Time (ms)')
if whichSubplot == 1
  ylabel('XCorr (1/mean)')
end
set(gca, 'XTick', -100:50:100);
set(gca, 'XTickLabel', {num2str(-100/sr), num2str(-50/sr), '0', num2str(50/sr), num2str(100/sr)})

if isPop
  uFig(whichSubplot) = 1;
else
  uFig(whichSubplot) = unit;
end
if u == nUnits || ~mod(u, UPF)
  suffix = '_u';
  for i = 1:whichSubplot
    suffix = [suffix sprintf('_%i',uFig(i))]; %#ok<AGROW>
  end
  uFig = [];
  if strcmpi(figsubdirname,'')
    figName = [entryName '_' figType '_PSD_sh' num2str(sh) '_p' num2str(ceil(u/UPF)) suffix];
  else
    figName = [figsubdirname filesep entryName '_' figType '_PSD_sh' num2str(sh) '_p' num2str(ceil(u/UPF)) suffix];
  end
  set(figH, 'Name',figName);
  hgsave(figH, figName);
  close(figH);
end

% Half phase and coherence comparisons
figTmp3 = figure('units', 'normalized', 'position', [0.002, .04, 1/2, .88/1.04], 'Visible', visibility);
subplot(2, 1, 1);
unitData.phase_halves(1,isnan(unitData.phase_conf_halves(1,:))) = NaN;
unitData.phase_halves(1,isnan(unitData.phase_conf_halves(2,:))) = NaN;
unitData.phase_halves(2,isnan(unitData.phase_conf_halves(3,:))) = NaN;
unitData.phase_halves(2,isnan(unitData.phase_conf_halves(4,:))) = NaN;
semilogx(unitData.coh_halves_freq, unitData.phase_halves(1,:), 'r');
hold on, semilogx(unitData.coh_halves_freq, unitData.phase_halves(2,:), 'm');
hold on, semilogx(unitData.coh_halves_freq, unitData.phase_conf_halves(1,:), '--r');
hold on, semilogx(unitData.coh_halves_freq, unitData.phase_conf_halves(2,:), '--r');
hold on, semilogx(unitData.coh_halves_freq, unitData.phase_conf_halves(3,:), '--m');
hold on, semilogx(unitData.coh_halves_freq, unitData.phase_conf_halves(4,:), '--m');
if isPop
  titleStr = [sprintf('Phase sh%i', sh)  '_' unit];
else
  titleStr = sprintf('Phase sh%i u%i', sh, unit);
end
title(titleStr, 'Interpreter', 'none')
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');
legend('1st half', '2nd half');
subplot(2, 1, 2);
semilogx(unitData.coh_halves_freq, unitData.coh_halves(1,:) .* unitData.rateadjust_kappa_halves(1,:), 'r');
hold on, semilogx(unitData.coh_halves_freq, unitData.coh_halves(2,:) .* unitData.rateadjust_kappa_halves(2,:), 'm');
hold on, semilogx(unitData.coh_halves_freq,...
  (unitData.coh_halves(1,:) + unitData.coh_conf_halves(1,:)) .* unitData.rateadjust_kappa_halves(1,:), '--r');
hold on, semilogx(unitData.coh_halves_freq,...
  (unitData.coh_halves(1,:) - unitData.coh_conf_halves(1,:)) .* unitData.rateadjust_kappa_halves(1,:), '--r');
hold on, semilogx(unitData.coh_halves_freq,...
  (unitData.coh_halves(2,:) + unitData.coh_conf_halves(2,:)) .* unitData.rateadjust_kappa_halves(2,:), '--m');
hold on, semilogx(unitData.coh_halves_freq,...
  (unitData.coh_halves(2,:) - unitData.coh_conf_halves(2,:)) .* unitData.rateadjust_kappa_halves(2,:), '--m');
if isPop
  titleStr = [sprintf('Coherence sh%i', sh)  '_' unit];
else
  titleStr = sprintf('Coherence sh%i u%i', sh, unit);
end
title(titleStr, 'Interpreter', 'none')
xlabel('Frequency (Hz)');
ylabel('Coherence');
legend('1st half', '2nd half');
if strcmpi(figsubdirname,'')
  if isPop
    figName = [entryName '_' figType '_half_phase_and_coherence_sh' num2str(sh) '_' unit];
  else
    figName = [entryName '_' figType '_half_phase_and_coherence_sh' num2str(sh) '_u' num2str(unit)];
  end
else
  if isPop
    figName = [figsubdirname filesep entryName '_' figType '_half_phase_and_coherence_sh' num2str(sh) '_' unit];
  else
    figName = [figsubdirname filesep entryName '_' figType '_half_phase_and_coherence_sh' num2str(sh) '_u' num2str(unit)];
  end
end
set(figTmp3, 'Name',figName);
hgsave(figTmp3, figName);
close(figTmp3);