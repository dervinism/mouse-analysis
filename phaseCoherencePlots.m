function phaseCoherencePlots(sh, u, UPF, unitData, sr, unit, nUnits, powerUnits, x_lim, figsubdirname, entryName, figType,...
  visibility, spkAllShanks, shank, spkAll, isPop)
% Draw phase and coherence figures. A helper function to AnPSD_units_figs
% and eyeAnalysis_figs.

persistent figH uFig

if nargin < 17
  isPop = false;
end

if (mod(u, UPF) == 1 || isPop) || ~exist('figH', 'var') || isempty(figH) || ~isgraphics(figH)
  figH = figure('units', 'normalized', 'position', [0.002, .04, 1, .88], 'Visible', visibility);
end
if nUnits < UPF
  whichSubplot = u;
else
  whichSubplot = mod(u, UPF); if ~whichSubplot; whichSubplot=UPF; end
end

% PSD
subplot(5, min(UPF,nUnits), whichSubplot);
hold on, plot(unitData.psd_halves_freq, unitData.psd_halves(1,:), 'r');
hold on, plot(unitData.psd_halves_freq, unitData.psd_halves(2,:), 'm');
if isPop
  titleStr = sprintf('sh%i_mfr1_%3.1f_mfr2_%3.1f', sh, unitData.mfr_1sthalf, unitData.mfr_2ndhalf);
else
  titleStr = sprintf('sh%i_u%i_mfr1_%3.1f_mfr2_%3.1f', sh, unit, unitData.mfr_1sthalf, unitData.mfr_2ndhalf);
end
title(titleStr, 'Interpreter', 'none')
if whichSubplot == 1
  ylabel(['Power (' powerUnits ')'])
end

% Coherence with MUA
subplot(5, min(UPF,nUnits), whichSubplot + min(UPF,nUnits));
if numel(spkAllShanks) > 0 && numel(shank) > 0 && ~isempty(spkAll)
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
  subplot(5, min(UPF,nUnits), whichSubplot + 2*min(UPF,nUnits));
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
    title([entryName '_' num2str(sh)], 'Interpreter', 'none')
  else
    title([entryName '_' num2str(sh) '_u' num2str(unit)], 'Interpreter', 'none')
  end
  if strcmpi(figsubdirname,'')
    if isPop
      figName = [entryName '_' figType '_PSDpolar_sh' num2str(sh)];
    else
      figName = [entryName '_' figType '_PSDpolar_sh' num2str(sh) '_u' num2str(unit)];
    end
  else
    if isPop
      figName = [figsubdirname filesep entryName '_' figType '_PSDpolar_sh' num2str(sh)];
    else
      figName = [figsubdirname filesep entryName '_' figType '_PSDpolar_sh' num2str(sh) '_u' num2str(unit)];
    end
  end
  if exist('figTmp','var')
    set(figTmp, 'Name',figName);
    hgsave(figTmp, figName);
    close(figTmp);
  end
else
  error('Data is empty or was not processed or loaded properly.')
end

subplot(5, min(UPF,nUnits), whichSubplot);
set(gca, 'XScale', 'log');
set(gca, 'XTick', [0.03 0.1 0.3 1 3 10 30]);
xlim(x_lim);

subplot(5, min(UPF,nUnits), whichSubplot + min(UPF,nUnits));
set(gca, 'XScale', 'log');
set(gca, 'XTick', [0.03 0.1 0.3 1 3 10 30]);
xlim(x_lim);
ylim([0 1]);

subplot(5, min(UPF,nUnits), whichSubplot + 2*min(UPF,nUnits));
set(gca, 'XScale', 'log');
set(gca, 'XTick', [0.03 0.1 0.3 1 3 10 30]);
xlim(x_lim);

% Mean vs var using different bin sizes
subplot(5, min(UPF,nUnits), whichSubplot + 3*min(UPF,nUnits));
hold on, plot(log(unitData.mv1(:,1)), log(unitData.mv1(:,2)) , '*-r');
hold on, plot(log(unitData.mv2(:,1)), log(unitData.mv2(:,2)) , '*-m');
hold on, plot([-2 10], [-2 10], 'g--');
xlabel('log(mean)')
if whichSubplot == 1
  ylabel('log(var)')
end

% stPR
subplot(5, min(UPF,nUnits), whichSubplot + 4*min(UPF,nUnits));
if isfield(unitData, 'stPRglobal')
  hold on, plot(-100:100, unitData.stPRglobal, 'r');
end
if isfield(unitData, 'stPRsh')
  hold on, plot(-100:100, unitData.stPRsh, 'k');
end
hold on, plot(-100:100, unitData.stPR(1,:), 'b');
hold on, plot(-100:100, unitData.stPR(2,:), 'b');
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

% Firing rate
figTmp2 = figure('units', 'normalized', 'position', [0.002, .04, 1/2, .88/2], 'Visible', visibility);
plot(1:numel(unitData.lfr1),unitData.lfr1);
%       hold on
%       plot((1:numel(unitData.lfr5))*5,unitData.lfr5);
%       hold off
if isPop
  titleStr = sprintf('sh%i_mfr1_%3.1f_mfr2_%3.1f', sh, unitData.mfr_1sthalf, unitData.mfr_2ndhalf);
else
  titleStr = sprintf('sh%i_u%i_mfr1_%3.1f_mfr2_%3.1f', sh, unit, unitData.mfr_1sthalf, unitData.mfr_2ndhalf);
end
title(titleStr, 'Interpreter', 'none')
xlabel('Time (min)')
ylabel('Firing rate (APs/min)')
if strcmpi(figsubdirname,'')
  if isPop
    figName = [entryName '_' figType '_lfr_sh' num2str(sh)];
  else
    figName = [entryName '_' figType '_lfr_sh' num2str(sh) '_u' num2str(unit)];
  end
else
  if isPop
    figName = [figsubdirname filesep entryName '_' figType '_lfr_sh' num2str(sh)];
  else
    figName = [figsubdirname filesep entryName '_' figType '_lfr_sh' num2str(sh) '_u' num2str(unit)];
  end
end
set(figTmp2, 'Name',figName);
hgsave(figTmp2, figName);
close(figTmp2);

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
  titleStr = sprintf('Phase sh%i mfr1_%3.1f mfr2_%3.1f', sh, unitData.mfr_1sthalf, unitData.mfr_2ndhalf);
else
  titleStr = sprintf('Phase sh%i u%i mfr1_%3.1f mfr2_%3.1f', sh, unit, unitData.mfr_1sthalf, unitData.mfr_2ndhalf);
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
  titleStr = sprintf('Coherence sh%i mfr1_%3.1f mfr2_%3.1f', sh, unitData.mfr_1sthalf, unitData.mfr_2ndhalf);
else
  titleStr = sprintf('Coherence sh%i u%i mfr1_%3.1f mfr2_%3.1f', sh, unit, unitData.mfr_1sthalf, unitData.mfr_2ndhalf);
end
title(titleStr, 'Interpreter', 'none')
xlabel('Frequency (Hz)');
ylabel('Coherence');
legend('1st half', '2nd half');
if strcmpi(figsubdirname,'')
  if isPop
    figName = [entryName '_' figType '_half_phase_and_coherence_sh' num2str(sh)];
  else
    figName = [entryName '_' figType '_half_phase_and_coherence_sh' num2str(sh) '_u' num2str(unit)];
  end
else
  if isPop
    figName = [figsubdirname filesep entryName '_' figType '_half_phase_and_coherence_sh' num2str(sh)];
  else
    figName = [figsubdirname filesep entryName '_' figType '_half_phase_and_coherence_sh' num2str(sh) '_u' num2str(unit)];
  end
end
set(figTmp3, 'Name',figName);
hgsave(figTmp3, figName);
close(figTmp3);