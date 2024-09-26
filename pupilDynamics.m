function [freq, psd] = pupilDynamics(areaFilt, sr, opt)

[freq, psd] = freqDependentWindowCoherenceMD(areaFilt-mean(areaFilt), [], 1/sr, [], opt);

loglog(freq,psd)
ylabel('Power ((Area units)^2/Hz)')
xlabel('Frequency (Hz)')
title('PSD_pupil', 'Interpreter','none');
hgsave(gcf, [filepath filesep 'PSD_pupil']);