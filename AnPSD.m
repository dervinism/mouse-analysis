% Run this script to perform coherence and phase analyses on mouse data and
% extract waveforms. Figure drawing scripts are also no longer used. They
% were replaced by analyses scripts written subsequently.

cleanUp
%parpool(32); % This line should be commented out if not running on NSG Portal


dataFile = 'M191106_MD.mat';
dbEntries = inf;
dbEntries_c = inf;
dbEntries_ca = inf;


params
lists
addDependencies
%downsamplingConv
makedb_new
AnPSD_load
AnPSD_load_waveforms
% eyeLoad
% eyeLoad_filter
% motionLoad
% runningSpeedLoad

% AnPSD_units
% compareHalves
% AnPSD_muas
% AnPSD_unitsLaminar
% eyeAnalysis
% eyeAnalysis_muas
% motionAnalysis

AnPSD_load_subpop_new
AnPSD_unitsLaminar_positive
AnPSD_unitsLaminar_negative

% AnPSD_units_ca
% compareHalves_ca
AnPSD_units_ca_positive
compareHalves_ca_positive
AnPSD_units_ca_negative
compareHalves_ca_negative

% AnPSD_units_ca_muas
% compareHalves_ca_muas
AnPSD_units_ca_muas_positive
compareHalves_ca_muas_positive
AnPSD_units_ca_muas_negative
compareHalves_ca_muas_negative


%% No longer used (but can be run if there is a need)

% AnPSD_units_figs
% compareHalves_figs
% eyeAnalysis_figs
% motionAnalysis_figs
% AnPSD_units_figs_ca
% compareHalves_figs_ca
% AnPSD_units_figs_ca_positive
% compareHalves_figs_ca_positive
% AnPSD_units_figs_ca_negative
% compareHalves_figs_ca_negative

% copyNSG
% AnPSD_segs
% AnPSD_segs_figs
% AnHT
% AnHT_plot
% initPCAm
% compareSeries
% compareSeries_figs
% compareSeriesMultiple
% AnPSD_load_subpopFilter

% AnPSD_units_ca_neutral
% AnPSD_units_figs_ca_neutral
% compareHalves_ca_neutral
% compareHalves_figs_ca_neutral

% lfpLoad
% lfpAnalysis
% lfpAnalysisMotion
% lfpAnalysisPR
% lfpAnalysisPR_ca

% PCAanalysis
% PCAanalysis2
% PCAanalysis2_noS1
% PCAanalysis2_noS1_noVB
% PCAanalysis2_noS1_noRSC
% PCAanalysis2_noS1_noRSC_noVB
% PCAanalysis3_pupil
% PCAanalysis3_motion
% PCAanalysis4