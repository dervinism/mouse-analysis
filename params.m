% Commonly used data analysis parameters

% input/output
uolSourceDirNeuronexus = 'S:\cortical_dynamics\Shared\CSN\Shared'; %'R:\CSN\Shared\Dynamics\Data';
uolSourceDirNeuropixels = 'S:\cortical_dynamics\Shared\CSN\Shared'; %'R:\Neuropix\Shared\Data';
allensdkSourceDir = 'E:\res\allensdk'; %'R:\Neuropix\Shared\Data\allensdk';
dataDir = 'D:\mouse_analysis_data'; 'R:\CSN\Shared\Dynamics\Data\mouse_analysis_data';
paDir = 'area_properties';
paDir_positive = 'area_properties_positive';
paDir_negative = 'area_properties_negative';
paDir_uol = 'area_properties_uol';
paDir_uol_positive = 'area_properties_uol_positive';
paDir_uol_negative = 'area_properties_uol_negative';
paDir_allensdk = 'area_properties_allensdk';
paDir_allensdk_positive = 'area_properties_allensdk_positive';
paDir_allensdk_negative = 'area_properties_allensdk_negative';
laDir = 'area2itself_comparisons';
laDir_uol = 'area2itself_comparisons_uol';
laDir_uol_positive = 'area2itself_comparisons_uol_positive';
laDir_uol_negative = 'area2itself_comparisons_uol_negative';
laDir_allensdk = 'area2itself_comparisons_allensdk';
laDir_allensdk_positive = 'area2itself_comparisons_allensdk_positive';
laDir_allensdk_negative = 'area2itself_comparisons_allensdk_negative';
caDir = 'area2area_comparisons';
caDir_positive = 'area2area_comparisons_positive';
caDir_negative = 'area2area_comparisons_negative';
caDir_uol = 'area2area_comparisons_uol';
caDir_uol_positive = 'area2area_comparisons_uol_positive';
caDir_uol_negative = 'area2area_comparisons_uol_negative';
caDir_allensdk = 'area2area_comparisons_allensdk';
caDir_allensdk_positive = 'area2area_comparisons_allensdk_positive';
caDir_allensdk_negative = 'area2area_comparisons_allensdk_negative';
clDir = 'layer2layer_comparisons';
clDir_positive = 'layer2layer_comparisons_positive';
clDir_negative = 'layer2layer_comparisons_negative';
clDir_uol = 'layer2layer_comparisons_uol';
clDir_uol_positive = 'layer2layer_comparisons_uol_positive';
clDir_uol_negative = 'layer2layer_comparisons_uol_negative';
clDir_allensdk = 'layer2layer_comparisons_allensdk';
clDir_allensdk_positive = 'layer2layer_comparisons_allensdk_positive';
clDir_allensdk_negative = 'layer2layer_comparisons_allensdk_negative';
pupilDir = 'pupilArea';
PRsFolder = 'PRs';
MUAsFolder = 'MUAs';
unitsFolder = 'units';
qualityUnitsFolder = 'unitsQuality';
significantUnitsFolder = 'unitsSignificant';
phaseFrequencyProfilesSubfolder = 'phaseFrequencyProfiles';
meanPhaseSumsSubfolder = 'meanPhaseSums';
phaseDiffProfilesSubfolder = 'phaseDiffProfiles';
phaseSumPredictionsSubfolder = 'phaseSumPredictions';
coherenceFrequencyProfilesSubfolder = 'coherenceFrequencyProfiles';
coherenceCorrelationsSubfolder = 'coherenceCorrelations';
mapsSubfolder = 'maps';
histosSubfolder = 'histograms';
phaseCorrelationsSubfolder = 'phaseCorrelations';
phaseLocationProfilesSubfolder = 'phaseLocationProfiles';
betaCorrelationsSubfolder = 'betaCorrelations';
acgdSubfolder = 'autocorrelationDecays';
waveformsSubfolder = 'waveforms';
firingRatesSubfolder = 'firingRates';
firingRatesSubfolder2 = 'firingRates2';
conditionsSubfolder = 'wakefulness_anaesthesia';
area2pupilDir = 'area2pupil_comparisons';
area2pupilDir_positive = 'area2pupil_comparisons_positive';
area2pupilDir_negative = 'area2pupil_comparisons_negative';
area2pupilDir_uol = 'area2pupil_comparisons_uol';
area2pupilDir_uol_positive = 'area2pupil_comparisons_uol_positive';
area2pupilDir_uol_negative = 'area2pupil_comparisons_uol_negative';
area2pupilDir_allensdk = 'area2pupil_comparisons_allensdk';
area2pupilDir_allensdk_positive = 'area2pupil_comparisons_allensdk_positive';
area2pupilDir_allensdk_negative = 'area2pupil_comparisons_allensdk_negative';
area2motionDir = 'area2motion_comparisons';
area2motionDir_uol = 'area2motion_comparisons_uol';
area2motionDir_allensdk = 'area2motion_comparisons_allensdk';
lfp2pupilDir = 'lfp2pupil_comparisons';
lfp2motionDir = 'lfp2motion_comparisons';
lfp2lfpDir = 'lfp2lfp_comparisons';
lfp2lfpDir_uol = 'lfp2lfp_comparisons_uol';
lfp2lfpDir_allensdk = 'lfp2lfp_comparisons_allensdk';

% data sampling
srData = 60; % data sampling rate. Not to be confused with recording sampling rate
srRecording = 3e4; % recording sampling rate

% coherence analyses parameters
fRef = 0.03; % Hz
maxFreq = 20; %120;
maxFreq_ca = 20;
maxFreq_pupil = 2.79396772384644;
maxFreq_motion = 2.79396772384644;
winfactor = 4; %10;
freqfactor = 1.6; %1.333;
tapers = 3; %5;
FOI = [20 15 10 8 6 5 4 3 2 1 0.7 0.5 0.3 0.2 0.1...
  0.07 0.05 0.03 0.02 0.01]; % frequencies of interest (Hz)

% exclusion radius around a unit when calculating local population firing rate
exclRad = 60; %um

% unit quality check parameters
refractCont = 0.2;
cluDist = 20;
minSpikeCount = 300;

% phase histograms
phaseCentre = 2*pi/4;
phaseLim = [-pi pi] + phaseCentre;
edges = (-pi:pi/8:pi) + phaseCentre;
phaseHistoBinCentres = (pi/8:pi/8:2*pi) - pi/16;

% frequency profiles
freqLimUOL = [10e-3 - 0.002   20];
freqLimAllen = [10e-3 - 0.002   20];

% unit spiking autocorrelations
acgdPeriod = 5; %s

% figure properties
figSize = 15; % cm

% significance level
alpha = 0.05;

% exclude periods when animal is running (relevant to Allen data only)
excludeRunning = true;