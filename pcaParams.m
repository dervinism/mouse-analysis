%% INITIALISE PARAMETERS
dssrLFPfinal = 10; % final down-sampled LFP sampling rate.
chunkSize = 4500000;
bandNames = {'n slow'; 'delta'; 'alpha'; 'theta'; 'beta'; 's gamma'; 'f gamma'; 'ripples/uf'};
LFPbands  = { [0.1 1];   [1 4];   [4 8];  [8 12]; [12 30]; [30 50];   [50 120];   [120 200]};
medianSubtracted = false;
deleteChans = [];
intermediateSaving = false;