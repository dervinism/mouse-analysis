% Run this script to perform pca analysis on mouse data.


%% Load pre-processed data
load(dataFile);


%% Resample data
% LOOP THROUGH DB ENTRIES
fnsData = fieldnames(dataStruct);
for dbCount = numel(fnsData):-1:1
  dbStruct = dataStruct.(fnsData{dbCount});
  shankIDs = fieldnames(dbStruct);

% LOOP THROUGH SHANKS
  for sh = 1:numel(shankIDs)
    shankStruct = dbStruct.(shankIDs{sh});
    fprintf('%s shank %d -------------------\n', fnsData{dbCount}, sh);
    
% LOAD THE CONTENTS OF THE SHANK STRUCTURE VARIABLE
    spk = shankStruct.spk;
    dbID = shankStruct.dbID;
    units = shankStruct.units;
    opt = shankStruct.opt;
    params = shankStruct.params;
    spkAll = shankStruct.spkAll;
    spkAllShanks = shankStruct.spkAllShanks;
    if ~iscell(spkAllShanks)
      spkAllShanks = {spkAllShanks};
    end
    currentGrp = shankStruct.currentGrp;
    load_opt = shankStruct.load_opt;
    data = shankStruct.Q;
    for dbCountShank = 1:numel(shankStruct.db)
      if shankStruct.db(dbCountShank).entryName == fnsData{dbCount}
        animal = shankStruct.db(dbCountShank).animal;
        series = shankStruct.db(dbCountShank).series;
        entryName = shankStruct.db(dbCountShank).entryName;
        shank = shankStruct.db(dbCountShank).shank;
        break
      end
    end
    
% INITIALISE RESAMPLED DATA STORAGE VARIABLES
    resampledFs = [100 10 1 0.1 0.05 0.035 0.03];
    resampledData100 = zeros(size(spk,1), ceil(size(spk,2)*(resampledFs(1)/params.Fs)));
    resampledData10 = zeros(size(spk,1), ceil(size(spk,2)*(resampledFs(2)/params.Fs)));
    resampledData1 = zeros(size(spk,1), ceil(size(spk,2)*(resampledFs(3)/params.Fs)));
    resampledData01 = zeros(size(spk,1), ceil(size(spk,2)*(resampledFs(4)/params.Fs)));
    resampledData005 = zeros(size(spk,1), ceil(size(spk,2)*(resampledFs(5)/params.Fs)));
    resampledData0035 = zeros(size(spk,1), ceil(size(spk,2)*(resampledFs(6)/params.Fs)));
    resampledData003 = zeros(size(spk,1), ceil(size(spk,2)*(resampledFs(7)/params.Fs)));

% LOOP THROUGH UNITS
    for u = 1:numel(load_opt.selectedUnits)
      if ~any(units == load_opt.selectedUnits(u))
        continue
      end
      fprintf('Started processing unit %i\n',load_opt.selectedUnits(u));
      
      spkOI = full(spk(u,:));
      resampledData100(u,:) = resampleData(spkOI, params.Fs, resampledFs(1));
      resampledData10(u,:) = resampleData(spkOI, params.Fs, resampledFs(2));
      resampledData1(u,:) = resampleData(spkOI, params.Fs, resampledFs(3));
      resampledData01(u,:) = resampleData(spkOI, params.Fs, resampledFs(4));
      resampledData005(u,:) = resampleData(spkOI, params.Fs, resampledFs(5));
      resampledData0035(u,:) = resampleData(spkOI, params.Fs, resampledFs(6));
      resampledData003(u,:) = resampleData(spkOI, params.Fs, resampledFs(7));
    end


%% PCA
    FOI = resampledFs;
    pcOI = 10;
    filt = 1; %0;
    [PCAdata.pcaCoef100, PCAdata.explained100, PCAdata.PCs100, PCAdata.PCrho100, PCAdata.PCpval100, PCAdata.PCrhoFilt100,...
      PCAdata.PCpvalFilt100, PCAdata.nPCs100, PCAdata.prob100, pcOI] = pcaAnalysis(resampledData100, resampledFs(1), FOI, pcOI, filt);
    [PCAdata.pcaCoef10, PCAdata.explained10, PCAdata.PCs10, PCAdata.PCrho10, PCAdata.PCpval10, PCAdata.PCrhoFilt10,...
      PCAdata.PCpvalFilt10, PCAdata.nPCs10, PCAdata.prob10] = pcaAnalysis(resampledData10, resampledFs(2), FOI, pcOI, filt);
    [PCAdata.pcaCoef1, PCAdata.explained1, PCAdata.PCs1, PCAdata.PCrho1, PCAdata.PCpval1, PCAdata.PCrhoFilt1,...
      PCAdata.PCpvalFilt1, PCAdata.nPCs1, PCAdata.prob1] = pcaAnalysis(resampledData1, resampledFs(3), FOI, pcOI, filt);
    [PCAdata.pcaCoef01, PCAdata.explained01, PCAdata.PCs01, PCAdata.PCrho01, PCAdata.PCpval01, PCAdata.PCrhoFilt01,...
      PCAdata.PCpvalFilt01, PCAdata.nPCs01, PCAdata.prob01] = pcaAnalysis(resampledData01, resampledFs(4), FOI, pcOI, filt);
    [PCAdata.pcaCoef005, PCAdata.explained005, PCAdata.PCs005, PCAdata.PCrho005, PCAdata.PCpval005, PCAdata.PCrhoFilt005,...
      PCAdata.PCpvalFilt005, PCAdata.nPCs005, PCAdata.prob005] = pcaAnalysis(resampledData005, resampledFs(5), FOI, pcOI, filt);
    [PCAdata.pcaCoef0035, PCAdata.explained0035, PCAdata.PCs0035, PCAdata.PCrho0035, PCAdata.PCpval0035, PCAdata.PCrhoFilt0035,...
      PCAdata.PCpvalFilt0035, PCAdata.nPCs0035, PCAdata.prob0035] = pcaAnalysis(resampledData0035, resampledFs(6), FOI, pcOI, filt);
    [PCAdata.pcaCoef003, PCAdata.explained003, PCAdata.PCs003, PCAdata.PCrho003, PCAdata.PCpval003, PCAdata.PCrhoFilt003,...
      PCAdata.PCpvalFilt003, PCAdata.nPCs003, PCAdata.prob003] = pcaAnalysis(resampledData003, resampledFs(7), FOI, pcOI, filt);


%% PCA correlation with phase
    try
      frequencies = data{1}.freq{1};
    catch
      frequencies = data{1}.freq;
    end
    [PCAdata.phase0_03, PCAdata.coherence0_03, ~, PCAdata.phaseSig0_03, PCAdata.phaseAbsSig0_03, PCAdata.coherenceSig0_03,...
      PCAdata.pcaSig0_03, PCAdata.indsSig0_03] = phaseCohPCA(data, frequencies, 0.03, PCAdata.pcaCoef003(:,1:pcOI));
    [PCAdata.phase0_035, PCAdata.coherence0_035, ~, PCAdata.phaseSig0_035, PCAdata.phaseAbsSig0_035, PCAdata.coherenceSig0_035,...
      PCAdata.pcaSig0_035, PCAdata.indsSig0_035] = phaseCohPCA(data, frequencies, 0.035, PCAdata.pcaCoef0035(:,1:pcOI));
    [PCAdata.phase0_05, PCAdata.coherence0_05, ~, PCAdata.phaseSig0_05, PCAdata.phaseAbsSig0_05, PCAdata.coherenceSig0_05,...
      PCAdata.pcaSig0_05, PCAdata.indsSig0_05] = phaseCohPCA(data, frequencies, 0.05, PCAdata.pcaCoef005(:,1:pcOI));
    [PCAdata.phase0_1, PCAdata.coherence0_1, ~, PCAdata.phaseSig0_1, PCAdata.phaseAbsSig0_1, PCAdata.coherenceSig0_1,...
      PCAdata.pcaSig0_1, PCAdata.indsSig0_1] = phaseCohPCA(data, frequencies, 0.1, PCAdata.pcaCoef01(:,1:pcOI));
    [PCAdata.phase1, PCAdata.coherence1, ~, PCAdata.phaseSig1, PCAdata.phaseAbsSig1, PCAdata.coherenceSig1, PCAdata.pcaSig1,...
      PCAdata.indsSig1] = phaseCohPCA(data, frequencies, 1, PCAdata.pcaCoef1(:,1:pcOI));

    [PCAdata.r0_03, PCAdata.p0_03, PCAdata.rAbs0_03, PCAdata.pAbs0_03] = phaseCorr(PCAdata.phaseSig0_03, PCAdata.phaseAbsSig0_03, PCAdata.pcaSig0_03);
    [PCAdata.r0_035, PCAdata.p0_035, PCAdata.rAbs0_035, PCAdata.pAbs0_035] = phaseCorr(PCAdata.phaseSig0_035, PCAdata.phaseAbsSig0_035, PCAdata.pcaSig0_035);
    [PCAdata.r0_05, PCAdata.p0_05, PCAdata.rAbs0_05, PCAdata.pAbs0_05] = phaseCorr(PCAdata.phaseSig0_05, PCAdata.phaseAbsSig0_05, PCAdata.pcaSig0_05);
    [PCAdata.r0_1, PCAdata.p0_1, PCAdata.rAbs0_1, PCAdata.pAbs0_1] = phaseCorr(PCAdata.phaseSig0_1, PCAdata.phaseAbsSig0_1, PCAdata.pcaSig0_1);
    [PCAdata.r1, PCAdata.p1, PCAdata.rAbs1, PCAdata.pAbs1] = phaseCorr(PCAdata.phaseSig1, PCAdata.phaseAbsSig1, PCAdata.pcaSig1);


%% PCA correlation with unit mean firing rate
    fr = sum(spk,2);
    PCAdata.mfr = fr./(size(spk,2)/1000);
    [PCAdata.rMFR003, PCAdata.pMFR003] = rateCorr(PCAdata.mfr, PCAdata.pcaCoef003(:,1:pcOI), 'Pearson');
    [PCAdata.rMFR0035, PCAdata.pMFR0035] = rateCorr(PCAdata.mfr, PCAdata.pcaCoef0035(:,1:pcOI), 'Pearson');
    [PCAdata.rMFR005, PCAdata.pMFR005] = rateCorr(PCAdata.mfr, PCAdata.pcaCoef005(:,1:pcOI), 'Pearson');
    [PCAdata.rMFR01, PCAdata.pMFR01] = rateCorr(PCAdata.mfr, PCAdata.pcaCoef01(:,1:pcOI), 'Pearson');
    [PCAdata.rMFR1, PCAdata.pMFR1] = rateCorr(PCAdata.mfr, PCAdata.pcaCoef1(:,1:pcOI), 'Pearson');

    [PCAdata.rhoMFR003, PCAdata.pMFR003_Spearman] = rateCorr(PCAdata.mfr, PCAdata.pcaCoef003(:,1:pcOI), 'Spearman');
    [PCAdata.rhoMFR0035, PCAdata.pMFR0035_Spearman] = rateCorr(PCAdata.mfr, PCAdata.pcaCoef0035(:,1:pcOI), 'Spearman');
    [PCAdata.rhoMFR005, PCAdata.pMFR005_Spearman] = rateCorr(PCAdata.mfr, PCAdata.pcaCoef005(:,1:pcOI), 'Spearman');
    [PCAdata.rhoMFR01, PCAdata.pMFR01_Spearman] = rateCorr(PCAdata.mfr, PCAdata.pcaCoef01(:,1:pcOI), 'Spearman');
    [PCAdata.rhoMFR1, PCAdata.pMFR1_Spearman] = rateCorr(PCAdata.mfr, PCAdata.pcaCoef1(:,1:pcOI), 'Spearman');


%% Phase correlation with mean firing rate
    [PCAdata.rMFRphase003, PCAdata.pMFRphase003] = rateCorr(PCAdata.mfr(PCAdata.indsSig0_03), PCAdata.phaseAbsSig0_03, 'Pearson');
    [PCAdata.rMFRphase0035, PCAdata.pMFRphase0035] = rateCorr(PCAdata.mfr(PCAdata.indsSig0_035), PCAdata.phaseAbsSig0_035, 'Pearson');
    [PCAdata.rMFRphase005, PCAdata.pMFRphase005] = rateCorr(PCAdata.mfr(PCAdata.indsSig0_05), PCAdata.phaseAbsSig0_05, 'Pearson');
    [PCAdata.rMFRphase01, PCAdata.pMFRphase01] = rateCorr(PCAdata.mfr(PCAdata.indsSig0_1), PCAdata.phaseAbsSig0_1, 'Pearson');
    [PCAdata.rMFRphase1, PCAdata.pMFRphase1] = rateCorr(PCAdata.mfr(PCAdata.indsSig1), PCAdata.phaseAbsSig1, 'Pearson');

    [PCAdata.rhoMFRphase003, PCAdata.pMFRphase003_Spearman] = rateCorr(PCAdata.mfr(PCAdata.indsSig0_03), PCAdata.phaseAbsSig0_03, 'Spearman');
    [PCAdata.rhoMFRphase0035, PCAdata.pMFRphase0035_Spearman] = rateCorr(PCAdata.mfr(PCAdata.indsSig0_035), PCAdata.phaseAbsSig0_035, 'Spearman');
    [PCAdata.rhoMFRphase005, PCAdata.pMFRphase005_Spearman] = rateCorr(PCAdata.mfr(PCAdata.indsSig0_05), PCAdata.phaseAbsSig0_05, 'Spearman');
    [PCAdata.rhoMFRphase01, PCAdata.pMFRphase01_Spearman] = rateCorr(PCAdata.mfr(PCAdata.indsSig0_1), PCAdata.phaseAbsSig0_1, 'Spearman');
    [PCAdata.rhoMFRphase1, PCAdata.pMFRphase1_Spearman] = rateCorr(PCAdata.mfr(PCAdata.indsSig1), PCAdata.phaseAbsSig1, 'Spearman');


%% Phase correlation with coherence
    [PCAdata.rCoh003, PCAdata.pCoh003] = cohCorr(PCAdata.phaseAbsSig0_03, PCAdata.coherenceSig0_03, 'Pearson');
    [PCAdata.rCoh0035, PCAdata.pCoh0035] = cohCorr(PCAdata.phaseAbsSig0_035, PCAdata.coherenceSig0_035, 'Pearson');
    [PCAdata.rCoh005, PCAdata.pCoh005] = cohCorr(PCAdata.phaseAbsSig0_05, PCAdata.coherenceSig0_05, 'Pearson');
    [PCAdata.rCoh01, PCAdata.pCoh01] = cohCorr(PCAdata.phaseAbsSig0_1, PCAdata.coherenceSig0_1, 'Pearson');
    [PCAdata.rCoh1, PCAdata.pCoh1] = cohCorr(PCAdata.phaseAbsSig1, PCAdata.coherenceSig1, 'Pearson');
    
    [PCAdata.rhoCoh003, PCAdata.pCoh003_Spearman] = cohCorr(PCAdata.phaseAbsSig0_03, PCAdata.coherenceSig0_03, 'Spearman');
    [PCAdata.rhoCoh0035, PCAdata.pCoh0035_Spearman] = cohCorr(PCAdata.phaseAbsSig0_035, PCAdata.coherenceSig0_035, 'Spearman');
    [PCAdata.rhoCoh005, PCAdata.pCoh005_Spearman] = cohCorr(PCAdata.phaseAbsSig0_05, PCAdata.coherenceSig0_05, 'Spearman');
    [PCAdata.rhoCoh01, PCAdata.pCoh01_Spearman] = cohCorr(PCAdata.phaseAbsSig0_1, PCAdata.coherenceSig0_1, 'Spearman');
    [PCAdata.rhoCoh1, PCAdata.pCoh1_Spearman] = cohCorr(PCAdata.phaseAbsSig1, PCAdata.coherenceSig1, 'Spearman');


%% Save the analysis results
    dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.PCAdata = PCAdata;'];
    eval(dataString);
    save(dataFile,'dataStruct','-v7.3');


%% Plot data
    %load(dataFile);
    freq = [0.03 0.035 0.05 0.1 1];
    prefix = [entryName '_' shankIDs{sh} '_'];

    % Coherence vs absolute phase
    plotCohPhase(prefix, PCAdata.coherenceSig0_03, PCAdata.phaseAbsSig0_03, freq(1), PCAdata.rCoh003, PCAdata.pCoh003, 'Pearson');
    plotCohPhase(prefix, PCAdata.coherenceSig0_035, PCAdata.phaseAbsSig0_035, freq(2), PCAdata.rCoh0035, PCAdata.pCoh0035, 'Pearson');
    plotCohPhase(prefix, PCAdata.coherenceSig0_05, PCAdata.phaseAbsSig0_05, freq(3), PCAdata.rCoh005, PCAdata.pCoh005, 'Pearson');
    plotCohPhase(prefix, PCAdata.coherenceSig0_1, PCAdata.phaseAbsSig0_1, freq(4), PCAdata.rCoh01, PCAdata.pCoh01, 'Pearson');
    plotCohPhase(prefix, PCAdata.coherenceSig1, PCAdata.phaseAbsSig1, freq(5), PCAdata.rCoh1, PCAdata.pCoh1, 'Pearson');
    
    plotCohPhase(prefix, PCAdata.coherenceSig0_03, PCAdata.phaseAbsSig0_03, freq(1), PCAdata.rhoCoh003, PCAdata.pCoh003_Spearman, 'Spearman');
    plotCohPhase(prefix, PCAdata.coherenceSig0_035, PCAdata.phaseAbsSig0_035, freq(2), PCAdata.rhoCoh0035, PCAdata.pCoh0035_Spearman, 'Spearman');
    plotCohPhase(prefix, PCAdata.coherenceSig0_05, PCAdata.phaseAbsSig0_05, freq(3), PCAdata.rhoCoh005, PCAdata.pCoh005_Spearman, 'Spearman');
    plotCohPhase(prefix, PCAdata.coherenceSig0_1, PCAdata.phaseAbsSig0_1, freq(4), PCAdata.rhoCoh01, PCAdata.pCoh01_Spearman, 'Spearman');
    plotCohPhase(prefix, PCAdata.coherenceSig1, PCAdata.phaseAbsSig1, freq(5), PCAdata.rhoCoh1, PCAdata.pCoh1_Spearman, 'Spearman');

    % mfr vs absolute phase
    plotMFRPhase(prefix, PCAdata.mfr(PCAdata.indsSig0_03), PCAdata.phaseAbsSig0_03, freq(1), PCAdata.rMFRphase003, PCAdata.pMFRphase003, 'Pearson');
    plotMFRPhase(prefix, PCAdata.mfr(PCAdata.indsSig0_035), PCAdata.phaseAbsSig0_035, freq(2), PCAdata.rMFRphase0035, PCAdata.pMFRphase0035, 'Pearson');
    plotMFRPhase(prefix, PCAdata.mfr(PCAdata.indsSig0_05), PCAdata.phaseAbsSig0_05, freq(3), PCAdata.rMFRphase005, PCAdata.pMFRphase005, 'Pearson');
    plotMFRPhase(prefix, PCAdata.mfr(PCAdata.indsSig0_1), PCAdata.phaseAbsSig0_1, freq(4), PCAdata.rMFRphase01, PCAdata.pMFRphase01, 'Pearson');
    plotMFRPhase(prefix, PCAdata.mfr(PCAdata.indsSig1), PCAdata.phaseAbsSig1, freq(5), PCAdata.rMFRphase1, PCAdata.pMFRphase1, 'Pearson');

    plotMFRPhase(prefix, PCAdata.mfr(PCAdata.indsSig0_03), PCAdata.phaseAbsSig0_03, freq(1), PCAdata.rhoMFRphase003, PCAdata.pMFRphase003_Spearman, 'Spearman');
    plotMFRPhase(prefix, PCAdata.mfr(PCAdata.indsSig0_035), PCAdata.phaseAbsSig0_035, freq(2), PCAdata.rhoMFRphase0035, PCAdata.pMFRphase0035_Spearman, 'Spearman');
    plotMFRPhase(prefix, PCAdata.mfr(PCAdata.indsSig0_05), PCAdata.phaseAbsSig0_05, freq(3), PCAdata.rhoMFRphase005, PCAdata.pMFRphase005_Spearman, 'Spearman');
    plotMFRPhase(prefix, PCAdata.mfr(PCAdata.indsSig0_1), PCAdata.phaseAbsSig0_1, freq(4), PCAdata.rhoMFRphase01, PCAdata.pMFRphase01_Spearman, 'Spearman');
    plotMFRPhase(prefix, PCAdata.mfr(PCAdata.indsSig1), PCAdata.phaseAbsSig1, freq(5), PCAdata.rhoMFRphase1, PCAdata.pMFRphase1_Spearman, 'Spearman');

    % PCA coefficients vs absolute phase
    plotPCAPhase(prefix, PCAdata.pcaCoef003(PCAdata.indsSig0_03,1:pcOI), PCAdata.phaseAbsSig0_03, freq(1), PCAdata.rAbs0_03, PCAdata.pAbs0_03);
    plotPCAPhase(prefix, PCAdata.pcaCoef0035(PCAdata.indsSig0_035,1:pcOI), PCAdata.phaseAbsSig0_035, freq(2), PCAdata.rAbs0_035, PCAdata.pAbs0_035);
    plotPCAPhase(prefix, PCAdata.pcaCoef005(PCAdata.indsSig0_05,1:pcOI), PCAdata.phaseAbsSig0_05, freq(3), PCAdata.rAbs0_05, PCAdata.pAbs0_05);
    plotPCAPhase(prefix, PCAdata.pcaCoef01(PCAdata.indsSig0_1,1:pcOI), PCAdata.phaseAbsSig0_1, freq(4), PCAdata.rAbs0_1, PCAdata.pAbs0_1);
    plotPCAPhase(prefix, PCAdata.pcaCoef1(PCAdata.indsSig1,1:pcOI), PCAdata.phaseAbsSig1, freq(5), PCAdata.rAbs1, PCAdata.pAbs1);

    % PCA coefficients vs mfr
    plotPCAmfr(prefix, PCAdata.pcaCoef003(:,1:pcOI), PCAdata.mfr, freq(1), PCAdata.rMFR003, PCAdata.pMFR003, 'Pearson');
    plotPCAmfr(prefix, PCAdata.pcaCoef0035(:,1:pcOI), PCAdata.mfr, freq(2), PCAdata.rMFR0035, PCAdata.pMFR0035, 'Pearson');
    plotPCAmfr(prefix, PCAdata.pcaCoef005(:,1:pcOI), PCAdata.mfr, freq(3), PCAdata.rMFR005, PCAdata.pMFR005, 'Pearson');
    plotPCAmfr(prefix, PCAdata.pcaCoef01(:,1:pcOI), PCAdata.mfr, freq(4), PCAdata.rMFR01, PCAdata.pMFR01, 'Pearson');
    plotPCAmfr(prefix, PCAdata.pcaCoef1(:,1:pcOI), PCAdata.mfr, freq(5), PCAdata.rMFR1, PCAdata.pMFR1, 'Pearson');

    plotPCAmfr(prefix, PCAdata.pcaCoef003(:,1:pcOI), PCAdata.mfr, freq(1), PCAdata.rhoMFR003, PCAdata.pMFR003_Spearman, 'Spearman');
    plotPCAmfr(prefix, PCAdata.pcaCoef0035(:,1:pcOI), PCAdata.mfr, freq(2), PCAdata.rhoMFR0035, PCAdata.pMFR0035_Spearman, 'Spearman');
    plotPCAmfr(prefix, PCAdata.pcaCoef005(:,1:pcOI), PCAdata.mfr, freq(3), PCAdata.rhoMFR005, PCAdata.pMFR005_Spearman, 'Spearman');
    plotPCAmfr(prefix, PCAdata.pcaCoef01(:,1:pcOI), PCAdata.mfr, freq(4), PCAdata.rhoMFR01, PCAdata.pMFR01_Spearman, 'Spearman');
    plotPCAmfr(prefix, PCAdata.pcaCoef1(:,1:pcOI), PCAdata.mfr, freq(5), PCAdata.rhoMFR1, PCAdata.pMFR1_Spearman, 'Spearman');
  end % Loop through shanks
end % Loop through db entries



%% Functions
function [pcaCoef, explained, PCs, PCrho, PCpval, PCrhoFilt, PCpvalFilt, nPCs, prob, pcOI] = pcaAnalysis(data, Fs, FOI, pcOI, filt)

biggestIndex = size(data,2);

[pcaCoef, PCs, ~, ~, explained] = pca(data');
explained = explained';
PCs = PCs';
[nPCs, prob] = barttest(pcaCoef,0.05);
if size(PCs,1) < pcOI
  pcOI = size(PCs,1);
end

pr = sum(data);

if filt
  % BAND-PASS FILTER THE POPULATION RATE AND APPLY HILBERT TRANSFORM
  % Design the filters
  Fstop1 = 03; %17.5; %28.5;
  Fpass1 = 23; %27.5; %29.5;
  Fpass2 = 30;
  Fpass3 = 37; %32.5; %30.5;
  Fstop2 = 57; %42.5; %31.5;
  Fpass = [Fstop1; Fpass1; Fpass2; Fpass3; Fstop2];
  FpassLog = log10(Fpass);
  FpassLog = FpassLog-FpassLog(3);
  
  FOIlog = log10(FOI);
  FpassLog = [FOIlog(1)+FpassLog FOIlog(2)+FpassLog FOIlog(3)+FpassLog ...
    FOIlog(4)+FpassLog FOIlog(5)+FpassLog FOIlog(6)+FpassLog FOIlog(7)+FpassLog];
  Fpass = 10.^FpassLog;
  
  Astop1 = 65;
  Apass  = 0.5;
  Astop2 = 65;
  
  %     d003 = designFilterLO(Fpass(:,1), Astop1, Apass, Astop2, Fs);
  %     d01  = designFilterLO(Fpass(:,2), Astop1, Apass, Astop2, Fs);
  %     d03  = designFilterLO(Fpass(:,3), Astop1, Apass, Astop2, Fs);
  %     d1   = designFilterLO(Fpass(:,4), Astop1, Apass, Astop2, Fs);
  %     d3   = designFilterLO(Fpass(:,5), Astop1, Apass, Astop2, Fs);
  %     d10  = designFilterLO(Fpass(:,6), Astop1, Apass, Astop2, Fs);
  %d30  = designFilterLO(Fpass(:,7), Astop1, Apass, Astop2, Fs);
  d003 = designFilter(Fpass(:,1), Astop1, Apass, Astop2, Fs);
  d01  = designFilter(Fpass(:,2), Astop1, Apass, Astop2, Fs);
  d03  = designFilter(Fpass(:,3), Astop1, Apass, Astop2, Fs);
  d1   = designFilter(Fpass(:,4), Astop1, Apass, Astop2, Fs);
  d3   = designFilter(Fpass(:,5), Astop1, Apass, Astop2, Fs);
  d10  = designFilter(Fpass(:,6), Astop1, Apass, Astop2, Fs);
  d30  = designFilter(Fpass(:,7), Astop1, Apass, Astop2, Fs);
  d = {d003 d01 d03 d1 d3 d10 d30};
  
  %     fvtool(d003, 'FrequencyScale','Log');
  
  % Filter the signal
  prFilt = zeros(7,biggestIndex);
  for f = 1:numel(FOI)
    prFilt(f,:) = filtfilt(d{f},pr);
  end
  
  % Visualise the data
  figH = zeros(1,length(FOI));
  for f = 1:numel(FOI)
    figH(f) = figure('units', 'normalized', 'position', [0.002, .04, 1, .88]);
    plot(pr, 'b', 'LineWidth', 1);
    hold on
    plot(prFilt(f,:), 'g', 'LineWidth', 1);
    hold off
    legend('population rate', 'filtered population rate');
    title(['Band-pass: ' num2str(FOI(f)) ' Hz']);
  end
  
  
  % Correlate PCs with pr
  PCrhoFilt = zeros(numel(FOI),pcOI);
  PCpvalFilt = zeros(numel(FOI),pcOI);
  for f = 1:numel(FOI)
    for p = 1:pcOI
      [rho, pval] = corr([prFilt(f,:)' PCs(p,:)']);
      PCrhoFilt(f,p) = rho(2);
      PCpvalFilt(f,p) = pval(2);
    end
  end
else
  PCrhoFilt = [];
  PCpvalFilt = [];
end

PCrho = zeros(1,pcOI);
PCpval = zeros(1,pcOI);
for p = 1:pcOI
  [rho, pval] = corr([pr' PCs(p,:)']);
  PCrho(p) = rho(2);
  PCpval(p) = pval(2);
end
end

function [phase, coherence, phaseSig, phaseLeadLagSig, phaseAbsSig, coherenceSig, pcaSig, indsSig] = phaseCohPCA(data, frequencies, frequency, pcaCoef)

[~, i] = min(abs(frequencies - frequency));

phase = NaN(1,numel(data));
coherence = NaN(1,numel(data));

phaseSig = [];
phaseLeadLagSig = [];
phaseAbsSig = [];
pcaSig = [];
coherenceSig = [];
indsSig = [];

for u = 1:numel(data)
  try
    phase(u) = data{u}.phase{1}(i);
  catch
    phase(u) = data{u}.phase(i);
  end
  try
    coherence(u) = data{u}.coherencyAdj{1}(i);
  catch
    try
      coherence(u) = data{u}.coherenceAdj{1}(i);
    catch
      coherence(u) = data{u}.coh(i);
    end
  end
  if ~isnan(phase(u))
    phaseSig = [phaseSig; phase(u)]; %#ok<*AGROW>
    phaseLeadLagSig = [phaseLeadLagSig; phaseLeadLag(phase(u))];
    phaseAbsSig = [phaseAbsSig; abs(phaseLeadLag(phase(u)))];
    coherenceSig = [coherenceSig; coherence(u)];
    pcaSig = [pcaSig; pcaCoef(u,:)];
    indsSig = [indsSig; u];
  end
end
end

function [phase] = phaseLeadLag(phase)

phase = rem(phase, 2*pi);
if phase > pi && phase < 2*pi
  phase = phase - 2*pi;
end
end

function [r, p, rAbs, pAbs] = phaseCorr(phase, phaseAbs, pcaCoef)

r = zeros(1,size(pcaCoef,2));
p = zeros(1,size(pcaCoef,2));
rAbs = zeros(1,size(pcaCoef,2));
pAbs = zeros(1,size(pcaCoef,2));
for i = 1:size(pcaCoef,2)
  [R, P] = corr([phase pcaCoef(:,i)]);
  r(i) = R(2);
  p(i) = P(2);
  
  [R, P] = corr([phaseAbs pcaCoef(:,i)]);
  rAbs(i) = R(2);
  pAbs(i) = P(2);
end
end

function [r, p] = rateCorr(mfr, pcaCoef, type)

r = zeros(1,size(pcaCoef,2));
p = zeros(1,size(pcaCoef,2));
for i = 1:size(pcaCoef,2)
  [R, P] = corr([full(mfr) full(pcaCoef(:,i))], 'type',type);
  r(i) = R(2);
  p(i) = P(2);
end
end

function [r, p] = cohCorr(phase, coherence, type)

i = isnan(coherence);
phase = phase(~i);
coherence = coherence(~i);
if isempty(coherence)
  r = 0; p = NaN;
else
  [R, P] = corr([phase coherence], 'type',type);
  r = R(2);
  p = P(2);
end
end

function figH = plotCohPhase(prefix, coherenceSig, phaseAbsSig, freq, r, p, type)

i = isnan(coherenceSig);
coherenceSig = coherenceSig(~i);
phaseAbsSig = phaseAbsSig(~i);
if isempty(coherenceSig) || isempty(phaseAbsSig)
  figH = [];
  return
end

figH = figure('units', 'normalized', 'position', [0.002, .04, 1, .88]);
plot(coherenceSig, phaseAbsSig, '.', 'markerSize', 10)
if strcmp(type, 'Pearson')
  coeffs = polyfit(coherenceSig, phaseAbsSig, 1);
  fittedX = linspace(min(coherenceSig), max(coherenceSig), 200);
  fittedY = polyval(coeffs, fittedX);
  hold on
  plot(fittedX, fittedY, 'r-', 'LineWidth', 1);
  hold off
end
xlabel('Coherence');
ylabel('abs(phase) (rad)');
if strcmp(type, 'Pearson')
  title(['Coherence vs phase (f = ' num2str(freq) ' Hz, ' type ' r = ' num2str(r) ', p = ' num2str(p) ')']);
  figFileName = [prefix 'coherence_vs_phase___f' num2str(freq) 'Hz___' type ' r' num2str(r) '___p' num2str(p)];
else
  title(['Coherence vs phase (f = ' num2str(freq) ' Hz, ' type ' rho = ' num2str(r) ', p = ' num2str(p) ')']);
  figFileName = [prefix 'coherence_vs_phase___f' num2str(freq) 'Hz___' type ' rho' num2str(r) '___p' num2str(p)];
end
figFileName = strrep(figFileName,'.','p');
saveas(figH, figFileName, 'png')
end

function figHr = plotMFRPhase(prefix, mfr, phaseAbsSig, freq, r, p, type)

figHr = figure('units', 'normalized', 'position', [0.002, .04, 1, .88]);
if strcmp(type, 'Pearson')
  plot(mfr, phaseAbsSig, '.', 'markerSize', 10)
  coeffs = polyfit(mfr, phaseAbsSig, 1);
  fittedX = linspace(min(mfr), max(mfr), 200);
  fittedY = polyval(coeffs, fittedX);
  hold on
  plot(fittedX, fittedY, 'r-', 'LineWidth', 1);
  hold off
else
  semilogx(mfr, phaseAbsSig, '.', 'markerSize', 10)
end
xlabel('Mean firing rate (APs/s)');
ylabel('abs(phase) (rad)');
if strcmp(type, 'Pearson')
  title(['Mean firing rate vs phase (f = ' num2str(freq) ' Hz, ' type ' r = ' num2str(r) ', p = ' num2str(p) ')']);
  figFileName = [prefix 'mfr_vs_phase___f' num2str(freq) 'Hz___' type '___r' num2str(r) '___p' num2str(p)];
else
  title(['Mean firing rate vs phase (f = ' num2str(freq) ' Hz, ' type ' rho = ' num2str(r) ', p = ' num2str(p) ')']);
  figFileName = [prefix 'mfr_vs_phase___f' num2str(freq) 'Hz___' type '___rho' num2str(r) '___p' num2str(p)];
end
figFileName = strrep(figFileName,'.','p');
saveas(figHr, figFileName, 'png')
end

function figH = plotPCAPhase(prefix, pcaCoef, phase, freq, r, p)

for i = 1:size(pcaCoef,2)
  figH = figure('units', 'normalized', 'position', [0.002, .04, 1, .88]);
  plot(phase, pcaCoef(:,i), '.', 'markerSize', 10)
  coeffs = polyfit(phase, pcaCoef(:,i), 1);
  fittedX = linspace(min(phase), max(phase), 200);
  fittedY = polyval(coeffs, fittedX);
  hold on
  plot(fittedX, fittedY, 'r-', 'LineWidth', 1);
  hold off
  %ylim([-1 1]);
  xlabel('abs(phase) (rad)');
  ylabel('PCA coefficients');
  title(['Phase vs PCA coefficients (f = ' num2str(freq) ' Hz, PC' num2str(i) ', r = ' num2str(r(i)) ', p = ' num2str(p(i)) ')']);
  figFileName = [prefix 'phase_vs_pca___f' num2str(freq) 'Hz___PC' num2str(i), '___r' num2str(r(i)) '___p' num2str(p(i))];
  figFileName = strrep(figFileName,'.','p');
  saveas(figH, figFileName, 'png')
end
end

function figH = plotPCAmfr(prefix, pcaCoef, mfr, freq, r, p, type)

for i = 1:size(pcaCoef,2)
  figH = figure('units', 'normalized', 'position', [0.002, .04, 1, .88]);
  if strcmp(type, 'Pearson')
    plot(mfr, pcaCoef(:,i), '.', 'markerSize', 10)
    coeffs = polyfit(mfr, pcaCoef(:,i), 1);
    fittedX = linspace(min(mfr), max(mfr), 200);
    fittedY = polyval(coeffs, fittedX);
    hold on
    plot(fittedX, fittedY, 'r-', 'LineWidth', 1);
    hold off
  else
    semilogx(mfr, pcaCoef(:,i), '.', 'markerSize', 10)
  end
  %ylim([-1 1]);
  xlabel('Mean firing rate (APs/s)');
  ylabel('PCA coefficients');
  if strcmp(type, 'Pearson')
    title(['Mean firing rate vs PCA coefficients (f = ' num2str(freq) ' Hz, PC' num2str(i) ', ' type ' r = ' num2str(r(i)) ', p = ' num2str(p(i)) ')']);
    figFileName = [prefix 'mfr_vs_pca___f' num2str(freq) 'Hz___PC' num2str(i), '___' type '___r' num2str(r(i)) '___p' num2str(p(i))];
  else
    title(['Mean firing rate vs PCA coefficients (f = ' num2str(freq) ' Hz, PC' num2str(i) ', ' type ' rho = ' num2str(r(i)) ', p = ' num2str(p(i)) ')']);
    figFileName = [prefix 'mfr_vs_pca___f' num2str(freq) 'Hz___PC' num2str(i), '___' type '___rho' num2str(r(i)) '___p' num2str(p(i))];
  end
  figFileName = strrep(figFileName,'.','p');
  saveas(figH, figFileName, 'png')
end
end