% Run this script to perform Hilbert transform analyses on mouse data.


% LOAD PRE-PROCESSED DATA
if ~exist('dataStruct','var')
  load(dataFile);
end


% INITIALISE PARALLEL POOL
% delete(gcp('nocreate'))
% pool = parpool;
% feature('numcores');


% LOOP THROUGH DB ENTRIES
fnsData = fieldnames(dataStruct);
for dbCount = numel(fnsData):-1:1
%for dbCount = 3
  dbStruct = dataStruct.(fnsData{dbCount});
  shankIDs = fieldnames(dbStruct);

% LOOP THROUGH SHANKS
  for sh = 1:numel(shankIDs)
    shankStruct = dbStruct.(shankIDs{sh});
    fprintf('Processing %s shank %d \n', fnsData{dbCount}, sh);
    
% LOAD THE CONTENTS OF THE SHANK STRUCTURE VARIABLE
    spkAll = shankStruct.spkAll;
    spk = shankStruct.spk;
    params = shankStruct.params;
    units = shankStruct.units;
    for dbCountShank = 1:numel(shankStruct.db)
      if shankStruct.db(dbCountShank).entryName == fnsData{dbCount}
        entryName = shankStruct.db(dbCountShank).entryName;
        break
      end
    end

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
    
    FOI = [0.03 0.1 0.3 1 3 10 30];
    FOIlog = log10(FOI);
    FpassLog = [FOIlog(1)+FpassLog FOIlog(2)+FpassLog FOIlog(3)+FpassLog ...
      FOIlog(4)+FpassLog FOIlog(5)+FpassLog FOIlog(6)+FpassLog FOIlog(7)+FpassLog];
    Fpass = 10.^FpassLog;
    
    Astop1 = 65;
    Apass  = 0.5;
    Astop2 = 65;
    Fs = params.Fs;

%     d003 = designFilterLO(Fpass(:,1), Astop1, Apass, Astop2, Fs);
%     d01  = designFilterLO(Fpass(:,2), Astop1, Apass, Astop2, Fs);
%     d03  = designFilterLO(Fpass(:,3), Astop1, Apass, Astop2, Fs);
%     d1   = designFilterLO(Fpass(:,4), Astop1, Apass, Astop2, Fs);
%     d3   = designFilterLO(Fpass(:,5), Astop1, Apass, Astop2, Fs);
%     d10  = designFilterLO(Fpass(:,6), Astop1, Apass, Astop2, Fs);
%     d30  = designFilterLO(Fpass(:,7), Astop1, Apass, Astop2, Fs);
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
    spkAllFilt = zeros(7,length(spkAll));
    for f = 1:numel(FOI)
      spkAllFilt(f,:) = filtfilt(d{f},spkAll);
    end
    %spkAllFilt(1,:) = spkAllFilt(1,:) - mean(spkAllFilt(1,:));
    %spkAllFilt = removeFrequencies(spkAll, params.Fs, 1.1, 1e3);
    %spkAllFilt = removeFrequencies(spkAllFilt, params.Fs, 0.001, 0.9);

    % Hilbert tranform
    [hPop, phasePop, envelopePop] = hilbertTransform(spkAllFilt);
    
    % Find the zero phase
    t = (1:length(spkAll))*(1/params.Fs);
    [t0fullPop, i0fullPop] = zeroPhase(FOI, t, phasePop);
    
    % Visualise the data
    figH = zeros(1,length(FOI));
    for f = 1:numel(FOI)
      figH(f) = figure('units', 'normalized', 'position', [0.002, .04, 1, .88]);
      plot(t,spkAll, 'b', 'LineWidth', 1);
      hold on
      plot(t,spkAllFilt(f,:), 'g', 'LineWidth', 1);
      plot(t,envelopePop(f,:), 'r--', 'LineWidth', 3);
      plot(t,phasePop(f,:), 'k');
      plot(t0fullPop{f},zeros(size(t0fullPop{f})), 'r.', 'MarkerSize',20);
      legend('population rate', 'filtered population rate', 'envelope', 'phase', 'zero phase');
      title(['Band-pass: ' num2str(FOI(f)) ' Hz']);
    end

% LOOP THROUGH UNITS
    nUnits = numel(units);
    meanPhaseRelation = zeros(length(FOI),nUnits);
    meanPhaseRelationIS = zeros(length(FOI)-1,nUnits);
    meanPhaseRelationSD = zeros(1,nUnits);
    meanPhaseRelationSS = zeros(1,nUnits);
    meanPhaseRelationSF = zeros(1,nUnits);
    nPairs = nchoosek(numel(FOI),2);
    pairs = nchoosek(1:numel(FOI),2);
    for u = 1:nUnits
      fprintf('Processing %s shank %d unit %d\n', fnsData{dbCount}, sh, u);
      spkU = full(spk(u,:));
      
% BAND-PASS FILTER THE UNIT RATE AND APPLY HILBERT TRANSFORM
      % Filter the signal
      spkFilt = zeros(7,length(spkU));
      for f = 1:numel(FOI)
        spkFilt(f,:) = filtfilt(d{f},spkU);
      end
      
      % Hilbert tranform
      [h, phase, envelope] = hilbertTransform(spkFilt);

      % Find the zero phase
      [t0full{u}, i0full{u}] = zeroPhase(FOI, t, phase);
      
      % Visualise the data
%       figH = zeros(1,length(FOI));
%       for f = 1:numel(FOI)
%         figH(f) = figure('units', 'normalized', 'position', [0.002, .04, 1, .88]);
%         plot(t,spkU, 'b', 'LineWidth', 1);
%         hold on
%         plot(t,spkFilt(f,:), 'g', 'LineWidth', 1);
%         plot(t,envelope(f,:), 'r--', 'LineWidth', 3);
%         plot(t,phase(f,:), 'k');
%         plot(t0full{u}{f},zeros(size(t0full{u}{f})), 'r.', 'MarkerSize',20);
%         legend('unit rate', 'filtered unit rate', 'envelope', 'phase', 'zero phase');
%         title(['Band-pass: ' num2str(FOI(f)) ' Hz']);
%       end

% COMPARE WITH THE POPULATION RATE
      edges = -pi:pi/8:pi;
      edgesC = -pi+pi/16:pi/8:pi-pi/16;
      
      % Same frequency
      [~,meanPhaseRelation(:,u),MPRCu{u},MPRCl{u},histPhaseRelation{u},~] = phaseRelations(i0full{u},t0full{u},t,phasePop,edges,edgesC);
      [individualPhaseRelationSpk{u},meanPhaseRelationSpk(:,u),MPRSpkCu{u},MPRSpkCl{u},histPhaseRelationSpk{u}] = phaseRelationsSpk(spkU,phasePop,edges,edgesC);
      PLV{u} = phaseLocking(spkU,phasePop);
      [MI{u}, MIRange] = modulationIndex(histPhaseRelationSpk{u});
      r{u} = zeros(nPairs,3);
      for p = 1:nPairs
%         if pairs(p,1) < pairs(p,2)
%           plot(individualPhaseRelationSpk{u}(pairs(p,1),:),individualPhaseRelationSpk{u}(pairs(p,2),:), '.k', 'MarkerSize',10);
%           xlabel(['Phase @ ' num2str(FOI(pairs(p,1))) ' Hz (rad)']);
%           ylabel(['Phase @ ' num2str(FOI(pairs(p,2))) ' Hz (rad)']);
%           title(['Phase @ ' num2str(FOI(pairs(p,1))) ' Hz   vs   phase @ ' num2str(FOI(pairs(p,2))) ' Hz']);
%         else
%           plot(individualPhaseRelationSpk{u}(pairs(p,2),:),individualPhaseRelationSpk{u}(pairs(p,1),:), '.k', 'MarkerSize',10);
%           xlabel(['Phase @ ' num2str(FOI(pairs(p,2))) ' Hz (rad)']);
%           ylabel(['Phase @ ' num2str(FOI(pairs(p,1))) ' Hz (rad)']);
%           title(['Phase @ ' num2str(FOI(pairs(p,2))) ' Hz   vs   phase @ ' num2str(FOI(pairs(p,1))) ' Hz']);
%         end
%         for i = 1:numel(individualPhaseRelationSpk{u}(pairs(p,1),:))
%           plot(individualPhaseRelationSpk{u}(pairs(p,2),i),individualPhaseRelationSpk{u}(pairs(p,1),i), '.k', 'MarkerSize',10);
%           hold on
%         end
        [R, P] = corrcoef(individualPhaseRelationSpk{u}(pairs(p,1),:),individualPhaseRelationSpk{u}(pairs(p,2),:));
        r{u}(p,1) = R(1,2);
        r{u}(p,2) = r{u}(p,1)^2;
        r{u}(p,3) = P(1,2);
%         figure(p)
%         hold on
%         plot(r{u}(p,1), '.', 'MarkerSize',10)
%         hold off
        [RD, PD] = corrcoef(abs(individualPhaseRelationSpk{u}(pairs(p,1),:)),abs(individualPhaseRelationSpk{u}(pairs(p,2),:)));
        rd{u}(p,1) = RD(1,2);
        rd{u}(p,2) = rd{u}(p,1)^2;
        rd{u}(p,3) = PD(1,2);
%         figure(p)
%         hold on
%         plot(rd{u}(p,1), '.', 'MarkerSize',10)
%         hold off
      end
      
%       % Infra-slow vs the rest
%       for f = 2:numel(FOI)
%         i0unit{f-1} = i0full{u}{f};
%         t0unit{f-1} = t0full{u}{f};
%         phasePopTemp(f-1,:) = phasePop(1,:);
%       end
%       [~,meanPhaseRelationIS(:,u),MPRISCu{u},MPRISCl{u},histPhaseRelationIS{u},~] = phaseRelations(i0unit,t0unit,t,phasePopTemp,edges,edgesC);
%       clear i0unit t0unit phasePopTemp
%       
%       % Slow vs delta
%       i0unit{1} = i0full{u}{5};
%       t0unit{1} = t0full{u}{5};
%       phasePopTemp = phasePop(3,:);
%       [~,meanPhaseRelationSD(:,u),MPRSDCu{u},MPRSDCl{u},histPhaseRelationSD{u},~] = phaseRelations(i0unit,t0unit,t,phasePopTemp,edges,edgesC);
%       clear i0unit t0unit phasePopTemp
%       
%       % Slow vs spindle
%       i0unit{1} = i0full{u}{6};
%       t0unit{1} = t0full{u}{6};
%       phasePopTemp = phasePop(3,:);
%       [~,meanPhaseRelationSS(:,u),MPRSSCu{u},MPRSSCl{u},histPhaseRelationSS{u},~] = phaseRelations(i0unit,t0unit,t,phasePopTemp,edges,edgesC);
%       clear i0unit t0unit phasePopTemp
%       
%       % Spindle vs fast
%       i0unit{1} = i0full{u}{7};
%       t0unit{1} = t0full{u}{7};
%       phasePopTemp = phasePop(6,:);
%       [~,meanPhaseRelationSF(:,u),MPRSFCu{u},MPRSFCl{u},histPhaseRelationSF{u},~] = phaseRelations(i0unit,t0unit,t,phasePopTemp,edges,edgesC);
%       clear i0unit t0unit phasePopTemp
    end
    
    % Histogram of unit means
    histMeanPhaseRelation = zeros(length(FOI), length(edges));
    for f = 1:numel(FOI)
      histMeanPhaseRelation(f,:) = [sum(isnan(meanPhaseRelation(f,:))) histcounts(meanPhaseRelation(f,:), edges)];
    end
%     histMeanPhaseRelationIS = zeros(length(FOI)-1, length(edges));
%     for f = 2:numel(FOI)
%       histMeanPhaseRelationIS(f-1,:) = [sum(isnan(meanPhaseRelationIS(f-1,:))) histcounts(meanPhaseRelationIS(f-1,:), edges)];
%     end
%     histMeanPhaseRelationSD = [sum(isnan(meanPhaseRelationSD)) histcounts(meanPhaseRelationSD, edges)];
%     histMeanPhaseRelationSS = [sum(isnan(meanPhaseRelationSS)) histcounts(meanPhaseRelationSS, edges)];
%     histMeanPhaseRelationSF = [sum(isnan(meanPhaseRelationSF)) histcounts(meanPhaseRelationSF, edges)];

    histMeanPhaseRelationSpk = zeros(length(FOI), length(edges));
    for f = 1:numel(FOI)
      histMeanPhaseRelationSpk(f,:) = [sum(isnan(meanPhaseRelationSpk(f,:))) histcounts(meanPhaseRelationSpk(f,:), edges)];
    end
    
% STORE DATA
    dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.FOI = FOI;'];
    eval(dataString);
    dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.Fpass = Fpass;'];
    eval(dataString);
    dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.edges = edges;'];
    eval(dataString);
    dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.phasePairsHT = pairs;'];
    eval(dataString);
    
    dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.histPhaseRelation = histPhaseRelation;'];
    eval(dataString);
    dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.meanPhaseRelation = meanPhaseRelation;'];
    eval(dataString);
    dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.MPRCu = MPRCu;'];
    eval(dataString);
    dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.MPRCl = MPRCl;'];
    eval(dataString);
    dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.histMeanPhaseRelation = histMeanPhaseRelation;'];
    eval(dataString);
    
%     dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.histPhaseRelationIS = histPhaseRelationIS;'];
%     eval(dataString);
%     dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.meanPhaseRelationIS = meanPhaseRelationIS;'];
%     eval(dataString);
%     dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.MPRISCu = MPRISCu;'];
%     eval(dataString);
%     dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.MPRISCl = MPRISCl;'];
%     eval(dataString);
%     dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.histMeanPhaseRelationIS = histMeanPhaseRelationIS;'];
%     eval(dataString);
%     
%     dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.histPhaseRelationSD = histPhaseRelationSD;'];
%     eval(dataString);
%     dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.meanPhaseRelationSD = meanPhaseRelationSD;'];
%     eval(dataString);
%     dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.MPRSDCu = MPRSDCu;'];
%     eval(dataString);
%     dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.MPRSDCl = MPRSDCl;'];
%     eval(dataString);
%     dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.histMeanPhaseRelationSD = histMeanPhaseRelationSD;'];
%     eval(dataString);
%     
%     dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.histPhaseRelationSS = histPhaseRelationSS;'];
%     eval(dataString);
%     dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.meanPhaseRelationSS = meanPhaseRelationSS;'];
%     eval(dataString);
%     dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.MPRSSCu = MPRSSCu;'];
%     eval(dataString);
%     dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.MPRSSCl = MPRSSCl;'];
%     eval(dataString);
%     dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.histMeanPhaseRelationSS = histMeanPhaseRelationSS;'];
%     eval(dataString);
%     
%     dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.histPhaseRelationSF = histPhaseRelationSF;'];
%     eval(dataString);
%     dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.meanPhaseRelationSF = meanPhaseRelationSF;'];
%     eval(dataString);
%     dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.MPRSFCu = MPRSFCu;'];
%     eval(dataString);
%     dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.MPRSFCl = MPRSFCl;'];
%     eval(dataString);
%     dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.histMeanPhaseRelationSF = histMeanPhaseRelationSF;'];
%     eval(dataString);

    dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.histPhaseRelationSpk = histPhaseRelationSpk;'];
    eval(dataString);
    dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.meanPhaseRelationSpk = meanPhaseRelationSpk;'];
    eval(dataString);
    dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.MPRSpkCu = MPRSpkCu;'];
    eval(dataString);
    dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.MPRSpkCl = MPRSpkCl;'];
    eval(dataString);
    dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.histMeanPhaseRelationSpk = histMeanPhaseRelationSpk;'];
    eval(dataString);
    
    dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.PLV = PLV;'];
    eval(dataString);
    dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.MI = MI;'];
    eval(dataString);
    dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.MIRange = MIRange;'];
    eval(dataString);
    
    dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.r = r;'];
    eval(dataString);
    dataString = ['dataStruct.' entryName '.' shankIDs{sh} '.rd = rd;'];
    eval(dataString);
    
    save(dataFile,'dataStruct');
  end
end



function d = designFilter(Fpass, Astop1, Apass, Astop2, Fs)
% A helper function of AnHT for designing Buterworth filters.

d = designfilt('bandpassiir', ...
	'StopbandFrequency1',Fpass(1,1), 'PassbandFrequency1',Fpass(2,1), ...
  'PassbandFrequency2',Fpass(4,1), 'StopbandFrequency2',Fpass(5,1), ...
  'StopbandAttenuation1',Astop1, 'PassbandRipple',Apass, ...
  'StopbandAttenuation2',Astop2, 'DesignMethod','butter', 'SampleRate',Fs);
end


function d = designFilterLO(Fpass, Astop1, Apass, Astop2, Fs) %#ok<DEFNU>
% A helper function of AnHT for designing low order filters.

d = designfilt('bandpassiir', ...
  'FilterOrder',2, ...
  'PassbandFrequency1', Fpass(2,1),'PassbandFrequency2', Fpass(4,1), ...
  'StopbandAttenuation1', Astop1, 'PassbandRipple', Apass, ...
  'StopbandAttenuation2', Astop2, ...
  'SampleRate', Fs);
end


function [t0full, i0full] = zeroPhase(FOI, t, phase)
% A helper function of AnHT for obtaining zero phase time indeces given the
% phase vector.

zeroPhaseTrans = phase;
zeroPhaseTrans(zeroPhaseTrans > 0) = 0;
zeroPhaseTrans(zeroPhaseTrans < 0) = 1;
for f = 1:numel(FOI)
  [~, zeroPhaseInd] = findpeaks(zeroPhaseTrans(f,:));
  t0 = zeros(size(zeroPhaseInd));
  for i = 1:numel(zeroPhaseInd)
    tZeroPhase = [t(zeroPhaseInd(i)-1) t(zeroPhaseInd(i))];
    phaseZero = [phase(f,zeroPhaseInd(i)-1) phase(f,zeroPhaseInd(i))];
    slope = (phaseZero(2)-phaseZero(1))/(tZeroPhase(2)-tZeroPhase(1));
    t0(i) = t(zeroPhaseInd(i)-1) + phase(f,zeroPhaseInd(i)-1)/abs(slope);
  end
  t0full{f} = t0; %#ok<AGROW>
  i0full{f} = zeroPhaseInd; %#ok<AGROW>
end
end


function [phaseRel, meanPhaseRel, MPRCu, MPRCl, histPhaseRel, histPhaseRelC] = phaseRelations(i0unit,t0unit,t,phasePop,edges,edgesC)
% A helper function of AnHT for calculating phase relations between a unit
% and the population rate at narrow frequency ranges: the Hilbert transform
% of the unit rate vs the Hilbert transform of the population rate.

nF = numel(i0unit);
meanPhaseRel = zeros(nF,1);
MPRCu = zeros(size(meanPhaseRel));
MPRCl = zeros(size(meanPhaseRel));
histPhaseRel = zeros(nF, length(edges)-1);
histPhaseRelC = zeros(nF, length(edgesC)-1);
for f = 1:nF
  
  % Single cycle phase relations
  for i = 1:numel(i0unit{f})
    ind = i0unit{f}(i);
    phaseMatch = phasePop(f,ind);
    phaseMatchPre = phasePop(f,ind-1);
    tMatch = t(ind);
    tMatchPre = t(ind-1);
    tMatchActual = t0unit{f}(i);
    phaseRel{f}(i) = phaseMatchPre + (phaseMatch-phaseMatchPre)*((tMatchActual-tMatchPre)/(tMatch-tMatchPre)); %#ok<AGROW>
  end
  
  % Single unit phase relations
  %histogram(phaseRelation, edgesC);
  histPhaseRel(f,:) = histcounts(phaseRel{f}, edges);
  histPhaseRelC(f,:) = histcounts(phaseRel{f}, edgesC);
  [meanPhaseRel(f), MPRCu(f), MPRCl(f)] = circ_mean(phaseRel{f}');
  if isnan(MPRCu(f))
    meanPhaseRel(f) = NaN;
  end
end
end


function [phaseRel, meanPhaseRel, MPRCu, MPRCl, histPhaseRel, histPhaseRelC] = phaseRelationsSpk(spk,phasePop,edges,edgesC)
% A helper function of AnHT for calculating phase relations between a unit
% and the population rate at narrow frequency ranges: Unit rate vs the
% Hilbert transform of the population rate.

nF = size(phasePop,1);
meanPhaseRel = zeros(nF,1);
MPRCu = zeros(size(meanPhaseRel));
MPRCl = zeros(size(meanPhaseRel));
histPhaseRel = zeros(nF, length(edges)-1);
histPhaseRelC = zeros(nF, length(edgesC)-1);
phaseRel = zeros(nF,sum(spk));
for f = 1:nF
  
  % Single cycle phase relations
  ind = 1:length(spk);
  ind = ind(logical(spk));
  indSpk = spk(ind);
  countSpk = 1;
  for i = 1:length(ind)
    for s = 1:indSpk(i)
      phaseRel(f,countSpk) = phasePop(f, ind(i));
      countSpk = countSpk + 1;
    end
  end
  
  % Single unit phase relations
  %histogram(phaseRelation, edgesC);
  histPhaseRel(f,:) = histcounts(phaseRel(f,:), edges);
  histPhaseRelC(f,:) = histcounts(phaseRel(f,:), edgesC);
  [meanPhaseRel(f), MPRCu(f), MPRCl(f)] = circ_mean(phaseRel(f,:)');
  if isnan(MPRCu(f))
    meanPhaseRel(f) = NaN;
  end
end
end


function [R, RRange] = phaseLocking(spk,phasePop)
% A helper function of AnHT for calculating phase locking between a unit
% and the population rate at narrow frequency ranges: Unit rate vs the
% Hilbert transform of the population rate.

nF = size(phasePop,1);
phaseRel = zeros(nF,sum(spk));
for f = 1:nF
  
  % Single cycle phase relations
  ind = 1:length(spk);
  ind = ind(logical(spk));
  indSpk = spk(ind);
  countSpk = 1;
  for i = 1:length(ind)
    for s = 1:indSpk(i)
      phaseRel(f,countSpk) = phasePop(f, ind(i));
      countSpk = countSpk + 1;
    end
  end
end

R = abs(sum(exp(1j*phaseRel),2)./(countSpk-1));

% Calculate range
RRange = [0 1];
end


function [MI, MIRange] = modulationIndex(histPhaseRel)
% A helper function of AnHT for calculating phase modulation index between
% a unit and the population rate at narrow frequency ranges: Unit rate vs
% the Hilbert transform of the population rate.

% Normalise
for f = 1:size(histPhaseRel,1)
  histPhaseRel(f,:) = histPhaseRel(f,:)./sum(histPhaseRel(f,:));
end
histPhaseRel(histPhaseRel == 0) = 1e-9;

% Calculate MI
entropy = -sum(histPhaseRel.*log(histPhaseRel),2);
KLDistance = log(size(histPhaseRel,2)) - entropy;
MI = KLDistance./log(size(histPhaseRel,2));

% Calculate range
MIRange = [0 log(size(histPhaseRel,2))];
end