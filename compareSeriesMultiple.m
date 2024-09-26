% Run this script to perform phase and coherence comparisons between
% vigilance states pooling over multiple series and animals.


% INITIALISE PARAMETERS
visibility = 'on';
dataFileList = {'R:\CSN\Shared\Dynamics\Data\mouse_analysis_data\M190128_A_MD\area_comparisons\ThVsHp_restrained_vs_anaesthesia\M190128_A_MD_sM190128_A_MD_s2019020810040424__M190128_A_MD_s2019020810040456_sM190128_A_MD_s2019020811053624__M190128_A_MD_s2019020811053656.mat';
                'R:\CSN\Shared\Dynamics\Data\mouse_analysis_data\M190128_A_MD\area_comparisons\ThVsHp_restrained_vs_anaesthesia\M190128_A_MD_sM190128_A_MD_s2019021209595024__M190128_A_MD_s2019021209595056_sM190128_A_MD_s2019021211150224__M190128_A_MD_s2019021211150256.mat';
                'R:\CSN\Shared\Dynamics\Data\mouse_analysis_data\M190128_A_MD\area_comparisons\ThVsHp_restrained_vs_anaesthesia\M190128_A_MD_sM190128_A_MD_s2019021511250724__M190128_A_MD_s2019021511250756_sM190128_A_MD_s2019021512353724__M190128_A_MD_s2019021512353756.mat';
                'R:\CSN\Shared\Dynamics\Data\mouse_analysis_data\M190128_B_MD\area_comparisons\ThVsHp_restrained_vs_anaesthesia\M190128_B_MD_sM190128_B_MD_s2019031115075056__M190128_B_MD_s2019031115075024_sM190128_B_MD_s2019031116185956__M190128_B_MD_s2019031116185924.mat';
                'R:\CSN\Shared\Dynamics\Data\mouse_analysis_data\M190128_C_MD\area_comparisons\ThVsHp_restrained_vs_anaesthesia\M190128_C_MD_sM190128_C_MD_s2019022510293924__M190128_C_MD_s2019022510293956_sM190128_C_MD_s2019022511490924__M190128_C_MD_s2019022511490956.mat'};
              
qualityFileList = {{'R:\CSN\Shared\Dynamics\Data\M190128_A_MD\201902081004042\continuous_probe2_swappedNoCAR.qua.1.mat';
                    'R:\CSN\Shared\Dynamics\Data\M190128_A_MD\201902081004043\continuous_probe3_swappedNoCAR.qua.1.mat';
                    'R:\CSN\Shared\Dynamics\Data\M190128_A_MD\201902081004044\continuous_probe4_swappedNoCAR.qua.1.mat'};
                   {'R:\CSN\Shared\Dynamics\Data\M190128_A_MD\201902120959502\continuous_probe2_swappedNoCAR.qua.1.mat';
                    'R:\CSN\Shared\Dynamics\Data\M190128_A_MD\201902120959503\continuous_probe3_swappedNoCAR.qua.1.mat';
                    'R:\CSN\Shared\Dynamics\Data\M190128_A_MD\201902120959504\continuous_probe4_swappedNoCAR.qua.1.mat'};
                   {'R:\CSN\Shared\Dynamics\Data\M190128_A_MD\201902121115022\continuous_probe2_swappedNoCAR_zeroedOut.qua.1.mat';
                    'R:\CSN\Shared\Dynamics\Data\M190128_A_MD\201902121115023\continuous_probe3_swappedNoCAR_zeroedOut.qua.1.mat';
                    'R:\CSN\Shared\Dynamics\Data\M190128_A_MD\201902121115024\continuous_probe4_swappedNoCAR_zeroedOut.qua.1.mat'};
                    'R:\CSN\Shared\Dynamics\Data\M190128_B_MD\2019031115075026\continuous_probe26_swappedNoCAR.qua.1.mat';
                    'R:\CSN\Shared\Dynamics\Data\M190128_C_MD\2019022510293926\continuous_probe26_swappedNoCAR.qua.1.mat'};

chColours = {{1:6; 7:12; 13:16}; {1:6; 7:12; 13:16}; {1:6; 7:12; 13:16}; {1:6; 7:12; 13:16}; {1:4; 5:9; 10:15}};
chAreas = {'VB'; 'Po'; 'LPM'};

% cluDist = 20;
% refractCont = 0.2;
cluDist = 0;
refractCont = inf;

% LOOP THROUGH FILES AND LOAD DATA
for iFile = 1:numel(dataFileList)
  
  if iscell(qualityFileList{iFile})
    unitQtemp = [];
    for iQS = 1:numel(qualityFileList{iFile})
      if ~isempty(qualityFileList{iFile}{iQS})
        load(qualityFileList{iFile}{iQS});
        unitQ(:,1) = unitQ(:,1)+1000*(iQS-1);
        unitQtemp = [unitQtemp; unitQ]; %#ok<*AGROW>
      end
    end
    unitQ = unitQtemp;
  else
    load(qualityFileList{iFile});
  end
  
  load(dataFileList{iFile});
  if iFile == 1
    for iFreq = 1:numel(FOI)
      phaseVec{iFreq} = [];
      cohVec{iFreq} = [];
      cohConfVec{iFreq} = [];
      uCentres{iFreq} = [];
      uSeries{iFreq} = [];
    end
  end
  for iFreq = 1:numel(FOI)
    for u = 1:numel(phaseCohSeries)
      if isfield(phaseCohSeries{u}, 'awake') && isfield(phaseCohSeries{u}, 'anaesthesia')
        ind = find(u == unitQ(:,1)); %#ok<*IDISVAR,*NODEF>
        if ind > size(unitQ,1)
          ind = ind - size(unitQ,1);
        end
        unitQwhich = unitQ;
        if unitQwhich(ind,2) >= cluDist && unitQwhich(ind,6) <= refractCont % quality check
          phaseVec{iFreq} = [phaseVec{iFreq} [phaseCohSeries{u}.awake.phase(iFreq);...
            phaseCohSeries{u}.anaesthesia.phase(iFreq)]];
          cohVec{iFreq} = [cohVec{iFreq} [phaseCohSeries{u}.awake.coh(iFreq);...
            phaseCohSeries{u}.anaesthesia.coh(iFreq)]];
          cohConfVec{iFreq} = [cohConfVec{iFreq} [phaseCohSeries{u}.awake.coh_conf(iFreq);...
            phaseCohSeries{u}.anaesthesia.coh_conf(iFreq)]];
          if isfield(phaseCohSeries{u}.awake,'uCentre')
            uCentre = phaseCohSeries{u}.awake.uCentre;
            if iscell(qualityFileList{iFile})
              if uCentre < 1000
                uCentre = chColours{iFile}{1}(uCentre);
              elseif uCentre < 2000
                uCentre = chColours{iFile}{2}(uCentre-1000);
              elseif uCentre < 3000
                uCentre = chColours{iFile}{3}(uCentre-2000);
              end
            end
            uCentres{iFreq} = [uCentres{iFreq} uCentre];
            uSeries{iFreq} = [uSeries{iFreq} iFile];
          end
        end
      end
    end
  end
end


% INITIALISE FIGURES
for j = 1:numel(FOI)
  figPhase(j) = figure('units', 'normalized', 'position', [0.002, .04, 1/2, .88/1.04], 'Visible', visibility); %#ok<*SAGROW>
  title(['Phase: ' num2str(FOI(j)) ' Hz']);
  xlabel('Wakefulness phase (rad)')
  ylabel('Anaesthesia phase (rad)')
  figCoherence(j) = figure('units', 'normalized', 'position', [0.002, .04, 1/2, .88/1.04], 'Visible', visibility);
  title(['Coherence: ' num2str(FOI(j)) ' Hz']);
  xlabel('Wakefulness coherence')
  ylabel('Anaesthesia coherence')
end


% DRAW FIGURES
follow = false;
for i = 1:numel(figPhase)
  for j = 1:numel(phaseVec{i}(1,:))
    if ~sum(isempty(phaseVec{i})) && ~sum(isempty(phaseVec{i}))
      figure(figPhase(i));
      hold on
      if isempty(uCentres{i})
        plot(phaseVec{i}(1,j),phaseVec{i}(2,j), '.', 'MarkerSize',20);
      else
        for k = 1:numel(chColours{uSeries(j)})
          if sum(uCentres{i}(j) == chColours{uSeries(j)}{k})
            if k == 1
              p1{i} = plot(phaseVec{i}(1,j),phaseVec{i}(2,j), 'r.', 'MarkerSize',20);
            elseif k == 2
              p2{i} = plot(phaseVec{i}(1,j),phaseVec{i}(2,j), 'b.', 'MarkerSize',20);
            else
              p3{i} = plot(phaseVec{i}(1,j),phaseVec{i}(2,j), 'g.', 'MarkerSize',20);
            end
            break
          end
        end
      end
      hold off
      follow = true;
    end
    if ~sum(isempty(cohVec{i})) && ~sum(isempty(cohVec{i}))
      figure(figCoherence(i));
      hold on
      if isempty(uCentres{i})
        plot(cohVec{i}(1,j),cohVec{i}(2,j), '.', 'MarkerSize',20);
      else
        for k = 1:numel(chColours{uSeries(j)})
          if sum(uCentres{i}(j) == chColours{uSeries(j)}{k})
            if k == 1
              p4{i} = plot(cohVec{i}(1,j),cohVec{i}(2,j), 'r.', 'MarkerSize',20);
            elseif k == 2
              p5{i} = plot(cohVec{i}(1,j),cohVec{i}(2,j), 'b.', 'MarkerSize',20);
            else
              p6{i} = plot(cohVec{i}(1,j),cohVec{i}(2,j), 'g.', 'MarkerSize',20);
            end
            break
          end
        end
      end
      hold off
      follow = true;
    end
  end
end
if exist('p1','var')
  p1{numel(figPhase)+1} = [];
end
if exist('p2','var')
  p2{numel(figPhase)+1} = [];
end
if exist('p3','var')
  p3{numel(figPhase)+1} = [];
end
if exist('p4','var')
  p4{numel(figPhase)+1} = [];
end
if exist('p5','var')
  p5{numel(figPhase)+1} = [];
end
if exist('p6','var')
  p6{numel(figPhase)+1} = [];
end
if exist('p1','var')
  for j = 1:numel(figPhase)
    figure(figPhase(j));
    if (exist('p1','var') && ~isempty(p1{j})) && (exist('p2','var') && ~isempty(p2{j})) && (exist('p3','var') && ~isempty(p3{j}))
      legend([p1{j} p2{j} p3{j}], chAreas{1}, chAreas{2}, chAreas{3}, 'Location','northwest');
    elseif (exist('p1','var') && ~isempty(p1{j})) && (exist('p2','var') && ~isempty(p2{j}))
      legend([p1{j} p2{j}], chAreas{1}, chAreas{2}, 'Location','northwest');
    elseif (exist('p1','var') && ~isempty(p1{j})) && (exist('p3','var') && ~isempty(p3{j}))
      legend([p1{j} p3{j}], chAreas{1}, chAreas{3}, 'Location','northwest');
    elseif (exist('p2','var') && ~isempty(p2{j})) && (exist('p3','var') && ~isempty(p3{j}))
      legend([p2{j} p3{j}], chAreas{2}, chAreas{3}, 'Location','northwest');
    elseif (exist('p1','var') && ~isempty(p1{j}))
      legend(p1{j}, chAreas{1}, 'Location','northwest');
    elseif (exist('p2','var') && ~isempty(p2{j}))
      legend(p2{j}, chAreas{2}, 'Location','northwest');
    elseif (exist('p3','var') && ~isempty(p3{j}))
      legend(p3{j}, chAreas{3}, 'Location','northwest');
    end
    figure(figCoherence(j));
    if (exist('p4','var') && ~isempty(p4{j})) && (exist('p5','var') && ~isempty(p5{j})) && (exist('p6','var') && ~isempty(p6{j}))
      legend([p4{j} p5{j} p6{j}], chAreas{1}, chAreas{2}, chAreas{3}, 'Location','northwest');
    elseif (exist('p4','var') && ~isempty(p4{j})) && (exist('p5','var') && ~isempty(p5{j}))
      legend([p4{j} p5{j}], chAreas{1}, chAreas{2}, 'Location','northwest');
    elseif (exist('p4','var') && ~isempty(p4{j})) && (exist('p6','var') && ~isempty(p6{j}))
      legend([p4{j} p6{j}], chAreas{1}, chAreas{3}, 'Location','northwest');
    elseif (exist('p5','var') && ~isempty(p5{j})) && (exist('p6','var') && ~isempty(p6{j}))
      legend([p5{j} p6{j}], chAreas{2}, chAreas{3}, 'Location','northwest');
    elseif (exist('p4','var') && ~isempty(p4{j}))
      legend(p4{j}, chAreas{1}, 'Location','northwest');
    elseif (exist('p5','var') && ~isempty(p5{j}))
      legend(p5{j}, chAreas{2}, 'Location','northwest');
    elseif (exist('p6','var') && ~isempty(p6{j}))
      legend(p6{j}, chAreas{3}, 'Location','northwest');
    end
  end
end


% CORRELATION ANALYSES
if follow
  [rPhase, rhoPhase, rhoCircPhase, rCoherence, rhoCoherence, rhoCircCoherence] = corrStates(...
    figPhase, figCoherence, phaseVec, cohVec, cohConfVec);
  
  
% SAVE FIGURES
  initFigName = '_ephys_vigilance_states_combined_';
  saveFigsState(initFigName, figPhase, figCoherence, FOI);
else
  rPhase = []; rhoPhase = []; rhoCircPhase = []; rCoherence = []; rhoCoherence = []; rhoCircCoherence = [];
end


% SAVE DATA
save('vigilance_states_combined.mat','phaseVec','cohVec','cohConfVec','rPhase','rhoPhase','rhoCircPhase',...
  'rCoherence','rhoCoherence','rhoCircCoherence','-v7.3');