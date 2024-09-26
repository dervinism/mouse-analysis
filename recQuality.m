% Run this script to examine the recording quality


fclose all;
close all
clear
clc


dataFile = 'R:\CSN\Shared\Dynamics\Data\mouse_analysis_data\M190523_A_MD\noCAR\M190523_A_MD.mat';

cx1Series = [1 10 19 28 37 46 55 64 73];
cx2Series = cx1Series+8;
hpSeries = cx1Series+7;
thSeries = cx1Series+4;

qualitySeries_cx1 = {'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906011331041\continuous_probe1_swappedNoCAR.qua.1.mat';
                     'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906022028431\continuous_probe1_swappedNoCAR.qua.1.mat';
                     'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906041747361\continuous_probe1_swappedNoCAR.qua.1.mat';
                     'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906042002581\continuous_probe1_swappedNoCAR.qua.1.mat';
                     'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906061429121\continuous_probe1_swappedNoCAR.qua.1.mat';
                     'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906061638161\continuous_probe1_swappedNoCAR.qua.1.mat';
                     'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906171841251\continuous_probe1_swappedNoCAR.qua.1.mat';
                     'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906172049561\continuous_probe1_swappedNoCAR.qua.1.mat';
                     'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906211143371\continuous_probe1_swappedNoCAR.qua.1.mat'};
qualitySeries_cx2 = {'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906011331042\continuous_probe2_swappedNoCAR.qua.1.mat';
                     'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906022028432\continuous_probe2_swappedNoCAR.qua.1.mat';
                     'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906041747362\continuous_probe2_swappedNoCAR.qua.1.mat';
                     'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906042002582\continuous_probe2_swappedNoCAR.qua.1.mat';
                     'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906061429122\continuous_probe2_swappedNoCAR.qua.1.mat';
                     'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906061638162\continuous_probe2_swappedNoCAR.qua.1.mat';
                     'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906171841252\continuous_probe2_swappedNoCAR.qua.1.mat';
                     'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906172049562\continuous_probe2_swappedNoCAR.qua.1.mat';
                     'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906211143372\continuous_probe2_swappedNoCAR.qua.1.mat'};
qualitySeries_hp = {'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906011331042\continuous_probe2_swappedNoCAR.qua.1.mat';
                    'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906022028432\continuous_probe2_swappedNoCAR.qua.1.mat';
                    'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906041747362\continuous_probe2_swappedNoCAR.qua.1.mat';
                    'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906042002582\continuous_probe2_swappedNoCAR.qua.1.mat';
                    'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906061429122\continuous_probe2_swappedNoCAR.qua.1.mat';
                    'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906061638162\continuous_probe2_swappedNoCAR.qua.1.mat';
                    'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906171841252\continuous_probe2_swappedNoCAR.qua.1.mat';
                    'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906172049562\continuous_probe2_swappedNoCAR.qua.1.mat';
                    'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906211143372\continuous_probe2_swappedNoCAR.qua.1.mat'};
qualitySeries_th = {'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906011331042\continuous_probe2_swappedNoCAR.qua.1.mat';
                    'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906022028432\continuous_probe2_swappedNoCAR.qua.1.mat';
                    'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906041747362\continuous_probe2_swappedNoCAR.qua.1.mat';
                    'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906042002582\continuous_probe2_swappedNoCAR.qua.1.mat';
                    'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906061429122\continuous_probe2_swappedNoCAR.qua.1.mat';
                    'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906061638162\continuous_probe2_swappedNoCAR.qua.1.mat';
                    'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906171841252\continuous_probe2_swappedNoCAR.qua.1.mat';
                    'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906172049562\continuous_probe2_swappedNoCAR.qua.1.mat';
                    'R:\CSN\Shared\Dynamics\Data\M190523_A_MD\201906211143372\continuous_probe2_swappedNoCAR.qua.1.mat'};

recDates = {'', '19-06-01', '19-06-02', '19-06-04:1', '19-06-04:2', '19-06-06:1', '19-06-06:2', '19-06-17:1', '19-06-17:2', '19-06-21', ''};

separators = [6.5];


load(dataFile)


[nUnitsDB_cx1, nGoodUnitsDB_cx1, mfrDB_cx1, mfrDB_mean_cx1, mfrMUADB_cx1] = qualityMeasures(dataStruct, cx1Series, qualitySeries_cx1);
[nUnitsDB_hp, nGoodUnitsDB_hp, mfrDB_hp, mfrDB_mean_hp, mfrMUADB_hp] = qualityMeasures(dataStruct, hpSeries, qualitySeries_hp);
[nUnitsDB_th, nGoodUnitsDB_th, mfrDB_th, mfrDB_mean_th, mfrMUADB_th] = qualityMeasures(dataStruct, thSeries, qualitySeries_th);
if ~isempty(cx2Series)
  [nUnitsDB_cx2, nGoodUnitsDB_cx2, mfrDB_cx2, mfrDB_mean_cx2, mfrMUADB_cx2] = qualityMeasures(dataStruct, cx2Series, qualitySeries_cx2);
  if numel(nUnitsDB_cx1) < numel(nUnitsDB_th)
    diff = abs(numel(nUnitsDB_cx1) - numel(nUnitsDB_th));
    nUnitsDB = [nUnitsDB_cx1 zeros(1,diff)] + nUnitsDB_cx2 + nUnitsDB_hp + nUnitsDB_th;
    nGoodUnitsDB = [nGoodUnitsDB_cx1 zeros(1,diff)] + nGoodUnitsDB_cx2 + nGoodUnitsDB_hp + nGoodUnitsDB_th;
  elseif numel(nUnitsDB_cx1) > numel(nUnitsDB_th)
    diff = abs(numel(nUnitsDB_cx1) - numel(nUnitsDB_th));
    nUnitsDB = nUnitsDB_cx1 + [nUnitsDB_cx2 zeros(1,diff)] + [nUnitsDB_hp zeros(1,diff)] + [nUnitsDB_th zeros(1,diff)];
    nGoodUnitsDB = nGoodUnitsDB_cx1 + [nGoodUnitsDB_cx2 zeros(1,diff)] + [nGoodUnitsDB_hp zeros(1,diff)] + [nGoodUnitsDB_th zeros(1,diff)];
  else
    nUnitsDB = nUnitsDB_cx1 + nUnitsDB_cx2 + nUnitsDB_hp + nUnitsDB_th;
    nGoodUnitsDB = nGoodUnitsDB_cx1 + nGoodUnitsDB_cx2 + nGoodUnitsDB_hp + nGoodUnitsDB_th;
  end
else
  nUnitsDB_cx2 = []; nGoodUnitsDB_cx2 = []; mfrDB_cx2 = []; mfrDB_mean_cx2 = []; mfrMUADB_cx2 = [];
  nUnitsDB = nUnitsDB_cx1 + nUnitsDB_hp + nUnitsDB_th;
  nGoodUnitsDB = nGoodUnitsDB_cx1 + nGoodUnitsDB_hp + nGoodUnitsDB_th;
end


colours = [000 000 255;       % blue
           000 255 255;       % cyan
           000 000 000;       % black
           128 128 128;       % grey
           255 000 000;       % red
           255 000 255;       % magenta
           000 128 000;       % green
           000 255 000]./255; % lime

fig1 = figure('units', 'normalized', 'position', [0.002, .04, 1, .88]);
for dbCount = 1:numel(nUnitsDB)
  plot(dbCount,nUnitsDB(dbCount)+0.05, '.', 'Color',colours(3,:), 'MarkerSize',20)
  if dbCount == 1
    hold on
  end
  plot(dbCount,nGoodUnitsDB(dbCount)+0.05, 'o', 'Color',colours(3,:), 'MarkerSize',10)
  try
    plot(dbCount-0.05,nUnitsDB_cx1(dbCount), '.', 'Color',colours(1,:), 'MarkerSize',20)
    plot(dbCount-0.05,nGoodUnitsDB_cx1(dbCount), 'o', 'Color',colours(1,:), 'MarkerSize',10)
  catch me
  end
  try
    plot(dbCount,nUnitsDB_hp(dbCount), '.', 'Color',colours(5,:), 'MarkerSize',20)
    plot(dbCount,nGoodUnitsDB_hp(dbCount), 'o', 'Color',colours(5,:), 'MarkerSize',10)
  catch me
  end
  try
    plot(dbCount+0.05,nUnitsDB_th(dbCount), '.', 'Color',colours(8,:), 'MarkerSize',20)
    plot(dbCount+0.05,nGoodUnitsDB_th(dbCount), 'o', 'Color',colours(8,:), 'MarkerSize',10)
  catch me
  end
  if ~isempty(cx2Series)
    try
      plot(dbCount,nUnitsDB_cx2(dbCount)-0.05, '.', 'Color',colours(2,:), 'MarkerSize',20)
      plot(dbCount,nGoodUnitsDB_cx2(dbCount)-0.05, 'o', 'Color',colours(2,:), 'MarkerSize',10)
    catch me
    end
  end
end
for s = 1:numel(separators)
  plot([separators(s) separators(s)],[0 max(nUnitsDB)+1], ':', 'Color',colours(4,:))
end
plot(separators,zeros(size(separators)), '.w', 'MarkerSize',20);
hold off
xlim([0 numel(cx1Series)+1]);
ylim([0 max(nUnitsDB)+1]);
title(['Unit count for animal ' dataFile(1:end-4)], 'Interpreter','None');
xlabel('Recording session');
ylabel('Unit count');
if ~isempty(cx2Series)
  legend('All','All quality', 'S1','S1 quality', 'Hp','Hp quality', 'Th','Th quality', 'RSC','RSC quality');
else
  legend('All','All quality', 'S1','S1 quality', 'Hp','Hp quality', 'Th','Th quality');
end
legend boxoff
xticks([0 1:numel(cx1Series)+1]);
xticklabels(recDates);
set(gcf,'color','w')
ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 18, 4/3, 2, [0.005 0], 'out',...
    'on', 'k', {'Recording session'}, [0 numel(cx1Series)+1], [0 1:numel(cx1Series)+1],...
    'on', 'k', {'Pupil area (a.u.)'}, [0 max(nUnitsDB)+1], 0:5:max(nUnitsDB)+1);
figName = ['Unit_count__' dataFile(1:end-4)];
set(fig1, 'Name',figName);
hgsave(fig1, figName);
close(fig1);

fig2 = figure('units', 'normalized', 'position', [0.002, .04, 1, .88]);
plot(1:numel(cx1Series),mfrMUADB_cx1, '-o', 'Color',colours(1,:), 'MarkerSize',5, 'MarkerFaceColor',colours(1,:))
hold on
plot(1:numel(cx1Series),mfrDB_mean_cx1, ':o', 'Color',colours(1,:), 'MarkerSize',10)
plot(1:numel(hpSeries),mfrMUADB_hp, '-o', 'Color',colours(5,:), 'MarkerSize',5, 'MarkerFaceColor',colours(5,:))
plot(1:numel(hpSeries),mfrDB_mean_hp, ':o', 'Color',colours(5,:), 'MarkerSize',10)
plot(1:numel(thSeries),mfrMUADB_th, '-o', 'Color',colours(8,:), 'MarkerSize',5, 'MarkerFaceColor',colours(8,:))
plot(1:numel(thSeries),mfrDB_mean_th, ':o', 'Color',colours(8,:), 'MarkerSize',10)
if ~isempty(cx2Series)
  plot(1:numel(cx2Series),mfrMUADB_cx2, '-o', 'Color',colours(2,:), 'MarkerSize',5, 'MarkerFaceColor',colours(2,:))
  plot(1:numel(cx2Series),mfrDB_mean_cx2, ':o', 'Color',colours(2,:), 'MarkerSize',10)
end
yRange = ylim;
for s = 1:numel(separators)
  plot([separators(s) separators(s)],[0.4 yRange(2)-0.4], 'w', 'LineWidth',90)
  plot([separators(s) separators(s)],[0 yRange(2)], ':', 'Color',colours(4,:))
end
plot(separators,zeros(size(separators)), '.w', 'MarkerSize',20);
hold off
xlim([0 numel(cx1Series)+1]);
ylim(yRange);
title(['Firing rates for animal ' dataFile(1:end-4)], 'Interpreter','None');
xlabel('Recording session');
ylabel('Firing rate / Mean firing rate (APs/second)');
if ~isempty(cx2Series)
  legend('Cx_1 MUA','Cx_1 units', 'Hp MUA','Hp units', 'Th MUA','Th units', 'Cx_2 MUA','Cx_2 units');
else
  legend('Cx_1 MUA','Cx_1 units', 'Hp MUA','Hp units', 'Th MUA','Th units');
end
xticks([0 1:numel(cx1Series)+1]);
xticklabels(recDates);
figName = ['Firing_rates__' dataFile(1:end-4)];
set(fig2, 'Name',figName);
hgsave(fig2, figName);
close(fig2);



function [nUnitsDB, nGoodUnitsDB, mfrDB, mfrDB_mean, mfrMUADB] = qualityMeasures(dataStruct, series, qualitySeries)

cluDist = 20;
refractCont = 0.2;

seriesCount = 0;
nUnitsDB = [];
nGoodUnitsDB = [];
mfrMUADB = [];
mfrDB_mean = [];

% LOOP THROUGH DB ENTRIES
fnsData = fieldnames(dataStruct.seriesData);
for dbCount = series
  seriesCount = seriesCount + 1;
  dbStruct = dataStruct.seriesData.(fnsData{dbCount});
  shankIDs = fieldnames(dbStruct.shankData);
  
% LOAD THE UNIT QUALITY FILE
  if iscell(qualitySeries{seriesCount})
    unitQtemp = [];
    for iQS = 1:numel(qualitySeries{seriesCount})
      if ~isempty(qualitySeries{seriesCount}{iQS})
        load(qualitySeries{seriesCount}{iQS});
        unitQ(:,1) = unitQ(:,1)+1000*(iQS-1);
        unitQtemp = [unitQtemp; unitQ];
      end
    end
    unitQ = unitQtemp;
  else
    if ~isempty(qualitySeries{seriesCount})
      load(qualitySeries{seriesCount});
    else
      unitQ = [];
    end
  end
  
  mfr = [];
  unitCount = 0;
  goodUnitCount = 0;
  goodUnits = [];
  
% LOOP THROUGH SHANKS
  for sh = 1:numel(shankIDs)
    shankStruct = eval(['dbStruct.shankData.' shankIDs{sh}]);
    fprintf('%s shank %d -------------------\n', fnsData{dbCount}, sh);
    
% LOAD THE CONTENTS OF THE SHANK STRUCTURE VARIABLE
    units = shankStruct.units;
    if isfield(shankStruct, 'phaseCoh')
      phaseCoh = shankStruct.phaseCoh;
    else
      mfrMUA = 0;
      break
    end
    
    for u = 1:numel(units)
      fprintf('Started processing unit %i\n',units(u));
      
      if ~isempty(unitQ)
        ind = find(units(u) == unitQ(:,1)); %#ok<*IDISVAR,*NODEF>
        if ind > size(unitQ,1)
          ind = ind - size(unitQ,1);
        end
        if unitQ(ind,2) >= cluDist && unitQ(ind,6) <= refractCont % quality check
          goodUnits = [goodUnits; units(u)]; %#ok<*AGROW>
        end
      end
      
      unitData = phaseCoh{u};
      mfr = [mfr; unitData.mfr]; %#ok<*AGROW>
      if sh == 1 && u == 1
        if isfield(unitData, 'mfrMUA')
          mfrMUA = unitData.mfr_mua;
        else
          mfrMUA = 0;
        end
      end
    end
    
    unitCount = unitCount + numel(units);
    goodUnitCount = goodUnitCount + numel(goodUnits);
  end
  
% SUMMARISE VARIABLES
  nUnitsDB = [nUnitsDB unitCount];
  nGoodUnitsDB = [nGoodUnitsDB goodUnitCount];
  mfrDB{seriesCount} = mfr; %#ok<*SAGROW>
  mfrDB_mean = [mfrDB_mean mean(mfr)];
  mfrMUADB = [mfrMUADB mfrMUA];
end
end