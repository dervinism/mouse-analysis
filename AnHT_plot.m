% Run this script to perform plot for Hilbert transform analyses of mouse
% data.


% LOAD PRE-PROCESSED DATA
if ~exist('dataStruct','var')
  load(dataFile);
end


% INITIALISE PARALLEL POOL
% delete(gcp('nocreate'))
% pool = parpool;
% feature('numcores');


% LOOP THROUGH DB ENTRIES
fWidth = 'wide';
visibility = 'on';
histMeanPhaseRelationTotal = [];
fnsData = fieldnames(dataStruct);
significance = 0.05;
for dbCount = numel(fnsData):-1:1
  dbStruct = dataStruct.(fnsData{dbCount});
  shankIDs = fieldnames(dbStruct);

% LOOP THROUGH SHANKS
  for sh = 1:numel(shankIDs)
    shankStruct = dbStruct.(shankIDs{sh});
    fprintf('Processing %s shank %d \n', fnsData{dbCount}, sh);
    
% LOAD THE CONTENTS OF THE SHANK STRUCTURE VARIABLE
    units = shankStruct.units;
    FOI = shankStruct.FOI;
    edges = shankStruct.edges;
    phasePairsHT = shankStruct.phasePairsHT;
    histPhaseRelation = shankStruct.histPhaseRelation;
    meanPhaseRelation = shankStruct.meanPhaseRelation;
    MPRCu = shankStruct.MPRCu;
    MPRCl = shankStruct.MPRCl;
    histMeanPhaseRelation = shankStruct.histMeanPhaseRelation;
%     histPhaseRelationIS = shankStruct.histPhaseRelationIS;
%     meanPhaseRelationIS = shankStruct.meanPhaseRelationIS;
%     MPRISCu = shankStruct.MPRISCu;
%     MPRISCl = shankStruct.MPRISCl;
%     histMeanPhaseRelationIS = shankStruct.histMeanPhaseRelationIS;
%     histPhaseRelationSD = shankStruct.histPhaseRelationSD;
%     meanPhaseRelationSD = shankStruct.meanPhaseRelationSD;
%     MPRSDCu = shankStruct.MPRSDCu;
%     MPRSDCl = shankStruct.MPRSDCl;
%     histMeanPhaseRelationSD = shankStruct.histMeanPhaseRelationSD;
%     histPhaseRelationSS = shankStruct.histPhaseRelationSS;
%     meanPhaseRelationSS = shankStruct.meanPhaseRelationSS;
%     MPRSSCu = shankStruct.MPRSSCu;
%     MPRSSCl = shankStruct.MPRSSCl;
%     histMeanPhaseRelationSS = shankStruct.histMeanPhaseRelationSS;
%     histPhaseRelationSF = shankStruct.histPhaseRelationSF;
%     meanPhaseRelationSF = shankStruct.meanPhaseRelationSF;
%     MPRSFCu = shankStruct.MPRSFCu;
%     MPRSFCl = shankStruct.MPRSFCl;
%     histMeanPhaseRelationSF = shankStruct.histMeanPhaseRelationSF;
    histPhaseRelationSpk = shankStruct.histPhaseRelationSpk;
    meanPhaseRelationSpk = shankStruct.meanPhaseRelationSpk;
    MPRSpkCu = shankStruct.MPRSpkCu;
    MPRSpkCl = shankStruct.MPRSpkCl;
    PLV = shankStruct.PLV;
    MI = shankStruct.MI;
    r = shankStruct.r;
    rd = shankStruct.rd;
    histMeanPhaseRelationSpk = shankStruct.histMeanPhaseRelationSpk;
    Q = shankStruct.Q;
    for dbCountShank = 1:numel(shankStruct.db)
      if shankStruct.db(dbCountShank).entryName == fnsData{dbCount}
        animal = shankStruct.db(dbCountShank).animal;
        entryName = shankStruct.db(dbCountShank).entryName;
        break
      end
    end

% LOOP THROUGH UNITS
    nUnits = numel(units);
    for u = 1:nUnits
      fprintf('Processing %s shank %d unit %d\n', fnsData{dbCount}, sh, u);
      q = Q{u};
      
% ASSIGN UNIT-RELATED VARIABLES
%       histPhaseRelation_u = histPhaseRelation{u};
%       histPhaseRelationIS_u = histPhaseRelationIS{u};
%       histPhaseRelationSD_u = histPhaseRelationSD{u};
%       histPhaseRelationSF_u = histPhaseRelationSF{u};
%       histPhaseRelationSS_u = histPhaseRelationSS{u};
      histPhaseRelationSpk_u = histPhaseRelationSpk{u};
      
%       meanPhaseRelation_u = meanPhaseRelation(:,u);
%       meanPhaseRelationIS_u = meanPhaseRelationIS(:,u);
%       meanPhaseRelationSD_u = meanPhaseRelationSD(:,u);
%       meanPhaseRelationSF_u = meanPhaseRelationSF(:,u);
%       meanPhaseRelationSS_u = meanPhaseRelationSS(:,u);
      meanPhaseRelationSpk_u = meanPhaseRelationSpk(:,u);
      
% PLOT UNIT-RELATED DATA AND SAVE GRAPHS
%       for f = 1:numel(FOI)
%         figFileName = [entryName '_' num2str(q.unit) '_' fWidth 'FreqPhaseRelation' '_same_freq_' num2str(FOI(f))...
%           '___mean_' num2str(meanPhaseRelation_u(f))];
%         figFileName = strrep(figFileName,'.','p');
%         h = phaseHistPlot(edges, histPhaseRelation_u(f,:), '# cycles', figFileName, visibility);
%         close(h);
%       end
      
%       for f = 1:numel(meanPhaseRelationIS_u)
%         figFileName = [entryName '_' num2str(q.unit) '_' fWidth 'FreqPhaseRelation' '_P_infra-slow_' num2str(FOI(1)) '_vs_U_' num2str(FOI(f+1))...
%           '___mean_' num2str(meanPhaseRelationIS_u(f))];
%         figFileName = strrep(figFileName,'.','p');
%         h = phaseHistPlot(edges, histPhaseRelationIS_u(f,:), '# cycles', figFileName, visibility);
%         close(h);
%       end
%       
%       for f = 1:numel(meanPhaseRelationSD_u)
%         figFileName = [entryName '_' num2str(q.unit) '_' fWidth 'FreqPhaseRelation' '_P_slow_' num2str(FOI(3)) '_vs_U_delta_' num2str(FOI(5))...
%           '___mean_' num2str(meanPhaseRelationSD_u(f))];
%         figFileName = strrep(figFileName,'.','p');
%         h = phaseHistPlot(edges, histPhaseRelationSD_u(f,:), '# cycles', figFileName, visibility);
%         close(h);
%       end
%       
%       for f = 1:numel(meanPhaseRelationSS_u)
%         figFileName = [entryName '_' num2str(q.unit) '_' fWidth 'FreqPhaseRelation' '_P_slow_' num2str(FOI(3)) '_vs_U_spindle_' num2str(FOI(6))...
%           '___mean_' num2str(meanPhaseRelationSS_u(f))];
%         figFileName = strrep(figFileName,'.','p');
%         h = phaseHistPlot(edges, histPhaseRelationSS_u(f,:), '# cycles', figFileName, visibility);
%         close(h);
%       end
%       
%       for f = 1:numel(meanPhaseRelationSF_u)
%         figFileName = [entryName '_' num2str(q.unit) '_' fWidth 'FreqPhaseRelation' '_P_spindle_' num2str(FOI(6)) '_vs_U_fast_' num2str(FOI(7))...
%           '___mean_' num2str(meanPhaseRelationSF_u(f))];
%         figFileName = strrep(figFileName,'.','p');
%         h = phaseHistPlot(edges, histPhaseRelationSF_u(f,:), '# cycles', figFileName, visibility);
%         close(h);
%       end

      for f = 1:numel(FOI)
        figFileName = [entryName '_' num2str(q.unit) '_' fWidth 'FreqPhaseRelation' 'Spk' '_same_freq_' num2str(FOI(f))...
          '___mean_' num2str(meanPhaseRelationSpk_u(f)) '_PLV_' num2str(PLV{u}(f)) '_MI_' num2str(MI{u}(f))];
        figFileName = strrep(figFileName,'.','p');
        h = phaseHistPlot(edges, histPhaseRelationSpk_u(f,:), '# spikes', figFileName, visibility);
        close(h);
      end
      
      for p = 1:size(phasePairsHT,1)
        if u == 1
          figFileName = ['Correlation: Phase @ ' num2str(FOI(phasePairsHT(p,1))) ' Hz   vs   phase @ ' num2str(FOI(phasePairsHT(p,2))) ' Hz'];
          figH(p) = figProperties(figFileName, 'normalized', [0, .005, .97, .90], 'w', visibility);
          ylabel('Correlation coef. r');
          title(['Correlation: Phase @ ' num2str(FOI(phasePairsHT(p,1))) ' Hz   vs   phase @ ' num2str(FOI(phasePairsHT(p,2))) ' Hz']);
        else
          set(0, 'CurrentFigure', figH(p))
        end
        if r{u}(p,3) < significance
          hold on
          plot(r{u}(p,1), '.', 'MarkerSize',10)
          hold off
        end
        
        if u == 1
          figFileName = ['Correlation: Absolute phase @ ' num2str(FOI(phasePairsHT(p,1))) ' Hz   vs   phase @ ' num2str(FOI(phasePairsHT(p,2))) ' Hz'];
          figHD(p) = figProperties(figFileName, 'normalized', [0, .005, .97, .90], 'w', visibility);
          ylabel('Correlation coef. r');
          title(['Correlation: Phase @ ' num2str(FOI(phasePairsHT(p,1))) ' Hz   vs   phase @ ' num2str(FOI(phasePairsHT(p,2))) ' Hz']);
        else
          set(0, 'CurrentFigure', figHD(p))
        end
        if rd{u}(p,3) < significance
          hold on
          plot(rd{u}(p,1), '.', 'MarkerSize',10)
          hold off
        end
      end
    end
    
% PLOT SHANK-RELATED DATA AND SAVE GRAPHS
%     for f = 1:numel(FOI)
%       figFileName = [entryName '_' fWidth 'FreqPhaseRelation' '_same_freq_' num2str(FOI(f))];
%       figFileName = strrep(figFileName,'.','p');
%       h = phaseHistPlot(edges, histMeanPhaseRelation(f,:), '# neurons', figFileName, visibility);
%       close(h);
%     end
    
%     for f = 1:size(histMeanPhaseRelationIS,1)
%       figFileName = [entryName '_' fWidth 'FreqPhaseRelation' '_infra-slow_' num2str(FOI(1)) '_vs_U_' num2str(FOI(f+1))];
%       figFileName = strrep(figFileName,'.','p');
%       h = phaseHistPlot(edges, histMeanPhaseRelationIS(f,:), '# neurons', figFileName, visibility);
%       close(h);
%     end
%     
%     for f = 1:size(histMeanPhaseRelationSD,1)
%       figFileName = [entryName '_' fWidth 'FreqPhaseRelation' '_slow_' num2str(FOI(3)) '_vs_delta_' num2str(FOI(5))];
%       figFileName = strrep(figFileName,'.','p');
%       h = phaseHistPlot(edges, histMeanPhaseRelationSD(f,:), '# neurons', figFileName, visibility);
%       close(h);
%     end
%     
%     for f = 1:size(histMeanPhaseRelationSS,1)
%       figFileName = [entryName '_' fWidth 'FreqPhaseRelation' '_slow_' num2str(FOI(3)) '_vs_spindle_' num2str(FOI(6))];
%       figFileName = strrep(figFileName,'.','p');
%       h = phaseHistPlot(edges, histMeanPhaseRelationSS(f,:), '# neurons', figFileName, visibility);
%       close(h);
%     end
%     
%     for f = 1:size(histMeanPhaseRelationSF,1)
%       figFileName = [entryName '_' fWidth 'FreqPhaseRelation' '_spindle_' num2str(FOI(6)) '_vs_fast_' num2str(FOI(7))];
%       figFileName = strrep(figFileName,'.','p');
%       h = phaseHistPlot(edges, histMeanPhaseRelationSF(f,:), '# neurons', figFileName, visibility);
%       close(h);
%     end

    for f = 1:numel(FOI)
      figFileName = [entryName '_' fWidth 'FreqPhaseRelation' 'Spk' '_same_freq_' num2str(FOI(f))];
      figFileName = strrep(figFileName,'.','p');
      h = phaseHistPlot(edges, histMeanPhaseRelationSpk(f,:), '# neurons', figFileName, visibility);
      close(h);
    end
    
%     for p = 1:size(phasePairsHT,1)
%       figFileName = [entryName '_' fWidth 'PhaseCorr' 'Spk' '_freq_' num2str(FOI(phasePairsHT(p,1))) '_vs_' num2str(FOI(phasePairsHT(p,2)))];
%       figFileName = strrep(figFileName,'.','p');
%       rPlot(figure(figH(p)),'Correlation coef. r',figFileName);
%       
%       figFileName = [entryName '_' fWidth 'AbsPhaseCorr' 'Spk' '_freq_' num2str(FOI(phasePairsHT(p,1))) '_vs_' num2str(FOI(phasePairsHT(p,2)))];
%       figFileName = strrep(figFileName,'.','p');
%       rPlot(figure(figHD(p)),'Correlation coef. r',figFileName);
%     end

    if dbCount == numel(fnsData)
%       histMeanPhaseRelationTotal = histMeanPhaseRelation;
      histMeanPhaseRelationTotalSpk = histMeanPhaseRelationSpk;
    else
%       histMeanPhaseRelationTotal = histMeanPhaseRelationTotal + histMeanPhaseRelation;
      histMeanPhaseRelationTotalSpk = histMeanPhaseRelationTotalSpk + histMeanPhaseRelationSpk;
    end
  end
end

% PLOT TOTAL CONCATENATED DATA AND SAVE GRAPHS
% for f = 1:numel(FOI)
%   figFileName = [animal '_' fWidth 'FreqPhaseRelation' '_same_freq_' num2str(FOI(f))];
%   figFileName = strrep(figFileName,'.','p');
%   h = phaseHistPlot(edges, histMeanPhaseRelationTotal(f,:), '# units', figFileName, visibility);
%   close(h);
% end

for f = 1:numel(FOI)
  figFileName = [animal '_' fWidth 'FreqPhaseRelation' 'Spk' '_same_freq_' num2str(FOI(f))];
  figFileName = strrep(figFileName,'.','p');
  h = phaseHistPlot(edges, histMeanPhaseRelationTotalSpk(f,:), '# units', figFileName, visibility);
  close(h);
end



function h = rPlot(h,yLabel,figFileName) %#ok<DEFNU>
% A helper function of AnHT_plot for plotting phase correlations.

yLim = ylim();
yTicks = get(gca,'yTick');
ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 25, 4/3, 2, [0.005 0], 'out',...
  'off', 'k', {}, [0 2], [],...
  'on', 'k', {yLabel}, yLim, yTicks);

label = [4 2.6];
margin = [0 1];
width = 2*15-label(1)-margin(1);
height = (2*15)/(2/1.25)-label(2)-margin(2);
paperSize = resizeFig(h, ax1, width, height, label, margin, 0);
exportFig(h, [figFileName '.png'],'-dpng','-r300', paperSize);
end