function [phaseHistos, distributionStats] = phaseHistosPlotMaster(drawPhaseHistos, areas, conditions, FOI, areaPhaseFOIindividual, edges, options)
% phaseHistosPlotMaster(drawPhaseHistos, areas, conditions, FOI, areaPhaseFOIindividual, edges, options)
%
% Function controls the display of phase frequency histograms for
% individual and multiple frequencies, as well as the display of phase
% frequency maps.
% Input: drawPhaseHistos - a triple logical array controlling the display
%                          of one-dimensional phase frequency histograms,
%                          summary of 1-D histograms, and the display of
%                          2-D histograms of phase frequency maps.
%        areas - area acronyms or their comparisons.
%        conditions - recording conditions corresponding to states of
%                     vigilance.
%        FOI - frequencies of interest. They could be a single vector or a
%              cell array containing multiple vectors.
%        areaPhaseFOIindividual - a cell array of phase frequency profiles.
%        edges - phase histogram bin edges.
%        options - a structure with the following fields:
%          mainFolder
%          histosSubfolder
%          mapsSubfolder
%          figSize
%          figTitle
%          freqLim
%          phaseLimHisto
%          phaseLimMap
%          limExpansion
%          edges
%          xLabelHist
%          mask
%          iAreasOI.


%% Parse user input
if nargin < 7
  options.mainFolder = pwd;
  options.histosSubfolder = 'phaseHistos';
  options.mapsSubfolder = 'mapsHistos';
  options.figSize = 15;
  options.figTitle = 'SPIKING';
  options.freqLim = [10e-3 - 0.002   30];
  options.phaseLimHisto = [-pi pi];
  options.phaseLimMap = [-pi pi];
  options.xLabelHist = '# recordings or units';
  options.mask = [];
  options.iAreasOI = 1:numel(areas);
else
  if ~isfield(options,'mainFolder')
    options.mainFolder = pwd;
  end
  if ~isfield(options,'histosSubfolder')
    options.histosSubfolder = 'phaseHistos';
  end
  if ~isfield(options,'mapsSubfolder')
    options.histosSubfolder = 'mapsHistos';
  end
  if ~isfield(options,'figSize')
    options.figSize = 15;
  end
  if ~isfield(options,'figTitle')
    options.figTitle = 'SPIKING';
  end
  if ~isfield(options,'freqLim')
    options.freqLim = [10e-3 - 0.002   30];
  end
  if ~isfield(options,'phaseLimHisto')
    options.phaseLimHisto = [-pi pi];
  end
  if ~isfield(options,'limExpansion')
    options.limExpansion = pi/4;
  end
  if ~isfield(options,'phaseLimMap')
    options.phaseLimMap = [-pi pi];
  end
  if ~isfield(options,'xLabelHist')
    options.xLabelHist = '# recordings or units';
  end
  if ~isfield(options,'mask')
    options.mask = [];
  end
  if ~isfield(options,'iAreasOI')
    options.iAreasOI = 1:numel(areas);
  end
end

if (drawPhaseHistos(1) || drawPhaseHistos(2)) && ~exist([options.mainFolder filesep options.histosSubfolder], 'file')
  mkdir([options.mainFolder filesep options.histosSubfolder]);
end

if drawPhaseHistos(3) && ~exist([options.mainFolder filesep options.mapsSubfolder], 'file')
  mkdir([options.mainFolder filesep options.mapsSubfolder]);
end


%% Generate figures
phaseHistos = {};
distributionStats = {};
for iCond = 1:min([2 numel(conditions)])
  for iArea = 1:numel(options.iAreasOI)
    disp(['Processing data for ' conditions{iCond} ' ' areas{options.iAreasOI(iArea)}...
      ' (comparison # ' num2str((iCond-1)*numel(areas) + options.iAreasOI(iArea)) '/' num2str(numel(conditions)*numel(areas)) ')']);
    
    % Get phase and frequencies
    if iscell(FOI)
      areaFOI = FOI{iCond}{options.iAreasOI(iArea)};
    else
      areaFOI = FOI;
    end
    phase = areaPhaseFOIindividual{iCond}{options.iAreasOI(iArea)};
    if ~isempty(phase)
      phase = recentrePhase(phase, mean(options.phaseLimHisto));
      [areaFOI, sortInds] = sort(areaFOI, 'descend'); % Sort in descending order as the function was originally coded to work with descending frequencies
      reverseSortInds = 1:numel(areaFOI);
      reverseSortInds = reverseSortInds(sortInds);
      phase = phase(:,sortInds);
      histPhaseOverF = zeros(numel(edges),numel(areaFOI));
      for f = 1:numel(areaFOI)
        histPhaseOverF(:,f) = [size(phase,1) - sum(~isnan(phase(:,f))) histcounts(phase(:,f), edges)];
      end
      
      if drawPhaseHistos(1) || drawPhaseHistos(2) || drawPhaseHistos(3)
        if strcmp(options.figTitle, 'SPIKING') || strcmp(options.figTitle, 'LFP')
          compStr = areas{options.iAreasOI(iArea)};
        elseif strcmp(options.figTitle, 'PUPIL')
          compStr = [areas{options.iAreasOI(iArea)} 'VsPupil'];
        elseif strcmp(options.figTitle, 'MOTION')
          compStr = [areas{options.iAreasOI(iArea)} 'VsMotion'];
        end
        if drawPhaseHistos(2)
          plotCount = 0;
          sbPhase = figure('units', 'normalized', 'position', [0.002, .04, 1, .88], 'Visible', 'on');
        end
        for f = 1:numel(areaFOI)
          
          % Individual histos
          if drawPhaseHistos(1)
            if areaFOI(f) == 0.01 || areaFOI(f) == 0.03 || areaFOI(f) == 0.05 || areaFOI(f) == 0.1 || areaFOI(f) == 0.3 ||...
                areaFOI(f) == 0.5 || areaFOI(f) == 1 || areaFOI(f) == 4 || areaFOI(f) == 10 || f == numel(FOI)
              phaseHistPlotUnadjusted(edges, histPhaseOverF(:,f), options.xLabelHist, [options.mainFolder filesep options.histosSubfolder filesep...
                options.figTitle '_' compStr '_' conditions{iCond} '_phase_at_' num2str(areaFOI(f)) '_Hz'], 'on', mean(options.phaseLimHisto));
              [distributionStats.pRayleigh{iCond}{options.iAreasOI(iArea)}{f}, distributionStats.zRayleigh{iCond}{options.iAreasOI(iArea)}{f}, distributionStats.pOmnibus{iCond}{options.iAreasOI(iArea)}{f},...
                distributionStats.mOmnibus{iCond}{options.iAreasOI(iArea)}{f}, distributionStats.pRao{iCond}{options.iAreasOI(iArea)}{f}, distributionStats.U_Rao{iCond}{options.iAreasOI(iArea)}{f},...
                distributionStats.pV{iCond}{options.iAreasOI(iArea)}{f}, distributionStats.vV{iCond}{options.iAreasOI(iArea)}{f}, distributionStats.pHR{iCond}{options.iAreasOI(iArea)}{f},...
                distributionStats.T_HR{iCond}{options.iAreasOI(iArea)}{f}, distributionStats.nModes{iCond}{options.iAreasOI(iArea)}{f}, distributionStats.excessMass{iCond}{options.iAreasOI(iArea)}{f},...
                distributionStats.U2_KDE{iCond}{options.iAreasOI(iArea)}{f}, distributionStats.pKDE{iCond}{options.iAreasOI(iArea)}{f}, distributionStats.modes{iCond}{options.iAreasOI(iArea)}{f},...
                distributionStats.dipHDT{iCond}{options.iAreasOI(iArea)}{f}, distributionStats.pHDT{iCond}{options.iAreasOI(iArea)}{f}] = phaseDistributionTests(...
                phase(:,f), edges, true); % Stats
            else
              distributionStats.pRayleigh{iCond}{options.iAreasOI(iArea)}{f} = [];
              distributionStats.zRayleigh{iCond}{options.iAreasOI(iArea)}{f} = [];
              distributionStats.pOmnibus{iCond}{options.iAreasOI(iArea)}{f} = [];
              distributionStats.mOmnibus{iCond}{options.iAreasOI(iArea)}{f} = [];
              distributionStats.pRao{iCond}{options.iAreasOI(iArea)}{f} = [];
              distributionStats.U_Rao{iCond}{options.iAreasOI(iArea)}{f} = [];
              distributionStats.pV{iCond}{options.iAreasOI(iArea)}{f} = [];
              distributionStats.vV{iCond}{options.iAreasOI(iArea)}{f} = [];
              distributionStats.pHR{iCond}{options.iAreasOI(iArea)}{f} = [];
              distributionStats.T_HR{iCond}{options.iAreasOI(iArea)}{f} = [];
              distributionStats.nModes{iCond}{options.iAreasOI(iArea)}{f} = [];
              distributionStats.excessMass{iCond}{options.iAreasOI(iArea)}{f} = [];
              distributionStats.U2_KDE{iCond}{options.iAreasOI(iArea)}{f} = [];
              distributionStats.pKDE{iCond}{options.iAreasOI(iArea)}{f} = [];
              distributionStats.modes{iCond}{options.iAreasOI(iArea)}{f} = [];
              distributionStats.dipHDT{iCond}{options.iAreasOI(iArea)}{f} = [];
              distributionStats.pHDT{iCond}{options.iAreasOI(iArea)}{f} = [];
            end
          end
          
          % Summary histos
          if drawPhaseHistos(2)
            if areaFOI(f) == 0.01 || areaFOI(f) == 0.03 || areaFOI(f) == 0.05 || areaFOI(f) == 0.1 || areaFOI(f) == 0.3 ||...
                areaFOI(f) == 0.5 || areaFOI(f) == 1 || areaFOI(f) == 4 || areaFOI(f) == 10
              figure(sbPhase);
              plotCount = plotCount + 1;
              subplot(3, 3, 10 - plotCount);
              if plotCount <= 2
                phaseHistSubplotUnadjusted(edges, histPhaseOverF(:,f), '', 'Phase (rad)', 'k', mean(options.phaseLimHisto));
                dispFOI(areaFOI(f), sum(histPhaseOverF(:,f)));
              elseif plotCount == 3
                phaseHistSubplotUnadjusted(edges, histPhaseOverF(:,f), options.xLabelHist, 'Phase (rad)', 'k', mean(options.phaseLimHisto));
                dispFOI(areaFOI(f), sum(histPhaseOverF(:,f)));
              elseif plotCount <= 5
                phaseHistSubplotUnadjusted(edges, histPhaseOverF(:,f), '', '', 'k', mean(options.phaseLimHisto));
                dispFOI(areaFOI(f), sum(histPhaseOverF(:,f)));
              elseif plotCount == 6
                phaseHistSubplotUnadjusted(edges, histPhaseOverF(:,f), options.xLabelHist, '', 'k', mean(options.phaseLimHisto));
                dispFOI(areaFOI(f), sum(histPhaseOverF(:,f)));
              elseif plotCount <= 8
                phaseHistSubplotUnadjusted(edges, histPhaseOverF(:,f), '', '', 'k', mean(options.phaseLimHisto));
                dispFOI(areaFOI(f), sum(histPhaseOverF(:,f)));
              elseif plotCount == 9
                phaseHistSubplotUnadjusted(edges, histPhaseOverF(:,f), options.xLabelHist, '', 'k', mean(options.phaseLimHisto));
                dispFOI(areaFOI(f), sum(histPhaseOverF(:,f)));
                set(gcf, 'Color','w');
                print(sbPhase, [options.mainFolder filesep options.histosSubfolder filesep options.figTitle '_' compStr '_' conditions{iCond}...
                  '_phase_at_p01_p03_p05_p1_p3_p5_1_4_10_Hz.png'], '-dpng', '-r300');
                hgsave([options.mainFolder filesep options.histosSubfolder filesep options.figTitle '_' compStr '_' conditions{iCond}...
                  '_phase_at_p01_p03_p05_p1_p3_p5_1_4_10_Hz.fig'])
                close(sbPhase);
              end
            end
          end
        end
        
        % Phase frequency maps
        if drawPhaseHistos(3)
          cm = 'cold';
          
          % Regular
          edgesExtended = [edges(1:end-1)-2*pi edges edges(2:end)+2*pi];
          histPhaseOverFExtended = fliplr([histPhaseOverF(2:end,:); histPhaseOverF(2:end,:); histPhaseOverF(2:end,:)]);
          h = phaseGraph([options.figTitle ' Phase map ' compStr ' ' conditions{iCond}], fliplr(areaFOI), options.freqLim,...
            [0.001 0.01 0.1 1 10 30 100 1000], edgesExtended,...
            [options.phaseLimMap(1)-options.limExpansion options.phaseLimMap(end)+options.limExpansion], [-2*pi -3*pi/2 -pi -pi/2 0 pi/2 pi 3*pi/2 2*pi],...
            histPhaseOverFExtended, cm, 'on', mean(options.phaseLimMap));
          if ~isempty(options.mask)
            hold on
            for iMask = 1:numel(options.mask)
              assert(options.mask{iMask}(1,1) == options.mask{iMask}(3,1) && options.mask{iMask}(2,1) == options.mask{iMask}(4,1));
              ciplot([options.mask{iMask}(3,2) options.mask{iMask}(4,2)], [options.mask{iMask}(1,2) options.mask{iMask}(2,2)],...
                [options.mask{iMask}(1,1) options.mask{iMask}(2,1)], 'w');
            end
            hold off
          end
          ax1.XTickLabel = {'0.001','0.01','0.1','1','10','30','100','1000'};
          set(h,'color','white')
          hgsave(h, [options.mainFolder filesep options.mapsSubfolder filesep options.figTitle '_Phase_map_' compStr '_' conditions{iCond}]);
          print(h, [options.mainFolder filesep options.mapsSubfolder filesep options.figTitle '_Phase_map_' compStr '_' conditions{iCond} '.png'],'-dpng','-r300');
          close(h);
          
          % Scaled
          h = phaseGraph([options.figTitle ' Scaled phase map ' compStr ' ' conditions{iCond}], fliplr(areaFOI), options.freqLim,...
            [0.001 0.01 0.1 1 10 30 100 1000], edgesExtended,...
            [options.phaseLimMap(1)-options.limExpansion options.phaseLimMap(end)+options.limExpansion], [-2*pi -3*pi/2 -pi -pi/2 0 pi/2 pi 3*pi/2 2*pi],...
            histPhaseOverFExtended, cm, 'on', mean(options.phaseLimMap), true);
          if ~isempty(options.mask)
            hold on
            for iMask = 1:numel(options.mask)
              assert(options.mask{iMask}(1,1) == options.mask{iMask}(3,1) && options.mask{iMask}(2,1) == options.mask{iMask}(4,1));
              ciplot([options.mask{iMask}(3,2) options.mask{iMask}(4,2)], [options.mask{iMask}(1,2) options.mask{iMask}(2,2)],...
                [options.mask{iMask}(1,1) options.mask{iMask}(2,1)], 'w');
            end
            hold off
          end
          ax1.XTickLabel = {'0.001','0.01','0.1','1','10','30','100','1000'};
          set(h,'color','white')
          hgsave(h, [options.mainFolder filesep options.mapsSubfolder filesep options.figTitle '_Scaled_phase_map_' compStr '_' conditions{iCond}]);
          print(h, [options.mainFolder filesep options.mapsSubfolder filesep options.figTitle '_Scaled_phase_map_' compStr '_' conditions{iCond} '.png'],'-dpng','-r300');
          close(h);
        end
        
        % Asign output
        phaseHistos{iCond}{options.iAreasOI(iArea)} = histPhaseOverF(:,reverseSortInds); %#ok<*AGROW>
        if drawPhaseHistos(1)
          distributionStats.pRayleigh{iCond}{options.iAreasOI(iArea)} = distributionStats.pRayleigh{iCond}{options.iAreasOI(iArea)}(reverseSortInds);
          distributionStats.zRayleigh{iCond}{options.iAreasOI(iArea)} = distributionStats.zRayleigh{iCond}{options.iAreasOI(iArea)}(reverseSortInds);
          distributionStats.pOmnibus{iCond}{options.iAreasOI(iArea)} = distributionStats.pOmnibus{iCond}{options.iAreasOI(iArea)}(reverseSortInds);
          distributionStats.mOmnibus{iCond}{options.iAreasOI(iArea)} = distributionStats.mOmnibus{iCond}{options.iAreasOI(iArea)}(reverseSortInds);
          distributionStats.pRao{iCond}{options.iAreasOI(iArea)} = distributionStats.pRao{iCond}{options.iAreasOI(iArea)}(reverseSortInds);
          distributionStats.U_Rao{iCond}{options.iAreasOI(iArea)} = distributionStats.U_Rao{iCond}{options.iAreasOI(iArea)}(reverseSortInds);
          distributionStats.pV{iCond}{options.iAreasOI(iArea)} = distributionStats.pV{iCond}{options.iAreasOI(iArea)}(reverseSortInds);
          distributionStats.vV{iCond}{options.iAreasOI(iArea)} = distributionStats.vV{iCond}{options.iAreasOI(iArea)}(reverseSortInds);
          distributionStats.pHR{iCond}{options.iAreasOI(iArea)} = distributionStats.pHR{iCond}{options.iAreasOI(iArea)}(reverseSortInds);
          distributionStats.T_HR{iCond}{options.iAreasOI(iArea)} = distributionStats.T_HR{iCond}{options.iAreasOI(iArea)}(reverseSortInds);
          distributionStats.nModes{iCond}{options.iAreasOI(iArea)} = distributionStats.nModes{iCond}{options.iAreasOI(iArea)}(reverseSortInds);
          distributionStats.excessMass{iCond}{options.iAreasOI(iArea)} = distributionStats.excessMass{iCond}{options.iAreasOI(iArea)}(reverseSortInds);
          distributionStats.U2_KDE{iCond}{options.iAreasOI(iArea)} = distributionStats.U2_KDE{iCond}{options.iAreasOI(iArea)}(reverseSortInds);
          distributionStats.pKDE{iCond}{options.iAreasOI(iArea)} = distributionStats.pKDE{iCond}{options.iAreasOI(iArea)}(reverseSortInds);
          distributionStats.modes{iCond}{options.iAreasOI(iArea)} = distributionStats.modes{iCond}{options.iAreasOI(iArea)}(reverseSortInds);
          distributionStats.dipHDT{iCond}{options.iAreasOI(iArea)} = distributionStats.dipHDT{iCond}{options.iAreasOI(iArea)}(reverseSortInds);
          distributionStats.pHDT{iCond}{options.iAreasOI(iArea)} = distributionStats.pHDT{iCond}{options.iAreasOI(iArea)}(reverseSortInds);
        end
      end
    end
  end
end