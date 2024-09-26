function [slicePhaseFig, rFOI, pvalFOI, nFOI, stats] = phaseLocationProfile(phase, edges, locations, FOI, nSlices, area, condition, displayGraphs, figPath, modeTest, prefix, suffix)
% [slicePhaseFig, rFOI, pvalFOI, nFOI] = phaseLocationProfile(phase, edges, locations, FOI, nSlices, area, condition, displayGraphs, figPath, modTest)
%
% Function performs correlations between unit probe location and phase in
% addition to producing subplot graphs with location phase profile
% histograms corresponding to location profile slices.
% Input: phase - phase vector for a frequency of interest. It can also be a
%                matrix. In that case columns of the matrix correspond to
%                different frequencies of interest.
%        edges for phase histogram binning (see histcounts).
%        locations - probe location vector corresponding to rows in the
%                    phase matrix.
%        FOI - a frequency or frequencies of interest.
%        nSlices - number of location profile slices. It can be a vector.
%                  In that case figures are produced for all different
%                  slicing densities.
%        area - area comparison string.
%        condition - condition string.
%        displayGraphs - if true, displays graphs. Otherwise only
%                        correlation analyses are carried out without
%                        graphs being produced. Default is true.
%        figPath - a path to a folder for saving figures. Default is the
%                  present working directory.
%        modeTest - if true, carries out data modality (number of modes)
%                   estimation and determines mode properties. Default is
%                   false.
%        prefix and suffix of figure file names.
% Output: slicePhaseFig - figure handles. If no figures were displayed, it
%                         would be empty.
%         rFOI - linear-circular correlation coefficients for
%                location-phase correlations for every frequency of
%                interest.
%         pvalFOI - corresponding p-values.
%         nFOI - count of units with significant phase values/ total number
%                of units for every frequency of interest.
%         stats - results of statistical distribution uniformity and
%                 modality tests. stats has the following fields:
%                 slices - profile slicing density;
%                 anatomicalRange - slice anatomical range in micrometers.
%                                   It is a 3D field with the first
%                                   dimension being slicing density, the
%                                   second being frequency, and the third
%                                   being slice number. The rest of the
%                                   fields of stats variable have the same
%                                   structure.
%                 pRayleigh
%                 zRayleigh
%                 pOmnibus
%                 mOmnibus
%                 pRao
%                 U_Rao
%                 pV
%                 vV
%                 pHR
%                 T_HR
%                 nModes
%                 excessMass
%                 U2_KDE
%                 pKDE
%                 modes
%                 pHDT
%                 dipHDT

if nargin < 12
  suffix = '';
end

if nargin < 11
  prefix = 'SPIKING';
end

if nargin < 10
  modeTest = false;
end

if nargin < 9 || isempty(figPath)
  figPath = pwd;
end

if nargin < 8
  displayGraphs = true;
end
 
phase = recentrePhase(phase, 0);

% Correlation analyses
[rFOI, pvalFOI] = corrMulti(phase', repmat(locations, 1, numel(FOI))', 'circlinearnp');
for f = 1:numel(FOI)
  nFOI{f} = [sum(~isnan(phase(:,f))) numel(phase(:,f))]; %#ok<*AGROW>
end

% Draw graphs
xLabel = 'Location (\mum)';
yLabel = 'Unit phase (rad)';
yTickLabel = {'-\pi','-\pi/2','0','\pi/2','\pi'};
if displayGraphs
  minLoc = min(locations);
  maxLoc = max(locations);
  for sliceSize = 1:numel(nSlices) % loop over increasing slicing density
    if sliceSize == 1
      stats.slices = nSlices;
    end
    for f = 1:numel(FOI) % loop over frequencies of interest
      % The phase location profile
      if nSlices(sliceSize) < 4
        subplotRows = 2;
        slicePhaseFig(f) = subplotFig(subplotRows);
      else
        subplotRows = 3;
        slicePhaseFig(f) = subplotFig(4);
      end
      subplot(subplotRows, 3, 1);
      [locations, sortOrder] = sort(locations);
      phaseFOI = phase(sortOrder,f);
      xySubplot(locations, phaseFOI, xLabel, [], [], [minLoc maxLoc], yLabel, [], yTickLabel, [], rFOI(f), pvalFOI(f), nFOI{f},...
        false, true, 'linear-circular')
      topAxes = gca;
      
      for slice = 1:nSlices(sliceSize) % loop over slices of the same size
        
        % The location range histograms
        locSlice = (maxLoc-minLoc)/nSlices(sliceSize);
        locSliceRange = [minLoc + locSlice*(slice-1) minLoc + locSlice*slice];
        if slice == nSlices(sliceSize)
          locSliceRange(2) = locSliceRange(2) + 1;
        end
        axes(topAxes); %#ok<*LAXES>
        hold on
        yLim = ylim;
        rangeColour = matlabColours(slice);
        ciplot([yLim(1) yLim(1)], [yLim(2) yLim(2)], locSliceRange, rangeColour, 0.2); % mark the location range in the main graph
        if slice == nSlices(sliceSize)
          dispFOI(FOI(f)); % display frequency of interest
        end
        hold off
        locSliceInds = logical(locations >= locSliceRange(1) & locations < locSliceRange(2));
        phaseSlice = phaseFOI(locSliceInds);
        histSlice = [sum(sum(isnan(phaseSlice))) histcounts(phaseSlice, edges)];
        subplot(subplotRows, 3, 3 + slice);
        if slice == 1 || slice == 4
          if slice + 3 > nSlices(sliceSize)
            phaseHistSubplotUnadjusted(edges, histSlice, '# units', yLabel, rangeColour); % draw the histograms
          else
            phaseHistSubplotUnadjusted(edges, histSlice, '# units', '', rangeColour);
          end
        else
          if slice + 3 > nSlices(sliceSize)
            phaseHistSubplotUnadjusted(edges, histSlice, '', yLabel, rangeColour); % draw the histograms
          else
            phaseHistSubplotUnadjusted(edges, histSlice, '', '', rangeColour);
          end
        end
        
        % Statistical tests for phase distribution uniformity and modality
        [pRayleigh, zRayleigh, pOmnibus, mOmnibus, pRao, U_Rao, pV, vV, pHR, T_HR,...
          nModes, excessMass, U2_KDE, pKDE, modes, dipHDT, pHDT] = phaseDistributionTests(phaseSlice, edges, modeTest);
        
        % Print out the statistical testing results
        if ~isnan(nModes)
          axes(topAxes)
          testMsg = sprintf('%.6g-%.6g\\mum: pRayleigh=%g(%.6g), pHermansRasson=%g(%.6g), pHartigan=%g(%.6g)',...
            round(locSliceRange(1),2), round(locSliceRange(2),2), round(pRayleigh,4), round(zRayleigh,2), round(pHR,4), round(T_HR,2),...
            round(pHDT,4), round(dipHDT,2));
          xLim = xlim;
          xAxisLength = xLim(2)-xLim(1);
          yLim = ylim;
          yAxisLength = yLim(2)-yLim(1);
          text(xLim(2)+xAxisLength*0.1, yLim(2)-yAxisLength*0.05-yAxisLength*0.2*(slice-1), testMsg, 'FontSize',16);
          if modeTest
            modeMsg = [];
            for iFit = 1:numel(nModes)
              if numel(nModes) == 1 || (numel(nModes) > 1 && pKDE(iFit) < 0.05)
                modeMsg = sprintf('#modes=%g, S_k=%.6g, pU2=%g(%.6g), modes:', nModes(iFit), round(excessMass(iFit),4), round(pKDE(iFit),4),...
                  round(U2_KDE(iFit),2));
                for iMode = 1:nModes(iFit)
                  modeMsg = [modeMsg ' ' num2str(round(modes{iFit}.peaks(iMode),2)) '(' num2str(round(modes{iFit}.centreMass(iMode),2)) ')'];
                end
                break
              end
            end
            if ~isempty(modeMsg)
              text(xLim(2)+xAxisLength*0.15, yLim(2)-yAxisLength*0.15-yAxisLength*0.2*(slice-1), modeMsg, 'FontSize',16);
            end
          end
        end
        
        % Compile stats output variable
        stats.anatomicalRange{sliceSize}{f}{slice} = locSliceRange;
        stats.pRayleigh{sliceSize}{f}{slice} = pRayleigh;
        stats.zRayleigh{sliceSize}{f}{slice} = zRayleigh;
        stats.pOmnibus{sliceSize}{f}{slice} = pOmnibus;
        stats.mOmnibus{sliceSize}{f}{slice} = mOmnibus;
        stats.pRao{sliceSize}{f}{slice} = pRao;
        stats.U_Rao{sliceSize}{f}{slice} = U_Rao;
        stats.pV{sliceSize}{f}{slice} = pV;
        stats.vV{sliceSize}{f}{slice} = vV;
        stats.pHR{sliceSize}{f}{slice} = pHR;
        stats.T_HR{sliceSize}{f}{slice} = T_HR;
        stats.nModes{sliceSize}{f}{slice} = nModes;
        stats.excessMass{sliceSize}{f}{slice} = excessMass;
        stats.U2_KDE{sliceSize}{f}{slice} = U2_KDE;
        stats.pKDE{sliceSize}{f}{slice} = pKDE;
        stats.modes{sliceSize}{f}{slice} = modes;
        stats.pHDT{sliceSize}{f}{slice} = pHDT;
        stats.dipHDT{sliceSize}{f}{slice} = dipHDT;
      end
      
      % Save the figures
      inds = ~isnan(locations) & ~isnan(phaseFOI);
      [~, slope] = fitLine(locations(inds), phaseFOI(inds), 'linear-circular');
      if slope < 0 && rFOI(f) > 0
        figTitle = ['Unit location and phase correlations and histograms for ' area ' ' condition ' ' num2str(FOI(f))...
          ' Hz ' num2str(nSlices(sliceSize)) ' slices: r=' num2str(-rFOI(f)) ' p=' num2str(pvalFOI(f))...
          ' n=' num2str(nFOI{f}(1)) '/' num2str(nFOI{f}(2))];
      else
        figTitle = ['Unit location and phase correlations and histograms for ' area ' ' condition ' ' num2str(FOI(f))...
          ' Hz ' num2str(nSlices(sliceSize)) ' slices: r=' num2str(rFOI(f)) ' p=' num2str(pvalFOI(f))...
          ' n=' num2str(nFOI{f}(1)) '/' num2str(nFOI{f}(2))];
      end
      set(slicePhaseFig(f), 'Name',figTitle);
      figName = [prefix '_unit_location_and_phase_correlations_and_histograms_for_' area suffix '_' condition '_'...
        num2str(FOI(f)) '_Hz_' num2str(nSlices(sliceSize)) '_slices'];
      set (slicePhaseFig(f), 'PaperUnits','centimeters' , 'PaperSize',[25 10])
      print(slicePhaseFig(f), [figPath filesep figName '.png'], '-dpng', '-r300');
      hgsave(slicePhaseFig(f), [figPath filesep figName '.fig']);
    end
  end
else
  slicePhaseFig = [];
end