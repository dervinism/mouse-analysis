function sumPhases(areas, conditions, FOI, areaPhaseFOIindividual, options)
% sumPhases(areas, conditions, FOI, areaPhaseFOIindividual, options)
%
% Function adds cross-area phase frequency profiles of intermediate areas
% in the information processing hierarchy. These sums should perdict
% cross-area phase frequency profiles of areas lying at the extremes of the
% information processing hierarchy.
% Input: areas - area acronyms or their comparisons.
%        conditions - recording conditions corresponding to states of
%                     vigilance.
%        FOI - frequencies of interest. They could be a single vector or a
%              cell array containing multiple vectors.
%        areaPhaseFOIindividual - a cell array of phase frequency profiles.
%        options - a structure with the following fields:
%          mainFolder
%          phaseSumPredictionsSubfolder
%          figSize
%          freqLim
%          phaseLim.

meanPhaseProfiles = {};
for iCond = 1:numel(conditions)
  for iArea = 1:numel(areas)
    if ~isempty(areaPhaseFOIindividual{iCond}{iArea})
      
      % Obtain mean phase values and 95% confidence intervals for area comparisons of interest
      if ~(strcmpi(areas{iArea}, 'VB1VsS1') || strcmpi(areas{iArea}, 'S1VsCA1') || strcmpi(areas{iArea}, 'VB1VsCA1') ||...
          strcmpi(areas{iArea}, 'VB2VsS1') || strcmpi(areas{iArea}, 'S1VsCA1') || strcmpi(areas{iArea}, 'VB2VsCA1') ||...
          strcmpi(areas{iArea}, 'VB1VsRSC') || strcmpi(areas{iArea}, 'RSCVsCA1') ||...
          strcmpi(areas{iArea}, 'VB2VsRSC') || strcmpi(areas{iArea}, 'RSCVsCA1'))
        continue
      end
      if iscell(FOI)
        areaFOI = FOI{iCond}{iArea};
      else
        areaFOI = FOI;
      end
      [phaseMean, phaseCI95] = datamean(areaPhaseFOIindividual, 'circularNP');
      phaseCI95(phaseCI95 == 0) = NaN;
      
      % Store the mean phase and confidence interval data
      if strcmpi(areas{iArea}, 'VB1VsS1')
        meanPhaseProfiles{iCond}.VB1S1.phaseComparison{1} = areas{iArea}; %#ok<*AGROW>
        meanPhaseProfiles{iCond}.VB1S1.phaseMean{1} = phaseMean;
        meanPhaseProfiles{iCond}.VB1S1.phaseStd{1} = phaseStd;
        meanPhaseProfiles{iCond}.VB1S1.phaseCI95{1} = phaseCI95;
      elseif strcmpi(areas{iArea}, 'S1VsCA1')
        meanPhaseProfiles{iCond}.VB1S1.phaseComparison{2} = areas{iArea};
        meanPhaseProfiles{iCond}.VB1S1.phaseMean{2} = phaseMean;
        meanPhaseProfiles{iCond}.VB1S1.phaseStd{2} = phaseStd;
        meanPhaseProfiles{iCond}.VB1S1.phaseCI95{2} = phaseCI95;
        meanPhaseProfiles{iCond}.VB2S1.phaseComparison{2} = areas{iArea};
        meanPhaseProfiles{iCond}.VB2S1.phaseMean{2} = phaseMean;
        meanPhaseProfiles{iCond}.VB2S1.phaseStd{2} = phaseStd;
        meanPhaseProfiles{iCond}.VB2S1.phaseCI95{2} = phaseCI95;
      elseif strcmpi(areas{iArea}, 'VB1VsCA1')
        meanPhaseProfiles{iCond}.VB1S1.phaseComparison{3} = areas{iArea};
        meanPhaseProfiles{iCond}.VB1S1.phaseMean{3} = phaseMean;
        meanPhaseProfiles{iCond}.VB1S1.phaseStd{3} = phaseStd;
        meanPhaseProfiles{iCond}.VB1S1.phaseCI95{3} = phaseCI95;
        meanPhaseProfiles{iCond}.VB1RSC.phaseComparison{3} = areas{iArea};
        meanPhaseProfiles{iCond}.VB1RSC.phaseMean{3} = phaseMean;
        meanPhaseProfiles{iCond}.VB1RSC.phaseStd{3} = phaseStd;
        meanPhaseProfiles{iCond}.VB1RSC.phaseCI95{3} = phaseCI95;
      elseif strcmpi(areas{iArea}, 'VB2VsS1')
        meanPhaseProfiles{iCond}.VB2S1.phaseComparison{1} = areas{iArea};
        meanPhaseProfiles{iCond}.VB2S1.phaseMean{1} = phaseMean;
        meanPhaseProfiles{iCond}.VB2S1.phaseStd{1} = phaseStd;
        meanPhaseProfiles{iCond}.VB2S1.phaseCI95{1} = phaseCI95;
      elseif strcmpi(areas{iArea}, 'VB2VsCA1')
        meanPhaseProfiles{iCond}.VB2S1.phaseComparison{3} = areas{iArea};
        meanPhaseProfiles{iCond}.VB2S1.phaseMean{3} = phaseMean;
        meanPhaseProfiles{iCond}.VB2S1.phaseStd{3} = phaseStd;
        meanPhaseProfiles{iCond}.VB2S1.phaseCI95{3} = phaseCI95;
        meanPhaseProfiles{iCond}.VB2RSC.phaseComparison{3} = areas{iArea};
        meanPhaseProfiles{iCond}.VB2RSC.phaseMean{3} = phaseMean;
        meanPhaseProfiles{iCond}.VB2RSC.phaseStd{3} = phaseStd;
        meanPhaseProfiles{iCond}.VB2RSC.phaseCI95{3} = phaseCI95;
      elseif strcmpi(areas{iArea}, 'VB1VsRSC')
        meanPhaseProfiles{iCond}.VB1RSC.phaseComparison{1} = areas{iArea};
        meanPhaseProfiles{iCond}.VB1RSC.phaseMean{1} = phaseMean;
        meanPhaseProfiles{iCond}.VB1RSC.phaseStd{1} = phaseStd;
        meanPhaseProfiles{iCond}.VB1RSC.phaseCI95{1} = phaseCI95;
      elseif strcmpi(areas{iArea}, 'RSCVsCA1')
        meanPhaseProfiles{iCond}.VB1RSC.phaseComparison{2} = areas{iArea};
        meanPhaseProfiles{iCond}.VB1RSC.phaseMean{2} = phaseMean;
        meanPhaseProfiles{iCond}.VB1RSC.phaseStd{2} = phaseStd;
        meanPhaseProfiles{iCond}.VB1RSC.phaseCI95{2} = phaseCI95;
        meanPhaseProfiles{iCond}.VB2RSC.phaseComparison{2} = areas{iArea};
        meanPhaseProfiles{iCond}.VB2RSC.phaseMean{2} = phaseMean;
        meanPhaseProfiles{iCond}.VB2RSC.phaseStd{2} = phaseStd;
        meanPhaseProfiles{iCond}.VB2RSC.phaseCI95{2} = phaseCI95;
      elseif strcmpi(areas{iArea}, 'VB2VsRSC')
        meanPhaseProfiles{iCond}.VB2RSC.phaseComparison{1} = areas{iArea};
        meanPhaseProfiles{iCond}.VB2RSC.phaseMean{1} = phaseMean;
        meanPhaseProfiles{iCond}.VB2RSC.phaseStd{1} = phaseStd;
        meanPhaseProfiles{iCond}.VB2RSC.phaseCI95{1} = phaseCI95;
      end
    end
  end
end

% Draw the figures
for iCond = 1:numel(meanPhaseProfiles)
  fnsData = fieldnames(meanPhaseProfiles{iCond});
  for iProfile = 1:numel(fnsData)
    delete pCross
    plotLegendsAdd{iProfile} = {};
    txtLegendsAdd{iProfile} = {};
    phaseProfile = meanPhaseProfiles{iCond}.(fnsData{iProfile});
    fH{iCond}{iArea} = figure;
    for subProfile = 1:numel(phaseProfile.phaseMean)+1
      if subProfile == 1
        colourStr = 'g';
      elseif subProfile == 2
        colourStr = 'b';
      elseif subProfile == 3
        colourStr = [13 152 186]./256;
      elseif subProfile == 4
        colourStr = 'r';
      end
      if subProfile < 4
        if subProfile == 1
          semilogx([10e-4 10e2],[0 0], 'k:'); hold on
        end
        inds = ~isnan(phaseProfile.phaseMean{subProfile}) & ~isnan(phaseProfile.phaseCI95{subProfile});
        if ~isempty(phaseProfile.phaseMean{subProfile}(inds))
          indsInterp = ~inds;
          phaseMeanInterp = interp1(areaFOI(inds), phaseProfile.phaseMean{subProfile}(inds), areaFOI, 'linear', 'extrap');
          meanPhaseProfiles{iCond}.(fnsData{iProfile}).phaseMeanInterp{subProfile} = phaseMeanInterp;
          phaseCI95Interp = interp1(areaFOI(inds), phaseProfile.phaseCI95{subProfile}(inds), areaFOI, 'linear', 'extrap');
          meanPhaseProfiles{iCond}.(fnsData{iProfile}).phaseCI95Interp{subProfile} = phaseCI95Interp;
          
          % Draw interpolated and extrapolated means
          p = semilogx(areaFOI, phaseMeanInterp,...
            '-', 'color',colourStr, 'LineWidth',2, 'MarkerSize',5, 'MarkerEdgeColor',colourStr, 'MarkerFaceColor',colourStr);
          uistack(p,'bottom');
          if isempty(plotLegendsAdd)
            plotLegendsAdd{iProfile} = p;
          else
            plotLegendsAdd{iProfile} = [plotLegendsAdd{iProfile} p];
          end
          txtLegendsAdd{iProfile}{numel(txtLegendsAdd{iProfile})+1} = ['Mean \Phi '...
            meanPhaseProfiles{iCond}.(fnsData{iProfile}).phaseComparison{subProfile}];
          
          % Mark interpolations and extrapolations
          if sum(indsInterp)
            pCross = semilogx(areaFOI(indsInterp), phaseMeanInterp(indsInterp),...
              'x', 'color',[0 0 0], 'MarkerSize',15, 'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[0 0 0]);
            semilogx(areaFOI(indsInterp), phaseMeanInterp(indsInterp),...
              'x', 'color',colourStr, 'MarkerSize',15, 'MarkerEdgeColor',colourStr, 'MarkerFaceColor',colourStr);
          end
          
          % Draw 95% confidence intervals
          pC1 = semilogx(areaFOI, phaseMeanInterp + phaseCI95Interp,...
            ':', 'color',colourStr, 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor',colourStr, 'MarkerFaceColor',colourStr);
          uistack(pC1,'bottom');
          
          pC2 = semilogx(areaFOI, phaseMeanInterp - phaseCI95Interp,...
            ':', 'color',colourStr, 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor',colourStr, 'MarkerFaceColor',colourStr);
          uistack(pC2,'bottom');
        end
        
        % Add interpolation sign to the legend
        if subProfile == 3 && exist('pCross', 'var')
          plotLegendsAdd{iProfile} = [plotLegendsAdd{iProfile} pCross];
          txtLegendsAdd{iProfile}{numel(txtLegendsAdd{iProfile})+1} = 'Mean \Phi interp.';
        end
      else
        if isfield(meanPhaseProfiles{iCond}.(fnsData{iProfile}), 'phaseMeanInterp')
          meanPhaseProfilesInterpSum = ...
            meanPhaseProfiles{iCond}.(fnsData{iProfile}).phaseMeanInterp{1} + meanPhaseProfiles{iCond}.(fnsData{iProfile}).phaseMeanInterp{2};
          p = semilogx(areaFOI, meanPhaseProfilesInterpSum,...
            '-', 'color',colourStr, 'LineWidth',2, 'MarkerSize',5, 'MarkerEdgeColor',colourStr, 'MarkerFaceColor',colourStr);
          uistack(p,'bottom');
          if isempty(plotLegendsAdd)
            plotLegendsAdd{iProfile} = p;
          else
            plotLegendsAdd{iProfile} = [plotLegendsAdd{iProfile} p];
          end
          txtLegendsAdd{iProfile}{numel(txtLegendsAdd{iProfile})+1} = ['Mean \Phi '...
            meanPhaseProfiles{iCond}.(fnsData{iProfile}).phaseComparison{1} ' + ' meanPhaseProfiles{iCond}.(fnsData{iProfile}).phaseComparison{2}];
        end
        hold off
        
        legend(plotLegendsAdd{iProfile},txtLegendsAdd{iProfile}, 'Interpreter','tex');
        legend boxoff
        title(['Sum of phase profiles: ' conditions{iCond} ' '...
          meanPhaseProfiles{iCond}.(fnsData{iProfile}).phaseComparison{1} ' '...
          meanPhaseProfiles{iCond}.(fnsData{iProfile}).phaseComparison{2} ' '...
          meanPhaseProfiles{iCond}.(fnsData{iProfile}).phaseComparison{3}])
        set(fH{iCond}{iArea},'color','w');
        ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 12, 4/3, 2, [0.005 0], 'out',...
          'on', 'k', {'Frequency (Hz)'}, options.freqLim, [0.001 0.01 0.1 1 10 30 100 1000],...
          'on', 'k', {'Phase (rad)'}, options.phaseLim, [-2*pi -3*pi/2 -pi -pi/2 0 pi/2 pi 3*pi/2 2*pi]);
        ax1.YTickLabel = {'-2\pi','-3\pi/2','-\pi','-\pi/2','0','\pi/2','\pi','3\pi/2','2\pi'};
        
        set(fH{iCond}{iArea}, 'Name',['Sum of phase profiles: ' conditions{iCond} ' '...
          meanPhaseProfiles{iCond}.(fnsData{iProfile}).phaseComparison{1} ' '...
          meanPhaseProfiles{iCond}.(fnsData{iProfile}).phaseComparison{2} ' '...
          meanPhaseProfiles{iCond}.(fnsData{iProfile}).phaseComparison{3}]);
        figFileName = [options.mainFolder filesep options.phaseSumPredictionsSubfolder...
          filesep 'sum_phase_profiles_' conditions{iCond} '_'...
          meanPhaseProfiles{iCond}.(fnsData{iProfile}).phaseComparison{1} '_'...
          meanPhaseProfiles{iCond}.(fnsData{iProfile}).phaseComparison{2} '_'...
          meanPhaseProfiles{iCond}.(fnsData{iProfile}).phaseComparison{3}];
        if ~exist([mainFolder filesep phaseSumPredictionsSubfolder], 'file')
          mkdir([options.mainFolder filesep options.phaseSumPredictionsSubfolder]);
        end
        hgsave(fH{iCond}{iArea}, figFileName);
        
        label = [2 1.6];
        margin = [0.3 0.55];
        width = 1*options.figSize-label(1)-margin(1);
        height = 1*options.figSize-label(2)-margin(2);
        paperSize = resizeFig(gcf, ax1, width, height, label, margin, 0);
        exportFig(fH{iCond}{iArea}, [figFileName '.png'],'-dpng','-r300', paperSize);
      end
    end
  end
end
close all