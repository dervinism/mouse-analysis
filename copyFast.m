% filenames = {                        'determineArea',                     'lfpLoad',                           ...
%                                      'determineAreaFromSeries',                                                ...
%                                      'determineAreaFromSeriesOld',        'lfpSaturations_allensdk',           ...
%                                      'determineAreaOld',                  'lfpSaturations_uol',                ...
% 'AnHT',                              'determineCondition',                'lists',                             ...
% 'AnHT_plot',                         'determineInds',                     'loadAll',                           ...
%                                      'determineSeries',                   'loadAllen',                         ...
% 'AnPSD_load',                        'dispFOI',                           'loadAsMUA_noResClu',                ...
% 'AnPSD_load_subpop',                 'dispPhaseHistoStats',               'loadAsRaster',                      ...
% 'AnPSD_load_subpopFilter',           'displayVideoFrames',                'loadAsRasterSparse',                ...
% 'AnPSD_load_subpop_new',             'duplicateComparison',               'loadSpikes_allensdk',               ...
% 'AnPSD_load_waveforms',              'duplicateTest',                     'loadUOL',                           ...
% 'AnPSD_muas',                        'electrodeMap',                                                           ...
% 'AnPSD_old',                         'ellipseVdlc',                                                            ...
% 'AnPSD_old_old',                     'exceptionTest',                     'meanPhaseProfilesSummary',          ...
% 'AnPSD_segs',                        'exclusionTest',                     'meanPhaseProfilesSummary2',         ...
% 'AnPSD_segs_figs',                   'expandFields_ca',                   'meanPhaseProfilesSummary3',         ...
% 'AnPSD_units',                       'eyeAnalysis',                       'meanPhaseProfilesSummary4',         ...
% 'AnPSD_unitsLaminar',                'eyeAnalysis_figs',                  'meanVarCalc',                       ...
% 'AnPSD_unitsLaminar_negative',       'eyeAnalysis_muas',                  'motionAnalysis',                    ...
% 'AnPSD_unitsLaminar_positive',                                            'motionAnalysis_figs',               ...
% 'AnPSD_units_ca',                    'eyeLoad',                                                                ...
% 'AnPSD_units_ca_negative',           'eyeLoad_allensdk',                  'motionFilt',                        ...
% 'AnPSD_units_ca_neutral',            'eyeLoad_filter',                    'motionLoad',                        ...
% 'AnPSD_units_ca_positive',           'firingRateTest',                    'moveAll',                           ...
% 'AnPSD_units_figs',                  'freqDependentWindowCoherenceMD',    'multiViolinPlots',                  ...
% 'AnPSD_units_figs_ca',               'getLog',                            'params',                            ...
% 'AnPSD_units_figs_ca_negative',      'get_dbStruct',                                                           ...
% 'AnPSD_units_figs_ca_neutral',       'get_shankStruct',                   'peakPhaseProfilesSummary',          ...
% 'AnPSD_units_figs_ca_positive',      'globalBetaAnalysis',                'phaseCohCalc',                      ...
% 'PCAanalysis',                       'globalCoherenceAnalysis',           'phaseCoherencePlots',               ...
% 'PCAanalysis2',                      'globalConditionsAnalysis',          'phaseCoherencePlots_ca',            ...
% 'PCAanalysis2_noS1',                 'globalConditionsEyeAnalysis',       'phaseComparisonStats',              ...
% 'PCAanalysis2_noS1_noRSC',           'globalEyeAnalysis',                 'phaseFigPairMouse',                 ...
% 'PCAanalysis2_noS1_noRSC_noVB',      'globalEyeAnalysisLFP',              'phaseFreqProfilePlotIndividual',    ...
% 'PCAanalysis2_noS1_noVB',            'globalEyeAnalysis_halves',          'phaseFreqProfilePlotMeans',         ...
% 'PCAanalysis3_motion',               'globalEyeAnalysis_units',           'phaseFreqProfilePlotUpdate',        ...
% 'PCAanalysis3_pupil',                'globalEyeAnalysis_unitsHalves',     'phaseFreqProfilePlotUpdate_units',  ...
% 'PCAanalysis4',                      'globalEyeAnalysis_unitsLocs',       'phaseGraph',                        ...
%                                      'globalFR',                          'phaseHistPlot',                     ...
% 'acgdPlotUnitMeans',                 'globalFRPhaseAnalysis',             'phaseHistPlotUnadjusted',           ...
% 'acgdPlotUnits',                     'globalFRPhaseAnalysis2',            'phaseHistSubplotUnadjusted',        ...
% 'addDependencies',                   'globalLfpAnalysis',                 'phaseHistosPlotMaster',             ...
% 'addParpoolDependencies',            'globalLfpAnalysisCoh',              'phaseLocationProfile',              ...
% 'adjustPhase',                       'globalLfpAnalysisMotion',           'phaseLocationProfileSummary',       ...
% 'adjustPi',                          'globalLfpAnalysisMotionCoh',        'phaseVphaseSubplot',                ...
% 'allensdkDB',                        'globalLfpAnalysisPR',               'pickChan',                          ...
% 'animalColours',                     'globalLfpAnalysisPR_ca_ripples',    'plotAssign',                        ...
% 'areaColours',                       'globalLfpAnalysisPRcoh',            'plotEyeOrMotion',                   ...
% 'areaColours2',                      'globalMotionAnalysis',              'plotFiringRate',                    ...
% 'areaEntry',                         'globalMotionAnalysis_halves',       'plotUpdateSeries',                  ...
% 'areaNames',                         'globalMotionAnalysis_units',        'psdCalc',                           ...
% 'bimodalitySort',                    'globalMotionAnalysis_unitsHalves',  'psdChannels',                       ...
% 'bimodalitySort2',                   'globalMotionAnalysis_unitsLocs',    'pupilCorrFiringRatePlot',           ...
% 'circPeak',                          'globalPCA',                         'pupilCorrFractionPlot',             ...
% 'cleanUp',                           'globalPCA2',                        'pupilCorrFractionPlot2',            ...
% 'cleanupExcept',                     'globalPCA2_noS1_noRSC',             'pupilCorrPhasePlot',                ...
% 'cohFreqProfilePlotIndividual',      'globalPCA4_FB',                     'pupilDynamics',                     ...
% 'cohFreqProfilePlotMeans',           'globalPCA4_FG',                     'pupilFilt',                         ...
% 'cohFreqProfilePlotUpdate',          'globalPCA4_HSB',                    'pupilMovieProcessingMD',            ...
% 'cohFreqProfilePlotUpdate_units',    'globalPCA4_LFP',                    'pupilPlotMulti',                    ...
% 'cohHistPlot',                       'globalPCA4_NSB',                    'qualityTest',                       ...
% 'cohMatrix',                         'globalPCA4_PR',                     'qualityTest2',                      ...
% 'cohVcohSubplot',                    'globalPCA4_SB',                     'rasterPlotMultiple',                ...
% 'combineData',                       'globalPCA4_SG',                     'rasterPlotOI',                      ...
% 'combineEyeData',                    'globalPCA4_UF',                     'rasterPlotSingle',                  ...
% 'combineMotionData',                 'globalPCA4_UFR',                    'rateCalc',                          ...
% 'combinePeriods',                    'globalPCA4_alpha',                  'recQuality',                        ...
% 'combineRuns',                       'globalPCA4_beta',                   'rescaleColour',                     ...
% 'commonDB',                          'globalPCA4_delta',                  'rippleTest',                        ...
% 'compareHalves',                     'globalPCA4_gamma',                  'runAll',                            ...
% 'compareHalves_ca',                  'globalPCA4_theta',                  'runAllen',                          ...
% 'compareHalves_ca_negative',         'globalPhaseAnalysis',               'runUOL',                            ...
% 'compareHalves_ca_neutral',          'globalPhaseAnalysisLFP',            'runUOLNp',                          ...
% 'compareHalves_ca_positive',         'globalPhaseAnalysis_halves',        'runUOLNx',                          ...
% 'compareHalves_figs',                'globalPhaseAnalysis_units',         'runningSpeedLoad',                  ...
% 'compareHalves_figs_ca',             'globalPhaseAnalysis_unitsHalves',   'saveFigsState',                     ...
% 'compareHalves_figs_ca_negative',    'globalPhaseAnalysis_unitsLaminar',  'saveFigsStateMFR',                  ...
% 'compareHalves_figs_ca_neutral',     'globalPhaseAnalysis_unitsLocal',    'searchStringArray',                 ...
% 'compareHalves_figs_ca_positive',    'globalPhaseAnalysis_unitsLocs',     'series2Comparison',                 ...
% 'compareSeries',                     'globalUnitsFigs',                   'series2compare',                    ...
% 'compareSeriesMultiple',             'globalWaveformAnalysis',            'series2condition',                  ...
% 'compareSeries_figs',                'halfCohCorr',                       'seriesFromEntry',                   ...
% 'comparison2areas',                  'halfCohCorrPlot_ca',                'seriesNames',                       ...
% 'concatenateSeries',                 'halfCohCorrSummary',                'shiftRadAxis',                      ...
% 'convertAll',                        'halfCorrPlot',                      'spkCoh',                            ...
% 'convertAllen',                      'halfCorrPlotMFR',                   'spkPhase',                          ...
% 'convertBetaWindows',                'halfCorrPlot_coh',                  'spkPhaseSignificant',               ...
% 'convertFormat',                     'halfCorrSummary_coh',               'stprCalc',                          ...
% 'convertUOL',                        'halfPhaseCorr',                     'stprHalfCalc',                      ...
% 'copyAll',                           'halfPhaseCorrPlot_ca',              'sumPhases',                         ...
% 'copyAllen',                         'halfPhaseCorrSummary',              'summaryFigs',                       ...
% 'copyFast',                          'histPlotFR',                        'syncFuncWrap',                      ...
% 'copyMouseAnalysisData',                                                  'tidySaveFig',                       ...
% 'copyNSG',                           'initFigsSeries',                    'unitCount',                         ...
% 'copyUOL',                           'initLegendsSeries',                 'unitMUARasterPupilPlot',            ...
% 'corrChannels',                      'initPCA',                           'unitMUARasterPupilPlotOI',          ...
% 'corrStates',                        'initPool',                          'unitMUApupilPlot',                  ...
% 'corrStatesMFR',                     'joinData',                          'unitMUApupilPlotOI',                ...
% 'correctEntries',                    'kappaCalc',                         'unitNamer',                         ...
% 'correctPhaseCoh',                   'lfpAnalysis',                       'unitPos',                           ...
% 'correctlfpLoad',                    'lfpAnalysisLFP_ca',                 'updateDB',                          ...
%                                      'lfpAnalysisMotion',                 'updateDataStruct',                  ...
%                                      'lfpAnalysisPR',                     'updateLegendsSeries',               ...
%                                      'lfpAnalysisPR_ca',                  'updateSeriesData',                  ...
% 'detectFrames',                      'lfpDataSplit',                      'xySubplot',                         ...
% 'detectFramesWrapMotion',            'lfpDataSplitFull',                                                       ...
% 'detectFramesWrapPupil',             'lfpDataSplitRipples',               'AnPSD_units_ca_muas',               ...
% 'AnPSD_units_ca_muas',               'AnPSD_units_ca_muas_positive',      'AnPSD_units_ca_muas_positive',      ...
% 'AnPSD_units_ca_muas_negative',      'AnPSD_units_ca_muas_negative',      'globalPhaseAnalysis_muas',          ...
% 'globalPhaseAnalysis_muasHalves',    'summaryFigsMUAs',                   'globalFRPhaseAnalysis_muas',        ...
% 'globalEyeAnalysis_muasHalves',      'globalEyeAnalysis_muas',            'globalPhaseAnalysis_MUAsLaminar',   ...
% 'summaryStatsPhase',                 'dispPhaseHistoStats_ca',            'phaseComparisonStats_ca',           ...
% 'summaryStatsPupilHalfPhase',        'summaryStatsPupilPhase',            'summaryStatsFRdots',                ...
% 'summaryStatsFR'};
filenames = {                        'determineArea',                     'lfpLoad',                           ...
                                     'determineAreaFromSeries',                                                ...
                                     'determineAreaFromSeriesOld',        'lfpSaturations_allensdk',           ...
                                     'determineAreaOld',                  'lfpSaturations_uol',                ...
'AnHT',                              'determineCondition',                'lists',                             ...
'AnHT_plot',                         'determineInds',                     'loadAll',                           ...
                                     'determineSeries',                   'loadAllen',                         ...
'AnPSD_load',                        'dispFOI',                           'loadAsMUA_noResClu',                ...
'AnPSD_load_subpop',                 'dispPhaseHistoStats',               'loadAsRaster',                      ...
'AnPSD_load_subpopFilter',           'displayVideoFrames',                'loadAsRasterSparse',                ...
'AnPSD_load_subpop_new',             'duplicateComparison',               'loadSpikes_allensdk',               ...
'AnPSD_load_waveforms',              'duplicateTest',                     'loadUOL',                           ...
'AnPSD_muas',                        'electrodeMap',                                                           ...
'AnPSD_old',                         'ellipseVdlc',                                                            ...
'AnPSD_old_old',                     'exceptionTest',                     'meanPhaseProfilesSummary',          ...
'AnPSD_segs',                        'exclusionTest',                     'meanPhaseProfilesSummary2',         ...
'AnPSD_segs_figs',                   'expandFields_ca',                   'meanPhaseProfilesSummary3',         ...
'AnPSD_units',                       'eyeAnalysis',                       'meanPhaseProfilesSummary4',         ...
'AnPSD_unitsLaminar',                'eyeAnalysis_figs',                  'meanVarCalc',                       ...
'AnPSD_unitsLaminar_negative',       'eyeAnalysis_muas',                  'motionAnalysis',                    ...
'AnPSD_unitsLaminar_positive',                                            'motionAnalysis_figs',               ...
'AnPSD_units_ca',                    'eyeLoad',                                                                ...
'AnPSD_units_ca_negative',           'eyeLoad_allensdk',                  'motionFilt',                        ...
'AnPSD_units_ca_neutral',            'eyeLoad_filter',                    'motionLoad',                        ...
'AnPSD_units_ca_positive',           'firingRateTest',                    'moveAll',                           ...
'AnPSD_units_figs',                  'freqDependentWindowCoherenceMD',    'multiViolinPlots',                  ...
'AnPSD_units_figs_ca',               'getLog',                                                                 ...
'AnPSD_units_figs_ca_negative',      'get_dbStruct',                                                           ...
'AnPSD_units_figs_ca_neutral',       'get_shankStruct',                   'peakPhaseProfilesSummary',          ...
'AnPSD_units_figs_ca_positive',      'globalBetaAnalysis',                'phaseCohCalc',                      ...
'PCAanalysis',                       'globalCoherenceAnalysis',           'phaseCoherencePlots',               ...
'PCAanalysis2',                      'globalConditionsAnalysis',          'phaseCoherencePlots_ca',            ...
'PCAanalysis2_noS1',                 'globalConditionsEyeAnalysis',       'phaseComparisonStats',              ...
'PCAanalysis2_noS1_noRSC',           'globalEyeAnalysis',                 'phaseFigPairMouse',                 ...
'PCAanalysis2_noS1_noRSC_noVB',      'globalEyeAnalysisLFP',              'phaseFreqProfilePlotIndividual',    ...
'PCAanalysis2_noS1_noVB',            'globalEyeAnalysis_halves',          'phaseFreqProfilePlotMeans',         ...
'PCAanalysis3_motion',               'globalEyeAnalysis_units',           'phaseFreqProfilePlotUpdate',        ...
'PCAanalysis3_pupil',                'globalEyeAnalysis_unitsHalves',     'phaseFreqProfilePlotUpdate_units',  ...
'PCAanalysis4',                      'globalEyeAnalysis_unitsLocs',       'phaseGraph',                        ...
                                     'globalFR',                          'phaseHistPlot',                     ...
'acgdPlotUnitMeans',                 'globalFRPhaseAnalysis',             'phaseHistPlotUnadjusted',           ...
'acgdPlotUnits',                     'globalFRPhaseAnalysis2',            'phaseHistSubplotUnadjusted',        ...
'addDependencies',                   'globalLfpAnalysis',                 'phaseHistosPlotMaster',             ...
'addParpoolDependencies',            'globalLfpAnalysisCoh',              'phaseLocationProfile',              ...
'adjustPhase',                       'globalLfpAnalysisMotion',           'phaseLocationProfileSummary',       ...
'adjustPi',                          'globalLfpAnalysisMotionCoh',        'phaseVphaseSubplot',                ...
'allensdkDB',                        'globalLfpAnalysisPR',               'pickChan',                          ...
'animalColours',                     'globalLfpAnalysisPR_ca_ripples',    'plotAssign',                        ...
'areaColours',                       'globalLfpAnalysisPRcoh',            'plotEyeOrMotion',                   ...
'areaColours2',                      'globalMotionAnalysis',              'plotFiringRate',                    ...
'areaEntry',                         'globalMotionAnalysis_halves',       'plotUpdateSeries',                  ...
'areaNames',                         'globalMotionAnalysis_units',        'psdCalc',                           ...
'bimodalitySort',                    'globalMotionAnalysis_unitsHalves',  'psdChannels',                       ...
'bimodalitySort2',                   'globalMotionAnalysis_unitsLocs',    'pupilCorrFiringRatePlot',           ...
'circPeak',                          'globalPCA',                         'pupilCorrFractionPlot',             ...
'cleanUp',                           'globalPCA2',                        'pupilCorrFractionPlot2',            ...
'cleanupExcept',                     'globalPCA2_noS1_noRSC',             'pupilCorrPhasePlot',                ...
'cohFreqProfilePlotIndividual',      'globalPCA4_FB',                     'pupilDynamics',                     ...
'cohFreqProfilePlotMeans',           'globalPCA4_FG',                     'pupilFilt',                         ...
'cohFreqProfilePlotUpdate',          'globalPCA4_HSB',                    'pupilMovieProcessingMD',            ...
'cohFreqProfilePlotUpdate_units',    'globalPCA4_LFP',                    'pupilPlotMulti',                    ...
'cohHistPlot',                       'globalPCA4_NSB',                    'qualityTest',                       ...
'cohMatrix',                         'globalPCA4_PR',                     'qualityTest2',                      ...
'cohVcohSubplot',                    'globalPCA4_SB',                     'rasterPlotMultiple',                ...
'combineData',                       'globalPCA4_SG',                     'rasterPlotOI',                      ...
'combineEyeData',                    'globalPCA4_UF',                     'rasterPlotSingle',                  ...
'combineMotionData',                 'globalPCA4_UFR',                    'rateCalc',                          ...
'combinePeriods',                    'globalPCA4_alpha',                  'recQuality',                        ...
'combineRuns',                       'globalPCA4_beta',                   'rescaleColour',                     ...
'commonDB',                          'globalPCA4_delta',                  'rippleTest',                        ...
'compareHalves',                     'globalPCA4_gamma',                  'runAll',                            ...
'compareHalves_ca',                  'globalPCA4_theta',                  'runAllen',                          ...
'compareHalves_ca_negative',         'globalPhaseAnalysis',               'runUOL',                            ...
'compareHalves_ca_neutral',          'globalPhaseAnalysisLFP',            'runUOLNp',                          ...
'compareHalves_ca_positive',         'globalPhaseAnalysis_halves',        'runUOLNx',                          ...
'compareHalves_figs',                'globalPhaseAnalysis_units',         'runningSpeedLoad',                  ...
'compareHalves_figs_ca',             'globalPhaseAnalysis_unitsHalves',   'saveFigsState',                     ...
'compareHalves_figs_ca_negative',    'globalPhaseAnalysis_unitsLaminar',  'saveFigsStateMFR',                  ...
'compareHalves_figs_ca_neutral',     'globalPhaseAnalysis_unitsLocal',    'searchStringArray',                 ...
'compareHalves_figs_ca_positive',    'globalPhaseAnalysis_unitsLocs',     'series2Comparison',                 ...
'compareSeries',                     'globalUnitsFigs',                   'series2compare',                    ...
'compareSeriesMultiple',             'globalWaveformAnalysis',            'series2condition',                  ...
'compareSeries_figs',                'halfCohCorr',                       'seriesFromEntry',                   ...
'comparison2areas',                  'halfCohCorrPlot_ca',                'seriesNames',                       ...
'concatenateSeries',                 'halfCohCorrSummary',                'shiftRadAxis',                      ...
'convertAll',                        'halfCorrPlot',                      'spkCoh',                            ...
'convertAllen',                      'halfCorrPlotMFR',                   'spkPhase',                          ...
'convertBetaWindows',                'halfCorrPlot_coh',                  'spkPhaseSignificant',               ...
'convertFormat',                     'halfCorrSummary_coh',               'stprCalc',                          ...
'convertUOL',                        'halfPhaseCorr',                     'stprHalfCalc',                      ...
'copyAll',                           'halfPhaseCorrPlot_ca',              'sumPhases',                         ...
'copyAllen',                         'halfPhaseCorrSummary',              'summaryFigs',                       ...
'copyFast',                          'histPlotFR',                        'syncFuncWrap',                      ...
'copyMouseAnalysisData',                                                  'tidySaveFig',                       ...
'copyNSG',                           'initFigsSeries',                    'unitCount',                         ...
'copyUOL',                           'initLegendsSeries',                 'unitMUARasterPupilPlot',            ...
'corrChannels',                      'initPCA',                           'unitMUARasterPupilPlotOI',          ...
'corrStates',                        'initPool',                          'unitMUApupilPlot',                  ...
'corrStatesMFR',                     'joinData',                          'unitMUApupilPlotOI',                ...
'correctEntries',                    'kappaCalc',                         'unitNamer',                         ...
'correctPhaseCoh',                   'lfpAnalysis',                       'unitPos',                           ...
'correctlfpLoad',                    'lfpAnalysisLFP_ca',                 'updateDB',                          ...
                                     'lfpAnalysisMotion',                 'updateDataStruct',                  ...
                                     'lfpAnalysisPR',                     'updateLegendsSeries',               ...
                                     'lfpAnalysisPR_ca',                  'updateSeriesData',                  ...
'detectFrames',                      'lfpDataSplit',                      'xySubplot',                         ...
'detectFramesWrapMotion',            'lfpDataSplitFull',                                                       ...
'detectFramesWrapPupil',             'lfpDataSplitRipples',               'AnPSD_units_ca_muas',               ...
'AnPSD_units_ca_muas',               'AnPSD_units_ca_muas_positive',      'AnPSD_units_ca_muas_positive',      ...
'AnPSD_units_ca_muas_negative',      'AnPSD_units_ca_muas_negative',      'globalPhaseAnalysis_muas',          ...
'globalPhaseAnalysis_muasHalves',    'summaryFigsMUAs',                   'globalFRPhaseAnalysis_muas',        ...
'globalEyeAnalysis_muasHalves',      'globalEyeAnalysis_muas',            'globalPhaseAnalysis_MUAsLaminar',   ...
'summaryStatsPhase',                 'dispPhaseHistoStats_ca',            'phaseComparisonStats_ca',           ...
'summaryStatsPupilHalfPhase',        'summaryStatsPupilPhase',            'summaryStatsFRdots',                ...
'summaryStatsFR'};
% filenames = {'compareHalves_figs_ca_neutral','compareHalves_ca_neutral','AnPSD_units_figs_ca_neutral','AnPSD_units_ca_neutral',...
%   'compareHalves_figs_ca_negative','compareHalves_ca_negative','AnPSD_units_figs_ca_negative','AnPSD_units_ca_negative',...
%   'compareHalves_figs_ca_positive','compareHalves_ca_positive','AnPSD_units_figs_ca_positive','AnPSD_units_ca_positive',...
%   'compareHalves_figs_ca','compareHalves_ca','AnPSD_units_figs_ca','AnPSD_units_ca','AnPSD_unitsLaminar_negative',...
%   'AnPSD_unitsLaminar_positive','AnPSD_load_subpopFilter','AnPSD_load_subpop_new','motionAnalysis_figs','motionAnalysis',...
%   'eyeAnalysis_muas','eyeAnalysis_figs','eyeAnalysis','compareSeries_figs','compareSeries','compareHalves_figs','compareHalves',...
%   'AnPSD_unitsLaminar','AnPSD_muas','AnPSD_units_figs','AnPSD_units','runningSpeedLoad','motionLoad','eyeLoad_filter',...
%   'eyeLoad','AnPSD_load_waveforms','AnPSD_load','AnHT_plot','AnHT','AnPSD_segs_figs','AnPSD_segs','phaseCohCalc','psdCalc',...
%   'psdChannels','pupilDynamics','freqDependentWindowCoherenceMD','params','lfpAnalysis','lfpAnalysisLFP_ca','lfpAnalysisMotion',...
%   'lfpAnalysisPR','lfpAnalysisPR_ca','psdChannels'};
% filenames = {'AnPSD_load_subpop_new'};
destination = 'S:\cortical_dynamics\User\md406\runNSG';

for iFile = 1:numel(filenames)
  disp(iFile);
  source = ['S:\cortical_dynamics\User\md406\code\mouse_analysis' filesep filenames{iFile}];
  if ~isfolder(source)
    source = [source '.m']; %#ok<*AGROW>
  end
  folders = dir(destination);
  for f = 3:numel(folders)
    if ~isfolder(source)
      copyfile(source, [destination filesep folders(f).name])
    else
      copyfile(source, [destination filesep folders(f).name filesep filenames{iFile}])
    end
  end
end