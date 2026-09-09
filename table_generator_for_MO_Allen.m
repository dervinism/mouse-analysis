% Table generator for MO: Allen data

clear all %#ok<*CLALL>
close all
clc

% IO
phaseFile = 'C:\Users\m329786\Data\infraslow-dynamics\04_data_analysis\006_analysis_results_allen\noRun\area2pupil_comparisons_allensdk\units\globalUnits_area2pupil.mat';
firingRatesFile = 'C:\Users\m329786\Data\infraslow-dynamics\04_data_analysis\006_analysis_results_allen\noRun\area2itself_comparisons_allensdk\units\globalUnits.mat';
outputFile = 'phase_firingRate_table_Allen.mat';

% Parameters
areasOI = {'Th','V1','V2','CA','DG'};
condInd = 1;

% Load the files
phaseData = load(phaseFile);
firingRatesData = load(firingRatesFile);

% Get other parameters
areas = phaseData.areas;
interpFreqs_phase = phaseData.areaFreqFullInterpIndividual;
slowFreqMask_phase = interpFreqs_phase <= 1;
slowFreqs = interpFreqs_phase(slowFreqMask_phase);
nFreqs = sum(slowFreqMask_phase);

% Create table column labels
columnLabels = {'Area', 'Rate', 'rSpearman', 'pvalSpearman'};
for iF = 1:nFreqs
  label = ['Coh_' num2str(slowFreqs(iF))];
  columnLabels{numel(columnLabels)+1} = label;
end
for iF = 1:nFreqs
  label = ['Phi_' num2str(slowFreqs(iF))];
  columnLabels{numel(columnLabels)+1} = label;
end

% Populate the table
nAreas = numel(areas);
for iArea = 1:nAreas
  if ismember(areas{iArea}, areasOI)
    areaPhase = phaseData.areaPhaseFullInterpIndividual{condInd}{iArea}(:,slowFreqMask_phase);
    areaCoh = phaseData.areaCohFullInterpIndividual{condInd}{iArea}(:,slowFreqMask_phase);
    firingRate = firingRatesData.areaFRIndividual{condInd}{iArea};
    rSpearman = firingRatesData.arearSpearman5secIndividual{condInd}{iArea};
    pvalSpearman = firingRatesData.areapvalSpearman5secIndividual{condInd}{iArea};
    unitAreas = repmat(areas(iArea), numel(rSpearman), 1);
    areaUnitTable = table(unitAreas, firingRate, rSpearman, pvalSpearman, ...
        'VariableNames', {'Area', 'Rate', 'rSpearman', 'pvalSpearman'});
    for iF = 1:nFreqs
      label = ['Coh_' strrep(num2str(round(slowFreqs(iF),4)), '.', 'p')];
      areaUnitTable.(label) = areaCoh(:,iF);
    end
    for iF = 1:nFreqs
      label = ['Phi_' strrep(num2str(round(slowFreqs(iF),4)), '.', 'p')];
      areaUnitTable.(label) = areaPhase(:,iF);
    end
    if ~exist('unitTable','var')
      unitTable = areaUnitTable;
    else
      unitTable = [unitTable; areaUnitTable]; %#ok<*AGROW>
    end
  end
end

% Save the table
save(outputFile, 'unitTable', '-v7.3');