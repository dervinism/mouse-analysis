% Run this script to perform analyses comparing LFP and pupil area data.

fclose all;
close all
clear
clc



%% INITIALISE PARAMETERS
dataDir = 'R:\CSN\Shared\Dynamics\Data\mouse_analysis_data';
lists



%% DISPLAY RIPPLE RATE FREQUENCY PHASE PROFILES FOR ALL ANIMALS AND ALL RECORDINGS
% LOOP THROUGH DB ENTRIES
recCount = 0;
recSignificant = 0;
for dbCount = 1:numel(fullSeries)
  strSep = strfind(fullSeries{dbCount},'s');
  series = fullSeries{dbCount}(strSep+1:end);
  animal = fullSeries{dbCount}(1:strSep-2);
  animalColour = animalColours(animal);
  if dbCount == 1 || ~strcmpi(animal,prevAnimal)
    load([dataDir filesep animal filesep 'CAR' filesep animal '.mat'])
    updateLegend = true;
  end
  prevAnimal = animal;
  
% LOAD THE CONTENTS OF THE DB STRUCTURE VARIABLE
  dbStruct = dataStruct.seriesData.(fullSeries{dbCount});
  if ~isfield(dbStruct,'lfpphaseCohDataPR')
    disp(['No phase data for ' fullSeries{dbCount} '. Skipping to the next db entry...']);
    continue
  end
  recCount = recCount + 1;
  
% PLOT THE DATA
  phase = bestUnwrap(dbStruct.lfpphaseCohDataPR.rippleRate{ch{dbCount}}.phase);
  phase(isnan(dbStruct.lfpphaseCohDataPR.rippleRate{ch{dbCount}}.phase_confU)) = NaN;
  phase(isnan(dbStruct.lfpphaseCohDataPR.rippleRate{ch{dbCount}}.phase_confL)) = NaN;
  phase(isnan(dbStruct.lfpphaseCohDataPR.rippleRate{ch{dbCount}}.coh_conf)) = NaN;
  freq = dbStruct.lfpphaseCohDataPR.rippleRate{ch{dbCount}}.freq;
  [~, endFreq] = min(abs(freq - 2));
%   phase1 = phase(1:endFreq);
%   if sum(isnan(phase1)) > 0.6*numel(phase1)
%     continue
%   end
  recSignificant = recSignificant + 1;
  if mean(phase, 'omitnan') > 2*pi
    phase = phase - 2*pi;
  end
  if mean(phase, 'omitnan') < -(0.15)*pi
    phase = phase + 2*pi;
  end
  if ~exist('figH','var')
    figH = figure('units', 'normalized', 'position', [0.002, .04, 1, .88], 'Visible', 'on');
    semilogx([10e-4 10e2],[0 0], 'k:'); hold on
    p = semilogx(freq,phase, 'Color',animalColour);
    pLegend = p;
    txtLegend = {animal};
    updateLegend = false;
  else
    if updateLegend
      p = semilogx(freq,phase, 'Color',animalColour);
      pLegend = [pLegend p]; %#ok<AGROW>
      txtLegend{numel(txtLegend)+1} = animal;
      updateLegend = false;
    else
      p = semilogx(freq,phase, 'Color',animalColour);
    end
  end
end


xlabel('Frequency (Hz)');
xlim([10e-3 - 0.002   10e1 + 200])
ylabel('Phase (rad)');
legend(pLegend,txtLegend, 'Interpreter','None');
title(['Ripple rate vs population firing rate.' ' Significant recordings: ' num2str(recSignificant) '/' num2str(recCount)]);
set(figH, 'Name','Ripple rate vs population firing rate');
hgsave(figH, 'rippleRate2mua');
close(figH);