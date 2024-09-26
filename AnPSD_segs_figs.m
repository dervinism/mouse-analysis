% A part of the old AnPSD script adapted to produce figures for coherence
% and phase analyses comparing segment spiking rate data.


% LOAD PRE-PROCESSED DATA
if ~exist('dataStruct','var')
  load(dataFile);
end


% INITIALISE PARAMETERS
figsubdirname = '';
UPF = 5; % Units Per Figure
x_lim = [0.02 50];


% PLOT PHASE PREFERENCE FOR SEGMENT COMPARISONS
colours = [[0 0 255];    %blue
           [255 0 0];    %red
           [0 255 0];    %green
           [0 255 255];  %cyan
           [205 205 0];  %yellow
           [255 0 255]]; %magenta
colours = colours./255;


% LOOP THROUGH DB ENTRIES
fnsData = fieldnames(dataStruct);
for dbCount = numel(fnsData):-1:1
  dbStruct = dataStruct.(fnsData{dbCount});
  shankIDs = fieldnames(dbStruct);

% LOOP THROUGH SHANKS
  for sh = 1:numel(shankIDs)
    shankStruct = dbStruct.(shankIDs{sh});
    fprintf('%s shank %d -------------------\n', fnsData{dbCount}, sh);
    
% LOAD THE CONTENTS OF THE SHANK STRUCTURE VARIABLE
    QP = shankStruct.QP;
    eMap = shankStruct.eMap;
    for dbCountShank = 1:numel(shankStruct.db)
      if shankStruct.db(dbCountShank).entryName == fnsData{dbCount}
        entryName = shankStruct.db(dbCountShank).entryName;
        break
      end
    end
    
    for p = 1:numel(QP)
      phase = QP{p}.phase;
      countSignificant = length(phase(~isnan(QP{p}.phase_confU)));
      countTotal = length(phase);
      titleStr = ['_' num2str(countSignificant) '_' num2str(countTotal) sprintf('_mfr1_%4.1f_mfr2_%4.1f', QP{p}.mfr1, QP{p}.mfr2)];
      figStr = ['_' num2str(countSignificant) '_' num2str(countTotal)];
      
      prefix = [entryName '_pair' num2str(QP{p}.pair)];
      figH = phaseFigPairMouse(QP{p}, eMap, colours, titleStr, figStr, prefix, 'off');
      
      close(figH)
    end
  end
end
  