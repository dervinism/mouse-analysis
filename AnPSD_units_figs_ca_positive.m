% A part of the old AnPSD script adapted to display and save figures of
% coherence and phase analyses for unit vs population spiking rate of a
% different area (within positive subpopulation only).


%% LOAD PRE-PROCESSED DATA
if ~exist('dataStruct','var')
  load(dataFile);
end


%% INITIALISE PARAMETERS
params
lists
qualityCheck = false;
figsubdirname = '';
UPF = 5; % Units Per Figure
x_lim = [0.01 120];
visibility = 'on';
figType = 'ephysUnit';


%% DRAW FIGURES
fnsData = fieldnames(dataStruct.seriesData_positive);
if (~isempty(dbEntries) && dbEntries(1) == inf) ||...
    (~isempty(dbEntries) && dbEntries(end) > numel(fnsData))
  dbEntries_positive = 1:numel(fnsData);
else
  dbEntries_positive = dbEntries;
end
for dbCount = dbEntries_positive % Loop through db entries
  
  % Load the contents of dbStruct
  [~, ~, ~, entryName, ~, ~, shankIDs,...
    ~, ~, ~, ~, srData, ~, ~, FOI] = get_dbStruct(dataStruct, dbCount, 'positive');
  
  fnsData_ca = fieldnames(dataStruct.seriesData_ca_positive);
  for sca = 1:numel(fnsData_ca) % Loop through comparison areas
    sca_name = fnsData_ca{sca};
    ind = strfind(sca_name, '__');
    if ~strcmpi(entryName, sca_name(1:ind-1))
      continue
    end
    dbStruct_ca = dataStruct.seriesData_ca_positive.(sca_name);
    if isfield (dbStruct_ca, 'shankData')
      shankIDs = fieldnames(dbStruct_ca.shankData);
    else
      continue
    end
    phaseFOI = dbStruct_ca.popData.phaseCoh.phaseFOI;
    cohFOI = dbStruct_ca.popData.phaseCoh.cohFOI;
    
    % Output folder
    strSep = strfind(sca_name,'__');
    sca_name1 = sca_name(1:strSep-1);
    strSep1 = strfind(sca_name1,'s');
    sca_name1 = sca_name1(strSep1+1:end);
    sca_name2 = sca_name(strSep+2:end);
    strSep2 = strfind(sca_name2,'s');
    sca_name2 = sca_name2(strSep2+1:end);
    figsubdirname = ['positive' filesep sca_name1 'v' sca_name2];
    if ~exist(figsubdirname,'dir')
      mkdir(figsubdirname)
    end
    
    for sh = 1:numel(shankIDs) % Loop through shanks
      shankStruct = eval(['dbStruct_ca.shankData.' shankIDs{sh}]);
      fprintf('%s shank %d -------------------\n', fnsData{dbCount}, sh);
      
      % Load the contents of shankStruct
      shank = shankStruct.shankID;
      units = shankStruct.units;
      if isempty(units)
        continue
      end
      unitMetadata = shankStruct.unitMetadata;
      phaseCoh = shankStruct.phaseCoh;
      
      nUnits = numel(units);
      for u = 1:nUnits % Loop through units
        fprintf('Processing %s shank %d unit %d\n', fnsData{dbCount}, sh, u);
        unitData_ca = phaseCoh{u};
        if ~sum(isnan(unitData_ca.freq))
          phaseCoherencePlots_ca(sh, u, UPF, unitData_ca, srData, units(u), nUnits, x_lim, figsubdirname, sca_name, figType, visibility, false);
        end
      end % loop over units
    end % loop over shanks

    % Compare population firing rates
    if ~sum(isnan(dbStruct_ca.popData.phaseCohPop.freq))
      phaseCoherencePlots_ca(0, 1, 1, dbStruct_ca.popData.phaseCohPop, srData, 'mua', 1, x_lim, figsubdirname, sca_name, figType, visibility, true);
    else
      continue
    end
    
    if ~isempty(dbStruct_ca.popData.phaseCohPopUnitOnly)
      if ~sum(isnan(dbStruct_ca.popData.phaseCohPopUnitOnly.freq))
        phaseCoherencePlots_ca(0, 1, 1, dbStruct_ca.popData.phaseCohPopUnitOnly, srData, 'mua_unitOnly', 1, x_lim, figsubdirname, sca_name, figType, visibility, true);
      end
    end
    
    % Check unit quality
    if qualityCheck
      qualityInds = unitMetadata(:,6) <= refractCont & unitMetadata(:,7) >= cluDist;
      qualityUnits = unitMetadata(qualityInds,1);
    else
      qualityInds = 1:numel(units);
      qualityUnits = units;
    end
    phaseFOI_reduced = phaseFOI(qualityInds,:);
    cohFOI_reduced = cohFOI(qualityInds,:);
    
    % Draw phase histograms for all units
    edges = -pi:pi/8:pi; %#ok<*UNRCH>
    edgesC = -pi+pi/16:pi/8:pi-pi/16;
    histPhaseOverF = zeros(numel(edges),numel(FOI));
    histPhaseOverFC = zeros(numel(edgesC),numel(FOI));
    for f = 1:numel(FOI)
      histPhaseOverF(:,f) = [sum(isnan(phaseFOI_reduced(:,f))) histcounts(phaseFOI_reduced(:,f), edges)];
      histPhaseOverFC(:,f) = [sum(isnan(phaseFOI_reduced(:,f))) histcounts(phaseFOI_reduced(:,f), edgesC)];
    end
    histPhaseOverF = [histPhaseOverF(1,:); adjustPi(histPhaseOverF(2:end,:))];
    histPhaseOverFC = [histPhaseOverFC(1,:); adjustPi(histPhaseOverFC(2:end,:))]; %#ok<*NASGU>
    for f = 1:numel(FOI)
      phaseHistPlot(edges, histPhaseOverF(:,f), '# neurons', [figsubdirname filesep sca_name '_' figType '_phase_at_' num2str(FOI(f)) '_Hz'], visibility);
      %phaseHistPlot(edgesC, histPhaseOverFC(:,f), '# neurons', [figsubdirname filesep 'phaseC_at_' num2str(FOI(f)) '_Hz'], visibility);
    end
    
    % Draw coherence histograms for all units
    edges1 = 0:0.05:1;
    edges2 = 0:0.025:0.5;
    separator = 12;
    histCohOverF1 = zeros(numel(edges1),numel(FOI(separator+1:end)));
    histCohOverF2 = zeros(numel(edges2),numel(FOI(1:separator)));
    for f = 1:size(histCohOverF1,2)
      histCohOverF1(:,f) = [sum(isnan(cohFOI_reduced(:,separator+f))) histcounts(cohFOI_reduced(:,separator+f), edges1)];
    end
    for f = 1:size(histCohOverF1,2)
      cohHistPlot(edges1, histCohOverF1(:,f), {'non-signif.','0','0.25','0.5','0.75','1'}, '# neurons',...
        [figsubdirname filesep sca_name '_' figType '_coherence_at_' num2str(FOI(separator+f)) '_Hz'], visibility);
    end
    for f = 1:size(histCohOverF2,2)
      histCohOverF2(:,f) = [sum(isnan(cohFOI_reduced(:,f))) histcounts(cohFOI_reduced(:,f), edges2)];
    end
    for f = 1:size(histCohOverF2,2)
      cohHistPlot(edges2, histCohOverF2(:,f), {'non-signif.','0','0.125','0.25','0.375','0.5'}, '# neurons',...
        [figsubdirname filesep sca_name '_' figType '_coherence_at_' num2str(FOI(f)) '_Hz'], visibility);
    end
    
  end % loop over comparison areas
end % loop over db entries

clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca dataStruct