function [MUAsAll, spkMUAs, muas, muaMap, MUAsAll_units, spkUnits, units, unitMap] = loadSpikes_allensdk(headerUnits, unitMetadata, headerSpikes, spikeData, area, refractCont, cluDist, count, sr)
% [MUAsAll, spkMUAs, muas, muaLoc, MUAsAll_units, spkUnits, units, unitMap] = loadSpikes_allensdk(headerUnits, unitMetadata, headerSpikes, spikeData, area, refractCont, cluDist, count, sr)
%
% Function loads allensdk spikes into sparse spiking matrices.
% Input: headerUnits - unit metadata header.
%        unitMetadata - unit metadata matrix.
%        headerSpikes - spiking data header.
%        spikeData - spike time data.
%        area - a brain area of interest. The function will only load
%               spikes that were recorded in this area. If this input is
%               left empty, the function will load all data irrespective of
%               brain area.
%        refractCont - ISI violations criterion for distinguishing between
%                      MUAs and units. Allensdk data comes without unit and
%                      MUA labels. Therefore, this criterion is used to
%                      help distinguish between units and MUAs. Other
%                      criteria are also used. For more details look at
%                      https://allensdk.readthedocs.io/en/latest/_static/examples/nb/ecephys_quality_metrics.html
%        cluDist - cluster isolation distance (an unit quality extra metric
%                  to distinguish between units and MUAs).
%        count - minimum spike count. Units with fewer spikes are excluded.
%        sr - sampling rate in Hz of the output vectors and matrices.
% Output: MUAsAll - population firing rate for your area of interest.
%         spkMUAs - a sparse matrix of all individual unit and MUA vectors.
%         muas - all corresponding unit and MUA IDs.
%         muaMap - a matrix with columns of unit and MUA IDs, types (1 for
%                  MUA and unit ID otherwise), probe channel indices, probe
%                  vertical and horizontal coordinates in um, ISI
%                  violations, and cluster isolation distance.
%         MUAsAll_units - a population firing rate vector constructed from
%                         pure units only.
%         spkUnits - a sparse matrix of individual unit vectors.
%         units - corresponding unit IDs.
%         unitMap - a matrix with columns of unit IDs, types, probe channel
%                   index, probe vertical and horizontal coordinates in um,
%                   ISI violations, and cluster isolation distance.


%% Determine relevant data columns
unitIDColumn = find(contains(cellstr(headerUnits),'unit_id'));
acronymsColumn = find(contains(cellstr(headerUnits),'ecephys_structure_acronym'));
isolationDistColumn = find(contains(cellstr(headerUnits),'isolation_distance'));
isiViolationsColumn = find(contains(cellstr(headerUnits),'isi_violations'));
%firingRateColumn = find(contains(cellstr(headerUnits),'firing_rate'));
vPosColumn = find(contains(cellstr(headerUnits),'probe_vertical_position'));
hPosColumn = find(contains(cellstr(headerUnits),'probe_horizontal_position'));
channelColumn = find(contains(cellstr(headerUnits),'channel_local_index'));
APCCFColumn = find(contains(cellstr(headerUnits),'anterior_posterior_ccf_coordinate'));
DVCCFColumn = find(contains(cellstr(headerUnits),'dorsal_ventral_ccf_coordinate'));
LRCCFColumn = find(contains(cellstr(headerUnits),'left_right_ccf_coordinate'));
unitIDColumn2 = find(contains(cellstr(headerSpikes),'unit_id'));
spikeTimeColumn = find(contains(cellstr(headerSpikes),'time_since_stimulus_presentation_onset'));


%% Reduce full unit data to relevant units
if strcmp(area, 'lTh')
  areas = {'lVB', 'lLGN', 'lTh'};
elseif strcmp(area, 'lCA')
  areas = {'lCA1', 'lCA3', 'lCA'};
elseif strcmp(area, 'lHp')
  areas = {'lDG', 'lCA1', 'lCA3', 'lCA', 'lHp'};
elseif strcmp(area, 'lVIS')
  areas = {'lV1', 'lV2'};
elseif strcmp(area, 'rTh')
  areas = {'rVB', 'rLGN', 'rTh'};
elseif strcmp(area, 'rCA')
  areas = {'rCA1', 'rCA3', 'rCA'};
elseif strcmp(area, 'rHp')
  areas = {'rDG', 'rCA1', 'rCA3', 'rCA', 'rHp'};
elseif strcmp(area, 'rVIS')
  areas = {'rV1', 'rV2'};
else
  areas = area;
end
if isempty(areas)
  muasOIInds = true(size(unitMetadata,1),1);
else
  muasOIInds = contains(deblank(unitMetadata(:,acronymsColumn)), areas); %#ok<*FNDSB>
end
muaMetadata = unitMetadata(muasOIInds,:);
muas = muaMetadata(:,unitIDColumn);
muaMap = [double(cell2mat(muas)) ones(size(muas)) double(cell2mat(muaMetadata(:,channelColumn))+1) ...
  double(cell2mat(muaMetadata(:,hPosColumn))) double(cell2mat(muaMetadata(:,vPosColumn)))...
  double(cell2mat(muaMetadata(:,isiViolationsColumn))) double(cell2mat(muaMetadata(:,isolationDistColumn)))...
  double(cell2mat(muaMetadata(:,APCCFColumn))) double(cell2mat(muaMetadata(:,DVCCFColumn))) double(cell2mat(muaMetadata(:,LRCCFColumn)))];
unitsOIInds = cell2mat(muaMetadata(:,isolationDistColumn)) >= cluDist & cell2mat(muaMetadata(:,isiViolationsColumn)) <= refractCont;
units = muas(unitsOIInds);
unitMap = muaMap(unitsOIInds,:);


%% Get the population firing rate
spikeTimeInds = ismember(spikeData(:,unitIDColumn2), cell2mat(muas));
muaSpikeTimes = spikeData(spikeTimeInds,spikeTimeColumn)';
muaSpikeInds = round(muaSpikeTimes*sr);
bins = 1:max(muaSpikeInds);
MUAsAll = hist(muaSpikeInds,bins); %#ok<*HIST>


%% Get spike matrices
spkMUAs = zeros(numel(muas),numel(MUAsAll));
spkUnits = zeros(numel(units),numel(MUAsAll));
for iUnit = 1:numel(muas)
  spikeTimeInds = ismember(spikeData(:,unitIDColumn2), muas{iUnit});
  unitSpikeTimes = spikeData(spikeTimeInds,spikeTimeColumn)';
  unitSpikeInds = round(unitSpikeTimes*sr);
  spkMUAs(iUnit,:) = hist(unitSpikeInds,bins);
  if ismember(muas{iUnit}, cell2mat(units))
    spkUnits(cell2mat(units) == muas{iUnit},:) = spkMUAs(iUnit,:);
  end
end


%% Remove units with a low spike count
unitCounts = sum(spkUnits,2);
lowCountUnits = unitCounts < count;
spkUnits = spkUnits(~lowCountUnits,:);
units = units(~lowCountUnits);
unitMap = unitMap(~lowCountUnits,:);
muaMap(ismember(cell2mat(muas),cell2mat(units)), 2) = cell2mat(units);


%% Get the population firing rate constructed from pure units only
MUAsAll_units = sum(spkUnits,1);
assert(sum(logical((MUAsAll-MUAsAll_units) >= 0)) == numel(MUAsAll) && sum(logical((MUAsAll-MUAsAll_units) >= 0)) == numel(MUAsAll_units));


%% Compress output variables
spkMUAs = sparse(spkMUAs);
spkUnits = sparse(spkUnits);
muas = double(cell2mat(muas));
units = double(cell2mat(units));
muaMap = double(muaMap);
unitMap = double(unitMap);