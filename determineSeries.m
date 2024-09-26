function [series, seriesV1, nSeries, selectedUnitMetadata] = determineSeries(headerUnits, unitMetadata)
% [series, seriesV1, nSeries, selectedUnitMetadata] = determineSeries(headerUnits, unitMetadata)
%
% Function determines brain areas corresponding to recording series, which
% series contain V1 data, a number of recording series, and taylors the
% unit metadata specifically to the recording series of interest given the
% unit metadata information derived from allensdk repository.
% Input: headerUnits - unit metadata header.
%        unitMetadata.
% Output: series - areas corresponding to recording series of interest in
%                  the given unit metadata.
%         seriesV1 - a logical vector indicating series containing V1 data.
%         nSeries - a number of existing recording series of interest in
%                   the given unit metadata.
%         selectedUnitMetadata - taylored unit metadata.
%
% We are interested in series with data in these areas: lVB, lLGN, lLP, lTh,
% lDG, lCA1, lCA3, lCA, lHp, lV1, lV2, lVIS, rVB, rLGN, rLP, rTh, rDG, rCA1,
% rCA3, rCA, rHp, rV1, rV2, rVIS.
% Data exists on the following thalamic nuclei:
%   LGd - lateral geniculate complex, dorsal part;
%   LD - lateral dorsal nucleus;
%   LP - lateral posterior nucleus
%   VPM - ventral posteromedial nucleus
%   MGm - medial geniculate complex, medial part
%   MGv - medial geniculate complex, ventral part
%   MGd - medial geniculate complex, dorsal part
%   PO - posterior complex
%   LGv -lateral geniculate complex, ventral part;
%   VL - ventral lateral nucleus
%   VPL - ventral posterior lateral nucleus
%   POL - posterior limiting nucleus
%   Eth - ethmoid nucleus
%   PoT - posterior triangular nucleus
%   PP - peripeduncular nucleus
%   PIL - posterior intralaminar nucleus
%   IntG - intermediate geniculate nucleus
%   IGL - intergeniculate leaflet of the lateral geniculate complex
%   SGN - suprageniculate nucleus
%   PF - parafascicular nucleus

leftRightSplitCoord = 5700; % left-right CCF coordinate dividing left and right brain hemispheres

acronymsColumn = find(contains(cellstr(headerUnits),'ecephys_structure_acronym'));
leftRightCCFColumn = find(contains(cellstr(headerUnits),'left_right_ccf_coordinate'));

% Brain structure acronyms corresponding to broader brain areas
VB  = 'VPM VPL';
LGN = 'LGd LGv';
LP = 'LP';
Th  = 'LGd LD LP VPM TH MGm MGv MGd PO LGv VL VPL POL Eth PoT PP PIL IntG IGL SGN PF';
DG  = 'DG';
CA1 = 'CA1';
CA3 = 'CA2 CA3';
CA = 'CA1 CA2 CA3'; %#ok<*NASGU>
Hp = 'DG CA1 CA2 CA3';
V1  = 'VISp';
V2  = 'VISl VISrl VISam VISpm VISal VISmma VISmmp VISli VIS';
VIS  = 'VISp VISl VISrl VISam VISpm VISal VISmma VISmmp VISli VIS'; % All visual cortical areas

seriesOI = {'lVB', 'lLGN', 'lLP', 'lTh', 'lDG', 'lCA1', 'lCA3', 'lCA', 'lHp',...
  'lV1', 'lV2', 'lVIS', 'rVB', 'rLGN', 'rLP', 'rTh', 'rDG', 'rCA1', 'rCA3', 'rCA', 'rHp',...
  'rV1', 'rV2', 'rVIS'};

% Identify all existing areas of interest
areas = cell(size(unitMetadata,1),1);
for iUnit = 1:size(unitMetadata,1)
  if contains(VB, deblank(unitMetadata(iUnit,acronymsColumn)))
    if unitMetadata{iUnit,leftRightCCFColumn} < leftRightSplitCoord
      areas{iUnit} = 'lVB';
    else
      areas{iUnit} = 'rVB';
    end
  elseif contains(LGN, deblank(unitMetadata(iUnit,acronymsColumn)))
    if unitMetadata{iUnit,leftRightCCFColumn} < leftRightSplitCoord
      areas{iUnit} = 'lLGN';
    else
      areas{iUnit} = 'rLGN';
    end
  elseif contains(LP, deblank(unitMetadata(iUnit,acronymsColumn)))
    if unitMetadata{iUnit,leftRightCCFColumn} < leftRightSplitCoord
      areas{iUnit} = 'lLP';
    else
      areas{iUnit} = 'rLP';
    end
  elseif contains(Th, deblank(unitMetadata(iUnit,acronymsColumn)))
    if unitMetadata{iUnit,leftRightCCFColumn} < leftRightSplitCoord
      areas{iUnit} = 'lTh';
    else
      areas{iUnit} = 'rTh';
    end
  elseif contains(DG, deblank(unitMetadata(iUnit,acronymsColumn)))
    if unitMetadata{iUnit,leftRightCCFColumn} < leftRightSplitCoord
      areas{iUnit} = 'lDG';
    else
      areas{iUnit} = 'rDG';
    end
  elseif contains(CA1, deblank(unitMetadata(iUnit,acronymsColumn)))
    if unitMetadata{iUnit,leftRightCCFColumn} < leftRightSplitCoord
      areas{iUnit} = 'lCA1';
    else
      areas{iUnit} = 'rCA1';
    end
  elseif contains(CA3, deblank(unitMetadata(iUnit,acronymsColumn)))
    if unitMetadata{iUnit,leftRightCCFColumn} < leftRightSplitCoord
      areas{iUnit} = 'lCA3';
    else
      areas{iUnit} = 'rCA3';
    end
  elseif contains(V1, deblank(unitMetadata(iUnit,acronymsColumn)))
    if unitMetadata{iUnit,leftRightCCFColumn} < leftRightSplitCoord
      areas{iUnit} = 'lV1';
    else
      areas{iUnit} = 'rV1';
    end
  elseif contains(V2, deblank(unitMetadata(iUnit,acronymsColumn)))
    if unitMetadata{iUnit,leftRightCCFColumn} < leftRightSplitCoord
      areas{iUnit} = 'lV2';
    else
      areas{iUnit} = 'rV2';
    end
  else
    areas{iUnit} = 'NaN';
  end
end

series = unique(areas);
series = series(~contains(series, 'NaN'));
if sum(contains(series, {'lVB', 'lLP', 'lLGN'})) && ~sum(contains(series, {'lTh'}))
  series{end+1} = 'lTh';
end
if sum(contains(series, {'rVB', 'rLP', 'rLGN'})) && ~sum(contains(series, {'rTh'}))
  series{end+1} = 'rTh';
end
if sum(contains(series, {'lCA1', 'lCA2', 'lCA3'}))
  series{end+1} = 'lCA';
end
if sum(contains(series, {'rCA1', 'rCA2', 'rCA3'}))
  series{end+1} = 'rCA';
end
if sum(contains(series, {'lDG', 'lCA1', 'lCA2', 'lCA3'})) && ~sum(contains(series, {'lHp'}))
  series{end+1} = 'lHp';
end
if sum(contains(series, {'rDG', 'rCA1', 'rCA2', 'rCA3'})) && ~sum(contains(series, {'rHp'}))
  series{end+1} = 'rHp';
end
if sum(contains(series, {'lV1', 'lV2'})) && ~sum(contains(series, {'lVIS'}))
  series{end+1} = 'lVIS';
end
if sum(contains(series, {'rV1', 'rV2'})) && ~sum(contains(series, {'rVIS'}))
  series{end+1} = 'rVIS';
end
nSeries = numel(series);

% Arrange series
seriesInds = endsWith(seriesOI, series);
series = seriesOI(seriesInds);
seriesV1 = contains(series, {'lV1', 'rV1'});

% Taylor unit metadata to specific series of interest
selectedUnitMetadata = unitMetadata(~contains(areas, 'NaN'),:);
selectedUnitMetadata(:,acronymsColumn) = areas(~contains(areas, 'NaN'));