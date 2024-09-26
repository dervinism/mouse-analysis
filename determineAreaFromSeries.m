function [areaCode, areaID, areaName, areaMatEntry, probeID, areaCodeAugmented] = determineAreaFromSeries(seriesName)
% [areaCode, areaID, areaName, areaMatEntry, probeID, areaCodeAugmented] = determineAreaFromSeries(seriesName)
%
% Determine the area the series belongs to.
% Input: seriesName - a series name in the form of long numeric identifier
%                     (14 digits long or more; in some rare instances can
%                     be shorter).
% OutputL: areaCode - a number that uniquely identifies the area.
%          areaID - area ID (seriesName(15:end).
%          areaName - area acronym (e.g., S1).
%          areaMatEntry - coherence matrix entry.
%          probeID - 1 (S1) or 2 (RSC) if it was recorded at UOL.
%          areaCodeAugmented - augmented area number.

areaCode = []; %#ok<*NASGU>
areaID = [];
areaName = [];
areaMatEntry = [];
areaCodeAugmented = [];

% Go through series names and determine areaID
if numel(seriesName) > 14
  if strcmpi(seriesName, '201811281820101') || strcmpi(seriesName, '201811282042541')
    areaID = '13'; % left V1
  elseif strcmpi(seriesName, '201811281820102') || strcmpi(seriesName, '201811282042542')
    areaID = '14'; % right V1
  else
    areaID = seriesName(15:end);
  end
elseif numel(seriesName) < 14
  if strcmpi(seriesName, '20180712') || strcmpi(seriesName, '20180713') || strcmpi(seriesName, '20180719')
    areaID = '13'; % left V1
  elseif strcmpi(seriesName, '20180711') || strcmpi(seriesName, '20180717')
    areaID = '14'; % right V1
  else
    error('Unknown area');
  end
elseif numel(seriesName) == 14
  if strcmpi(seriesName, '20190326130324') || strcmpi(seriesName, '20190328102204') || strcmpi(seriesName, '20190328105028') ||...
      strcmpi(seriesName, '20190328114041') || strcmpi(seriesName, '20190511141202') || strcmpi(seriesName, '20190514115945') ||...
      strcmpi(seriesName, '20190514140920') || strcmpi(seriesName, '20190515122749') || strcmpi(seriesName, '20190515124626') ||...
      strcmpi(seriesName, '20190515135343') || strcmpi(seriesName, '20190627160609') || strcmpi(seriesName, '20190627161713') ||...
      strcmpi(seriesName, '20190627170248') || strcmpi(seriesName, '20190708155832') || strcmpi(seriesName, '20190708161317') ||...
      strcmpi(seriesName, '20190708171625')
    areaID = '1'; % left S1
  elseif strcmpi(seriesName, '20181219182543') || strcmpi(seriesName, '20181221145727') || strcmpi(seriesName, '20190102134858') ||...
      strcmpi(seriesName, '20190106174515') || strcmpi(seriesName, '20181221165147') || strcmpi(seriesName, '20190102162201') ||...
      strcmpi(seriesName, '20190106200423')
    areaID = '28'; % left VB
  elseif strcmpi(seriesName, '20180313172457') || strcmpi(seriesName, '20181005134915') ||...
      strcmpi(seriesName, '20181008164836') || strcmpi(seriesName, '20181019140918') || strcmpi(seriesName, '20181008184928')
    areaID = '12'; % left mPFC
  elseif strcmpi(seriesName, '20181128182010') || strcmpi(seriesName, '20181128204254') ||...
      strcmpi(seriesName, '20181129171747') || strcmpi(seriesName, '20181129195448')
    areaID = '13'; % left V1
  else
    error('Unknown area');
  end
end

% Determine areaCode and areaName
[areaCode, areaID, areaName, areaMatEntry, probeID, areaCodeAugmented] = determineArea(areaID);
end