function [areaCode, areaID, areaName, areaMatEntry] = determineAreaOld(area)
% [areaCode, areaID, areaName, areaMatEntry] = determineArea(area)
%
% Determine the area given its property.
% Input: area - an area property. This could be any of the following:
%                     numeric series name, area code, area ID, or areaName.
% OutputL: areaCode - a number that uniquely identifies the area.
%          areaID - area ID (seriesName(15:end).
%          areaName - area acronym (e.g., S1).
%          areaMatEntry - coherence matrix entry.

areaCode = []; %#ok<*NASGU>
areaID = [];
areaName = [];
areaMatEntry = [];

if ischar(area) && numel(area) > 4
  [areaCode, areaID, areaName, areaMatEntry] = determineAreaFromSeriesOld(area);
elseif (ischar(area) && (strcmpi(area, 'S1') || strcmpi(area, '1'))) || (isnumeric(area) && area == 1)
  areaCode = 1; areaID = '1'; areaName = 'S1'; areaMatEntry = 5;
elseif (ischar(area) && (strcmpi(area, 'VB1') || strcmpi(area, '2'))) || (isnumeric(area) && area == 2)
  areaCode = 2; areaID = '2'; areaName = 'VB1'; areaMatEntry = 1;
elseif (ischar(area) && (strcmpi(area, 'Th1') || strcmpi(area, '24'))) || (isnumeric(area) && area == 3)
  areaCode = 3; areaID = '24'; areaName = 'Th1'; areaMatEntry = 3;
elseif (ischar(area) && (strcmpi(area, 'DG') || strcmpi(area, '5'))) || (isnumeric(area) && area == 4)
  areaCode = 4; areaID = '5'; areaName = 'DG'; areaMatEntry = 7;
elseif (ischar(area) && (strcmpi(area, 'CA1') || strcmpi(area, '6'))) || (isnumeric(area) && area == 5)
  areaCode = 5; areaID = '6'; areaName = 'CA1'; areaMatEntry = 9;
elseif (ischar(area) && (strcmpi(area, 'Hp') || strcmpi(area, '56'))) || (isnumeric(area) && area == 6)
  areaCode = 6; areaID = '56'; areaName = 'Hp'; areaMatEntry = 10;
elseif (ischar(area) && (strcmpi(area, 'RSC') || strcmpi(area, '7'))) || (isnumeric(area) && area == 7)
  areaCode = 7; areaID = '7'; areaName = 'RSC'; areaMatEntry = 6;
elseif (ischar(area) && (strcmpi(area, 'VB2') || strcmpi(area, '8'))) || (isnumeric(area) && area == 8)
  areaCode = 8; areaID = '8'; areaName = 'VB2'; areaMatEntry = 2;
elseif (ischar(area) && (strcmpi(area, 'Th2') || strcmpi(area, '810'))) || (isnumeric(area) && area == 9)
  areaCode = 9; areaID = '810'; areaName = 'Th2'; areaMatEntry = 4;
elseif (ischar(area) && (strcmpi(area, 'CA3') || strcmpi(area, '11'))) || (isnumeric(area) && area == 10)
  areaCode = 10; areaID = '11'; areaName = 'CA3'; areaMatEntry = 8;
elseif (ischar(area) && (strcmpi(area, 'mPFC') || strcmpi(area, '12'))) || (isnumeric(area) && area == 11)
  areaCode = 11; areaID = '12'; areaName = 'mPFC';
elseif (ischar(area) && (strcmpi(area, 'lV1') || strcmpi(area, '13'))) || (isnumeric(area) && area == 12)
  areaCode = 12; areaID = '13'; areaName = 'lV1';
elseif (ischar(area) && (strcmpi(area, 'rV1') || strcmpi(area, '14'))) || (isnumeric(area) && area == 13)
  areaCode = 13; areaID = '14'; areaName = 'rV1';
elseif (ischar(area) && (strcmpi(area, 'Po') || strcmpi(area, '3'))) || (isnumeric(area) && area == 14)
  areaCode = 14; areaID = '3'; areaName = 'Po';
elseif (ischar(area) && (strcmpi(area, 'LP1') || strcmpi(area, '4'))) || (isnumeric(area) && area == 15)
  areaCode = 15; areaID = '4'; areaName = 'LP1';
elseif (ischar(area) && (strcmpi(area, 'LP2') || strcmpi(area, '9'))) || (isnumeric(area) && area == 16)
  areaCode = 16; areaID = '9'; areaName = 'LP2';
elseif (ischar(area) && (strcmpi(area, 'DLGN') || strcmpi(area, '10'))) || (isnumeric(area) && area == 17)
  areaCode = 17; areaID = '10'; areaName = 'DLGN';
elseif (ischar(area) && (strcmpi(area, 'VB') || strcmpi(area, '15'))) || (isnumeric(area) && area == 18)
  areaCode = 18; areaID = '15'; areaName = 'VB'; areaMatEntry = 11;
elseif (ischar(area) && (strcmpi(area, 'Th') || strcmpi(area, '16'))) || (isnumeric(area) && area == 19)
  areaCode = 19; areaID = '16'; areaName = 'Th'; areaMatEntry = 12;
elseif (ischar(area) && (strcmpi(area, 'CA') || strcmpi(area, '17'))) || (isnumeric(area) && area == 20)
  areaCode = 20; areaID = '17'; areaName = 'CA'; areaMatEntry = 13;
elseif (ischar(area) && (strcmpi(area, 'Cx') || strcmpi(area, '18'))) || (isnumeric(area) && area == 21)
  areaCode = 21; areaID = '18'; areaName = 'Cx';
elseif ~ischar(area) && ~isnumeric(area)
  error('Function determineArea only recognises string or numeric input values.')
else
  error('Unknown area.')
end