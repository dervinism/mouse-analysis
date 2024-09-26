function [areaCode, areaID, areaName, areaMatEntry, probeID, areaCodeAugmented] = determineArea(area)
% [areaCode, areaID, areaName, areaMatEntry, probeID, areaCodeAugmented] = determineArea(area)
%
% Determine the area given its property.
% Input: area - an area property. This could be any of the following:
%                 numeric series name, area code, area ID, or areaName.
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
probeID = [];
areaCodeAugmented = [];

if ischar(area) && numel(area) > 5
  [areaCode, areaID, areaName, areaMatEntry, probeID, areaCodeAugmented] = determineAreaFromSeries(area);
elseif (ischar(area) && (strcmpi(area, 'lS1') || strcmpi(area, '1'))) || (isnumeric(area) && area == 1)
  areaCode = [1 42 46]; areaID = '1'; areaName = 'lS1'; areaMatEntry = 4; probeID = 1; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'lVB1') || strcmpi(area, '2'))) || (isnumeric(area) && area == 2)
  areaCode = 2; areaID = '2'; areaName = 'lVB1'; probeID = 2; areaCodeAugmented = [areaCode 33];
elseif (ischar(area) && (strcmpi(area, 'lPo') || strcmpi(area, '3'))) || (isnumeric(area) && area == 3)
  areaCode = [3 52]; areaID = '3'; areaName = 'lPo'; probeID = 2; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'lLP1') || strcmpi(area, '4'))) || (isnumeric(area) && area == 4)
  areaCode = 4; areaID = '4'; areaName = 'lLP1'; probeID = 2; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'lTh1') || strcmpi(area, '24'))) || (isnumeric(area) && area == 5)
  areaCode = 5; areaID = '24'; areaName = 'lTh1'; probeID = 2; areaCodeAugmented = [areaCode 35];
elseif (ischar(area) && (strcmpi(area, 'lDG') || strcmpi(area, '5'))) || (isnumeric(area) && area == 6)
  areaCode = [6 36]; areaID = '5'; areaName = 'lDG'; probeID = 2; areaMatEntry = 7; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'lCA1') || strcmpi(area, '6'))) || (isnumeric(area) && area == 7)
  areaCode = [7 37]; areaID = '6'; areaName = 'lCA1'; probeID = 2; areaMatEntry = 8; areaCodeAugmented = [areaCode 39];
elseif (ischar(area) && (strcmpi(area, 'lHp1') || strcmpi(area, '56'))) || (isnumeric(area) && area == 8)
  areaCode = 8; areaID = '56'; areaName = 'lHp1'; probeID = 2; areaCodeAugmented = [areaCode 40];
elseif (ischar(area) && (strcmpi(area, 'lRSC') || strcmpi(area, '7'))) || (isnumeric(area) && area == 9)
  areaCode = [9 42 47]; areaID = '7'; areaName = 'lRSC'; probeID = 2; areaMatEntry = 6; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'lVB2') || strcmpi(area, '8'))) || (isnumeric(area) && area == 10)
  areaCode = 10; areaID = '8'; areaName = 'lVB2'; probeID = 1; areaCodeAugmented = [areaCode 33];
elseif (ischar(area) && (strcmpi(area, 'lLP2') || strcmpi(area, '9'))) || (isnumeric(area) && area == 11)
  areaCode = 11; areaID = '9'; areaName = 'lLP2'; probeID = 1; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'lLGN') || strcmpi(area, '10'))) || (isnumeric(area) && area == 12)
  areaCode = [12 34]; areaID = '10'; areaName = 'lLGN'; areaMatEntry = 2; probeID = 1; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'lTh2') || strcmpi(area, '810'))) || (isnumeric(area) && area == 13)
  areaCode = 13; areaID = '810'; areaName = 'lTh2'; probeID = 1; areaCodeAugmented = [areaCode 35];
elseif (ischar(area) && (strcmpi(area, 'lCA3') || strcmpi(area, '11'))) || (isnumeric(area) && area == 14)
  areaCode = [14 38]; areaID = '11'; areaName = 'lCA3'; probeID = 1; areaMatEntry = 9; areaCodeAugmented = [areaCode 39];
elseif (ischar(area) && (strcmpi(area, 'lVB') || strcmpi(area, '28'))) || (isnumeric(area) && area == 15)
  areaCode = [15 33]; areaID = '28'; areaName = 'lVB'; areaMatEntry = 1; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'lLP') || strcmpi(area, '49'))) || (isnumeric(area) && area == 16)
  areaCode = [16 48]; areaID = '49'; areaName = 'lLP'; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'lTh') || strcmpi(area, '210'))) || (isnumeric(area) && area == 17)
  areaCode = [17 35]; areaID = '210'; areaName = 'lTh'; areaMatEntry = 3; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'lCA') || strcmpi(area, '611'))) || (isnumeric(area) && area == 18)
  areaCode = [18 39]; areaID = '611'; areaName = 'lCA'; areaMatEntry = 10; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'lHp') || strcmpi(area, '511'))) || (isnumeric(area) && area == 19)
  areaCode = [19 40]; areaID = '511'; areaName = 'lHp'; areaMatEntry = 11; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'lmPFC') || strcmpi(area, '12'))) || (isnumeric(area) && area == 20)
  areaCode = 20; areaID = '12'; areaName = 'lmPFC'; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'lV1') || strcmpi(area, '13'))) || (isnumeric(area) && area == 21)
  areaCode = [21 41]; areaID = '13'; areaName = 'lV1'; areaMatEntry = 5; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'rV1') || strcmpi(area, '14'))) || (isnumeric(area) && area == 22)
  areaCode = [22 41]; areaID = '14'; areaName = 'rV1'; areaMatEntry = 5; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'rVB') || strcmpi(area, '15'))) || (isnumeric(area) && area == 23)
  areaCode = [23 33]; areaID = '15'; areaName = 'rVB'; areaMatEntry = 1; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'rPo') || strcmpi(area, '16'))) || (isnumeric(area) && area == 24)
  areaCode = [24 52]; areaID = '16'; areaName = 'rPo'; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'rLP') || strcmpi(area, '17'))) || (isnumeric(area) && area == 25)
  areaCode = [25 48]; areaID = '17'; areaName = 'rLP'; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'rLGN') || strcmpi(area, '18'))) || (isnumeric(area) && area == 26)
  areaCode = [26 34]; areaID = '18'; areaName = 'rLGN'; areaMatEntry = 2; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'rTh') || strcmpi(area, '1518'))) || (isnumeric(area) && area == 27)
  areaCode = [27 35]; areaID = '1518'; areaName = 'rTh'; areaMatEntry = 3; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'rDG') || strcmpi(area, '19'))) || (isnumeric(area) && area == 28)
  areaCode = [28 36]; areaID = '19'; areaName = 'rDG'; areaMatEntry = 7; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'rCA1') || strcmpi(area, '20'))) || (isnumeric(area) && area == 29)
  areaCode = [29 37]; areaID = '20'; areaName = 'rCA1'; areaMatEntry = 8; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'rCA3') || strcmpi(area, '21'))) || (isnumeric(area) && area == 30)
  areaCode = [30 38]; areaID = '21'; areaName = 'rCA3'; areaMatEntry = 9; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'rCA') || strcmpi(area, '2021'))) || (isnumeric(area) && area == 31)
  areaCode = [31 39]; areaID = '2021'; areaName = 'rCA'; areaMatEntry = 10; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'rHp') || strcmpi(area, '1921'))) || (isnumeric(area) && area == 32)
  areaCode = [32 40]; areaID = '1921'; areaName = 'rHp'; areaMatEntry = 11; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'VB') || strcmpi(area, '22'))) || (isnumeric(area) && area == 33)
  areaCode = 33; areaID = '22'; areaName = 'VB'; areaCodeAugmented = areaCode; %areaMatEntry = 1;
elseif (ischar(area) && (strcmpi(area, 'LGN') || strcmpi(area, '23'))) || (isnumeric(area) && area == 34)
  areaCode = 34; areaID = '23'; areaName = 'LGN'; areaCodeAugmented = areaCode; %areaMatEntry = 2;
elseif (ischar(area) && (strcmpi(area, 'Th') || strcmpi(area, '24'))) || (isnumeric(area) && area == 35)
  areaCode = 35; areaID = '24'; areaName = 'Th'; areaCodeAugmented = areaCode; %areaMatEntry = 3;
elseif (ischar(area) && (strcmpi(area, 'DG') || strcmpi(area, '25'))) || (isnumeric(area) && area == 36)
  areaCode = 36; areaID = '25'; areaName = 'DG'; areaCodeAugmented = areaCode; %areaMatEntry = 7;
elseif (ischar(area) && (strcmpi(area, 'CA1') || strcmpi(area, '26'))) || (isnumeric(area) && area == 37)
  areaCode = 37; areaID = '26'; areaName = 'CA1'; areaCodeAugmented = areaCode; %areaMatEntry = 8;
elseif (ischar(area) && (strcmpi(area, 'CA3') || strcmpi(area, '27'))) || (isnumeric(area) && area == 38)
  areaCode = 38; areaID = '27'; areaName = 'CA3'; areaCodeAugmented = areaCode; %areaMatEntry = 9;
elseif (ischar(area) && (strcmpi(area, 'CA') || strcmpi(area, '28'))) || (isnumeric(area) && area == 39)
  areaCode = 39; areaID = '28'; areaName = 'CA'; areaCodeAugmented = areaCode; %areaMatEntry = 10;
elseif (ischar(area) && (strcmpi(area, 'Hp') || strcmpi(area, '29'))) || (isnumeric(area) && area == 40)
  areaCode = 40; areaID = '29'; areaName = 'Hp'; areaCodeAugmented = areaCode; %areaMatEntry = 11;
elseif (ischar(area) && (strcmpi(area, 'V1') || strcmpi(area, '1314'))) || (isnumeric(area) && area == 41)
  areaCode = 41; areaID = '1314'; areaName = 'V1'; areaCodeAugmented = areaCode; %areaMatEntry = 5;
elseif (ischar(area) && (strcmpi(area, 'Cx') || strcmpi(area, '30'))) || (isnumeric(area) && area == 42)
  areaCode = 42; areaID = '30'; areaName = 'Cx'; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'lVIS') || strcmpi(area, '31'))) || (isnumeric(area) && area == 43)
  areaCode = [43 45]; areaID = '31'; areaName = 'lVIS'; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'rVIS') || strcmpi(area, '32'))) || (isnumeric(area) && area == 44)
  areaCode = [44 45]; areaID = '32'; areaName = 'rVIS'; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'VIS') || strcmpi(area, '33'))) || (isnumeric(area) && area == 45)
  areaCode = 45; areaID = '33'; areaName = 'VIS'; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'S1') || strcmpi(area, '34'))) || (isnumeric(area) && area == 46)
  areaCode = 46; areaID = '34'; areaName = 'S1'; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'RSC') || strcmpi(area, '35'))) || (isnumeric(area) && area == 47)
  areaCode = 47; areaID = '35'; areaName = 'RSC'; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'LP') || strcmpi(area, '36'))) || (isnumeric(area) && area == 48)
  areaCode = 48; areaID = '36'; areaName = 'LP'; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'lV2') || strcmpi(area, '37'))) || (isnumeric(area) && area == 49)
  areaCode = [49 51]; areaID = '37'; areaName = 'lV2'; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'rV2') || strcmpi(area, '38'))) || (isnumeric(area) && area == 50)
  areaCode = [50 51]; areaID = '38'; areaName = 'rV2'; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'V2') || strcmpi(area, '39'))) || (isnumeric(area) && area == 51)
  areaCode = 51; areaID = '39'; areaName = 'V2'; areaCodeAugmented = areaCode;
elseif (ischar(area) && (strcmpi(area, 'Po') || strcmpi(area, '40'))) || (isnumeric(area) && area == 52)
  areaCode = 52; areaID = '40'; areaName = 'Po'; areaCodeAugmented = areaCode;
elseif ~ischar(area) && ~isnumeric(area)
  error('Function determineArea only recognises string or numeric input values.')
else
  error('Unknown area.')
end