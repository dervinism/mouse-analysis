function [areaName1, areaName2] = comparison2areas(comparison)
% [areaName1, areaName2] = comparison2areas(comparison)
%
% Function returns area names given a comparison string (e.g., 'VBVsS1').

iStr = strfind(comparison, 'Vs');
areaName1 = comparison(1:iStr-1);
areaName2 = comparison(iStr+2:end);