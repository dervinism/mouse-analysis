function [seriesName1, seriesName2] = seriesNames(seriesName)
% A helper function that returns comparison series names

strSep = strfind(seriesName,'__');
strSep2 = strfind(seriesName,'s');
seriesName1 = seriesName(strSep2(1)+1:strSep-1);
seriesName2 = seriesName(strSep2(2)+1:end);