function [seriesName, animal] = seriesFromEntry(entryName)
% Given the entry name, determine the series name

strSep = strfind(entryName, 's');
seriesName = entryName(strSep+1:end);
animal = entryName(1:strSep-2);
end