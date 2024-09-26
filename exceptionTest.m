function [breakClause, iExcept] = exceptionTest(exceptionsList, seriesName1, seriesName2)
% [breakClause, iExcept] = exceptionTest(exceptionsList, seriesName1, seriesName2)
%
% Function tests if given series are exceptions that should be skipped
% over.
% Input: exceptionsList - a cell array with a list of series that are
%                         exceptions.
%        seriesName1 - the name of the series to be tested.
%        seriesName2 - the name of the second series to be tested
%                      (optional).
% Output: breakClause - true if any of the given series names are
%                       excpetions, false otherwise.
%         iExcept - a list index where the exception occured.

if nargin < 3
  seriesName2 = [];
end

for iExcept = 1:numel(exceptionsList)
  if strcmpi(seriesName1, exceptionsList{iExcept})
    breakClause = true;
    break
  elseif ~isempty(seriesName2) && strcmpi(seriesName2, exceptionsList{iExcept})
    breakClause = true;
    break
  else
    breakClause = false;
  end
end