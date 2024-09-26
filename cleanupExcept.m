function cleanupExcept(varargin)
% Clean-up before running scripts except specified variables

if isempty(varargin)
  clear
else
  for iVar = 1:numel(varargin)
    strEval = ['clearvars -except varargin ' varargin{iVar}];
    eval(strEval);
  end
  clear varargin
end

fclose all;
close all
clc