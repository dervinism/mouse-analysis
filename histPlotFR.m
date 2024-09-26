function [h, stats, options] = histPlotFR(edges, data, options)
% [h, stats] = histPlotFR(edges, data, options)
%
% Function histPlot draws a histogram.
% Input: edges
%        data - log firing rates data column matrix or a cell array.
%        options - a structure variable with fields:
%          figName - a figure name and filename string.
%          figFolder
%          xLabel
%          xLim
%          yLabel
%          yLim
%          legendLabels
%          figSize - a figure window size.
%          saveFig - a logical for saving the figure. Default is true.
%          plotType - a string variable to switch between two figure types.
%                     If set to 'extremes', would grey out intervening
%                     histograms between the first and the last histograms.
%                     Default is 'regular'.
%          colours - a row matrix of RGB colour codes corresponding to data
%                    columns. If left empty, matlabColours function will
%                    generate colour codes.
%          lineStyles - a cell array of line styles (e.g., '--', ':').
%          markerStyles - a cell array of marker styles (e.g., 'o', '^').
%          stats - if this field is left non-empty, stats will be displayed
%                  comparing variances of two data columns with the largest
%                  difference in their variances. The fields of stats
%                  structure variable are described in the varTest function.
%          binnedData - a column matrix with the original data already
%                       binned. This is useful when data needs to be
%                       scaled. By default the matrix is left empty.
%          markerVPos - marker vertical position. Default is 5% higher than
%                       the largest bin count.
%          fH - figure handle. If supplied, draws in the supplied figure
%               window.
%          lineWidths - a vector with line widths. Default is 2.
%          displayMeans - true (default) or false for marking means.
%          displayVars - true or false (default) for marking variances.
%          pdf - true or false (default) for drawing a probability density
%                function instead of histogram counts.
%          
% Output: h - figure handle.
%         stats - paired sample t-test statistics.


%% Parse user input
if isempty(edges)
  error('The first input variable is empty.');
end
h = [];
stats = [];
if isempty(data)
  return
end
if nargin < 3 || ~isfield(options, 'figName') || isempty(options.figName)
  options.figName = 'countHistogram';
end
if nargin < 3 || ~isfield(options, 'figFolder') || isempty(options.figFolder)
  options.figFolder = pwd;
end
if nargin < 3 || ~isfield(options, 'xLabel') || isempty(options.xLabel)
  options.xLabel = 'X';
end
if nargin < 3 || ~isfield(options, 'xLim') || isempty(options.xLim)
  options.xLim = [];
end
if nargin < 3 || ~isfield(options, 'yLabel') || isempty(options.yLabel)
  options.yLabel = 'Count';
end
if nargin < 3 || ~isfield(options, 'yLim') || isempty(options.yLim)
  options.yLim = [];
end
if nargin < 3 || ~isfield(options, 'legendLabels') || isempty(options.legendLabels)
  options.legendLabels = [];
end
if nargin < 3 || ~isfield(options, 'figSize') || isempty(options.figSize)
  options.figSize = 15;
end
if nargin < 3 || ~isfield(options, 'saveFig') || isempty(options.saveFig)
  options.saveFig = true;
end
if nargin < 3 || ~isfield(options, 'plotType') || isempty(options.plotType)
  options.plotType = 'regular';
end
if nargin < 3 || ~isfield(options, 'colours') || isempty(options.colours)
  options.colours = [];
end
if nargin < 3 || ~isfield(options, 'lineStyles') || isempty(options.lineStyles)
  options.lineStyles = [];
end
if nargin < 3 || ~isfield(options, 'markerStyles') || isempty(options.markerStyles)
  options.markerStyles = [];
end
if nargin < 3 || ~isfield(options, 'plotType') || isempty(options.plotType)
  options.plotType = 'regular';
end
if nargin < 3 || ~isfield(options, 'stats') || isempty(options.stats)
  options.stats = [];
end
if nargin < 3 || ~isfield(options, 'binnedData') || isempty(options.binnedData)
  options.binnedData = [];
end
if nargin < 3 || ~isfield(options, 'markerVPos') || isempty(options.markerVPos)
  options.markerVPos = [];
end
if nargin < 3 || ~isfield(options, 'fH') || isempty(options.fH)
  options.fH = [];
end
if nargin < 3 || ~isfield(options, 'lineWidths') || isempty(options.lineWidths)
  options.lineWidths = 2;
end
if nargin < 3 || ~isfield(options, 'displayMeans') || isempty(options.displayMeans)
  options.displayMeans = true;
end
if nargin < 3 || ~isfield(options, 'displayVars') || isempty(options.displayVars)
  options.displayVars = false;
end
if nargin < 3 || ~isfield(options, 'pdf') || isempty(options.pdf)
  options.pdf = false;
end


%% Estimate data mean and transform the data into a histogram
if iscell(data)
  nCol = numel(data);
  dataMean = zeros(1,nCol);
  dataCI95 = zeros(2,nCol);
  bars = zeros(nCol, numel(edges)-1);
  for col = 1:nCol
    try
      [dataMean(col), dataCI95(:,col)] = datamean(data{col});
    catch
      dataMean(col) = datamean(data{col});
    end
    if isempty(options.binnedData)
      bars(col,:) = histcounts(data{col}, edges);
    else
      bars(col,:) = options.binnedData(col,:);
    end
  end
else
  nCol = size(data,2);
  [dataMean, dataCI95] = datamean(data);
  if isempty(options.binnedData)
    bars = zeros(size(data,2), numel(edges)-1);
    for col = 1:nCol
      bars(col,:) = histcounts(data(:,col), edges);
    end
  else
    bars = options.binnedData;
  end
end

if options.pdf
  barsPDF = zeros(size(bars));
  binSize = edges(2)-edges(1);
  nRows2sum = round(size(bars,1)/2);
  for iRow = 1:size(bars,1)
    if iRow <= nRows2sum
      barsPDF(iRow,:) = bars(iRow,:)./(sum(sum(bars(1:nRows2sum,:),2))*binSize);
    else
      barsPDF(iRow,:) = bars(iRow,:)./(sum(sum(bars(nRows2sum+1:end,:),2))*binSize);
    end
  end
  bars = barsPDF;
end


%% Draw the histogram
if isempty(options.fH)
  h = figProperties(options.figName, 'normalized', [0, .005, .97, .90], 'w', 'on');
else
  h = figure(options.fH); hold on
end
if ~isempty(options.legendLabels)
  l = legend;
  nLegendEntries = numel(l.PlotChildren);
end

% Histos
centres = edges(2:end) - (edges(2)-edges(1))/2;
interveningLineColour = [200 200 200]./255;
p = [];
if numel(options.lineWidths) == 1
  options.lineWidths = repmat(options.lineWidths,1,nCol);
end
for col = 1:nCol
  if isempty(options.colours)
    lineColour = matlabColours(col);
  else
    lineColour = options.colours(col,:);
  end
  if isempty(options.lineStyles)
    lineStyle = '-';
  else
    lineStyle = options.lineStyles{col};
  end
  if strcmpi(options.plotType, 'extremes') && col > 1 && col < nCol
    p0 = plot(centres, bars(col,:), 'Color',interveningLineColour, 'LineWidth',options.lineWidths(col), 'LineStyle',lineStyle);
    uistack(p0,'bottom');
  else
    p0 = plot(centres, bars(col,:), 'Color',lineColour, 'LineWidth',options.lineWidths(col), 'LineStyle',lineStyle);
    p = [p p0]; %#ok<*AGROW>
  end
  if col == 1
    hold on
  end
end

% Mean markers
if options.displayMeans
  if isempty(options.markerVPos)
    markerVPos = 1.05*max(max(bars));
    options.markerVPos = markerVPos;
  else
    markerVPos = options.markerVPos;
  end
  for col = 1:nCol
    if isempty(options.colours)
      lineColour = matlabColours(col);
    else
      lineColour = options.colours(col,:);
    end
    if isempty(options.markerStyles)
      markerStyle = 'o';
    else
      markerStyle = options.markerStyles{col};
    end
    if strcmpi(options.plotType, 'extremes') && col > 1 && col < nCol
      if isempty(dataCI95)
        p0 = plot(dataMean(col), markerVPos, 'o', 'MarkerSize',15, 'MarkerEdgeColor',interveningLineColour, 'MarkerFaceColor','None');
      else
        p0 = errorbar(dataMean(col), markerVPos, abs(dataCI95(1,col)), abs(dataCI95(2,col)), 'horizontal', 'Marker',markerStyle,...
          'CapSize',0, 'MarkerSize',15, 'MarkerEdgeColor',interveningLineColour, 'MarkerFaceColor','None', 'Color',interveningLineColour);
      end
      uistack(p0,'bottom');
    else
      if isempty(dataCI95)
        plot(dataMean(col), markerVPos, 'o', 'MarkerSize',15, 'MarkerEdgeColor',lineColour, 'MarkerFaceColor','None')
      else
        errorbar(dataMean(col), markerVPos, abs(dataCI95(1,col)), abs(dataCI95(2,col)), 'horizontal', 'Marker',markerStyle,...
          'CapSize',0, 'MarkerSize',15, 'MarkerEdgeColor',lineColour, 'MarkerFaceColor','None', 'Color',lineColour);
      end
    end
  end
end

% Variance markers
if options.displayVars
  vars = zeros(1,nCol);
  for col = 1:nCol
    if isempty(options.colours)
      lineColour = matlabColours(col);
    else
      lineColour = options.colours(col,:);
    end
    vars(col) = var(data(:,col), 'omitnan');
    plot([dataMean(col)-vars(col)/2 dataMean(col)+vars(col)/2],...
      [0.2*max(max(bars)) 0.2*max(max(bars))]-0.02*max(max(bars))*(col-1),...
      'Color',lineColour, 'LineWidth',2, 'LineStyle','-')
  end
end
hold off

% Axes limits
if isempty(options.xLim)
  xLim = xlim;
else
  xLim = options.xLim;
end
if isempty(options.yLim)
  yLim = ylim;
else
  yLim = options.yLim;
end
ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 25, 4/3, 2, [0.005 0], 'out',...
  'on', 'k', {options.xLabel}, xLim, xticks,...
  'on', 'k', {options.yLabel}, yLim, yticks);
set(h, 'Name',options.figName);


%% Add a legend
if ~isempty(options.legendLabels)
  figureObjs = h.Children;
  legendExists = false;
  for iObj = 1:numel(figureObjs)
    if strcmpi(figureObjs(iObj).Type,'legend')
      legendExists = true;
      break
    end
  end
  if legendExists
    l = legend;
    entries2remove = [];
    for iPlot = nLegendEntries+1:numel(l.PlotChildren)
      tf = false;
      for iPlot2 = 1:numel(p)
        tf = eq(l.PlotChildren(iPlot), p(iPlot2));
        if tf
          break
        end
      end
      if tf
        continue
      else
        entries2remove = [entries2remove iPlot];
      end
    end
    for iPlot2 = 1:numel(p)
      l.String{nLegendEntries+iPlot2} = options.legendLabels{iPlot2};
    end
    entries2keep = true(1,numel(l.PlotChildren));
    entries2keep(entries2remove) = false;
    legend(l.PlotChildren(entries2keep), l.String(entries2keep), 'Location','NorthWest')
    legend boxoff
  else
    legend(p, options.legendLabels, 'Location','NorthWest')
    legend boxoff
  end
end


%% Display stats
if ~isempty(options.stats)
  % Var test
  stats = struct2table(options.stats, 'AsArray',true);
  varMat = [stats.var1 stats.var2];
  if strcmpi(options.plotType, 'extremes')
    iComp = stats.iCol1 == 1 & stats.iCol2 == 2;
  elseif strcmpi(options.plotType, 'regular')
    iComp = min(varMat,[],2) == min(min(varMat)) & max(varMat,[],2) == max(max(varMat));
    if sum(iComp) > 1
      iComp0 = find(iComp,1);
      iComp = false(size(iComp));
      iComp(iComp0) = true;
    end
  end
  if isempty(options.legendLabels)
    dist1 = 'dist1';
    dist2 = 'dist2';
  else
    dist1 = options.legendLabels{stats.iCol1(iComp)};
  	dist2 = options.legendLabels{stats.iCol2(iComp)};
  end
  var1 = round(stats.var1(iComp),2);
  var2 = round(stats.var2(iComp),2);
  pval = stats.p(iComp);
  
  xAxisLength = xLim(2)-xLim(1);
  yAxisLength = yLim(2)-yLim(1);
  textStr = [dist1 ' (\sigma^2=' num2str(var1) ') vs ' dist2 ' (\sigma^2=' num2str(var2) ') p=' num2str(pval)];
  text(xLim(2)-xAxisLength*0.85, yLim(2)-yAxisLength*0.01, textStr, 'FontSize',12);
  
  % Means t-test
  if strcmpi(options.plotType, 'extremes')
    minCol = 1;
    maxCol = nCol;
  elseif strcmpi(options.plotType, 'regular')
    [~, minCol] = min(dataMean);
    [~, maxCol] = max(dataMean);
  end
  [~, pval, ~, stats] = ttest(data(:,minCol), data(:,maxCol));
  stats.p = pval;
  if isempty(options.legendLabels)
    dist1 = 'dist1';
    dist2 = 'dist2';
  else
    if strcmpi(options.plotType, 'extremes')
      dist1 = options.legendLabels{1};
      dist2 = options.legendLabels{end};
    elseif strcmpi(options.plotType, 'regular')
      dist1 = options.legendLabels{minCol};
      dist2 = options.legendLabels{maxCol};
    end
  end
  mean1 = dataMean(minCol);
  mean2 = dataMean(maxCol);
  
  textStr = [dist1 ' (\mu=' num2str(mean1) ') vs ' dist2 ' (\mu=' num2str(mean2) ') p=' num2str(pval)];
  text(xLim(2)-xAxisLength*0.85, yLim(2)+yAxisLength*0.025, textStr, 'FontSize',12);
end


%% Save the figure
label = [3.5 3];
if ~isempty(options.stats)
  margin = [0.55 0.75];
else
  margin = [0.55 0.55];
end
width = options.figSize-label(1)-margin(1);
height = options.figSize-label(2)-margin(2);
paperSize = resizeFig(h, ax1, width, height, label, margin, 0);
if options.saveFig
  hgsave(h, [options.figFolder filesep options.figName '.fig']);
  exportFig(h, [options.figFolder filesep options.figName '.png'],'-dpng','-r300', paperSize);
  close(h);
end
end